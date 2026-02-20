// ============================================================
//  modules/quality_control.nf
//  Processes: CHOPPER, EXTRACT_READS
// ============================================================

// ── CHOPPER ──────────────────────────────────────────────────
process CHOPPER {
    tag "${sample_id}"
    label 'high_cpu'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.fq.gz"), emit: chopped

    script:
    def pre = reads.name.endsWith('.gz') ? "gunzip -c ${reads} -q |" : "cat ${reads} |"
    """
    ${pre} chopper \\
        -q ${params.phred} \\
        -l ${params.min_length} \\
        --maxlength ${params.max_length} \\
        -t ${params.threads} \\
        2> chopper_${sample_id}.log \\
        | gzip > ${sample_id}.fq.gz
    """
}

// ── EXTRACT_READS ─────────────────────────────────────────────
// Uses QIIME2 / RESCRIPt to extract amplicons by primers
process EXTRACT_READS {
    tag "${sample_id}"
    label 'high_cpu'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: fasta

    script:
    // Convert fastq → fasta → import to QIIME2 → extract reads → export
    def reads_input = reads.name.endsWith('.gz') ? reads : reads
    """
    # Convert FASTQ to FASTA
    python3 -c "
from Bio import SeqIO
import gzip, os

infile = '${reads}'
outfile = '${sample_id}.fa'

opener = gzip.open if infile.endswith('.gz') else open
with opener(infile, 'rt') as fh, open(outfile, 'w') as out:
    SeqIO.convert(fh, 'fastq', out, 'fasta')
print('Converted to FASTA')
"

    # Import into QIIME2
    qiime tools import \\
        --type 'FeatureData[Sequence]' \\
        --input-path ${sample_id}.fa \\
        --output-path ${sample_id}.qza

    # Extract reads by primers
    qiime feature-classifier extract-reads \\
        --i-sequences ${sample_id}.qza \\
        --p-f-primer ${params.primer_F} \\
        --p-r-primer ${params.primer_R} \\
        --p-read-orientation 'both' \\
        --p-min-length ${params.min_length} \\
        --p-max-length ${params.max_length} \\
        --p-identity ${params.primer_PI} \\
        --p-n-jobs ${params.threads} \\
        --o-reads ${sample_id}-extr.qza

    # Export to fasta
    qiime tools export \\
        --input-path ${sample_id}-extr.qza \\
        --output-path extracted_${sample_id}

    mv extracted_${sample_id}/dna-sequences.fasta ${sample_id}.fasta

    # Clean up
    rm -f ${sample_id}.fa ${sample_id}.qza ${sample_id}-extr.qza
    """
}
