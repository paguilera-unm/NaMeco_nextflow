// ============================================================
//  modules/read_correction.nf
//  Processes: READ_CORRECTION, COLLECT_POLISHED
// ============================================================

// ── READ_CORRECTION ───────────────────────────────────────────
// One process execution per cluster (fully parallelised)
process READ_CORRECTION {
    tag "${cluster_id}"
    label 'high_cpu'

    input:
    // [ cluster_id, consensus_fa, reads_fa ]
    tuple val(cluster_id), path(consensus_fa), path(reads_fa)

    output:
    tuple val(cluster_id), path("${cluster_id}_polished.fa"), emit: polished

    script:
    def N = params.n_polish
    """
    # Iterative polishing with minimap2 + racon
    current="${consensus_fa}"

    for ROUND in \$(seq 0 ${N - 1}); do
        next="${cluster_id}_racon\${ROUND}.fa"
        minimap2 -ax map-ont -t ${params.threads} \$current ${reads_fa} -o mapping.sam 2>> racon.log
        racon -m 8 -x -6 -g -8 -t ${params.threads} \\
              ${reads_fa} mapping.sam \$current > \$next 2>> racon.log
        rm -f mapping.sam
        current=\$next
    done

    cp \$current ${cluster_id}_polished.fa
    """
}

// ── COLLECT_POLISHED ──────────────────────────────────────────
// Gather all per-cluster polished fastas → single rep_seqs.fasta
process COLLECT_POLISHED {
    label 'low_cpu'

    input:
    path polished_files   // collected list: [ [cluster_id, polished_fa], ... ]
    path counts           // cluster_counts.tsv (to enforce order)

    output:
    path "rep_seqs.fasta", emit: rep_seqs

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd, os

    counts   = pd.read_csv("${counts}", sep='\\t', index_col=0)
    clusters = counts.index.tolist()

    # polished_files are staged flat in work dir as <cluster_id>_polished.fa
    with open('rep_seqs.fasta', 'w') as out:
        for cid in clusters:
            fname = f'{cid}_polished.fa'
            if os.path.exists(fname):
                with open(fname) as fh:
                    lines = fh.read().strip().split('\\n')
                    # Rename header to cluster ID for clarity
                    out.write(f'>{cid}\\n{lines[1]}\\n')
            else:
                print(f'WARNING: {fname} not found, skipping.')

    n = sum(1 for l in open('rep_seqs.fasta') if l.startswith('>'))
    print(f"rep_seqs.fasta contains {n} representative sequences")
    """
}
