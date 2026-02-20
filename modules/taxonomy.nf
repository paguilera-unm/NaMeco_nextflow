// ============================================================
//  modules/taxonomy.nf
//  Processes: TAXONOMY_DB, TAXONOMY_ANNOTATE, COLLAPSE_TAXA
// ============================================================

// ── TAXONOMY_DB ───────────────────────────────────────────────
// Download & build GTDB BLAST database (runs once, cached)
process TAXONOMY_DB {
    label 'high_cpu'
    storeDir "${params.db_path ?: "${params.out_dir}/GTDB-${params.db_version}/${params.db_type}"}"

    output:
    path ".", emit: db_dir

    script:
    """
    # Download GTDB data via RESCRIPt
    qiime rescript get-gtdb-data \\
        --p-db-type '${params.db_type}' \\
        --p-version ${params.db_version} \\
        --o-gtdb-taxonomy taxa.qza \\
        --o-gtdb-sequences seqs.qza

    qiime tools export --input-path seqs.qza  --output-path .
    qiime tools export --input-path taxa.qza  --output-path .

    # Build BLAST database (without the raw fasta to save space)
    makeblastdb -in dna-sequences.fasta -parse_seqids -dbtype nucl
    rm -f dna-sequences.fasta    # kept as .ndb/.nin/... files by makeblastdb
    """
}

// ── TAXONOMY_ANNOTATE ──────────────────────────────────────────
process TAXONOMY_ANNOTATE {
    label 'high_cpu'

    input:
    path rep_seqs
    path db_dir
    path counts

    output:
    path "Taxonomy.tsv",    emit: taxonomy
    path "Taxonomy-q2.tsv", emit: taxonomy_q2, optional: true

    script:
    def mask_flag = params.mask_taxa ? "True" : "False"
    """
    #!/usr/bin/env python3
    import os, subprocess
    import pandas as pd

    Q       = "${rep_seqs}"
    DB      = "${db_dir}"
    GAP     = ${params.gap}
    FRAC    = ${params.min_fraction}
    MASK    = ${mask_flag} == 'True'
    THREADS = ${params.threads}

    THRESHOLDS = {
        'Domain': 65, 'Phylum': 75, 'Class': 78.5,
        'Order': 82, 'Family': 86.5, 'Genus': 94.5, 'Species': 97
    }

    queries = [l[1:].split()[0].rstrip()
               for l in open(Q) if l.startswith('>')]

    # ── BLAST ─────────────────────────────────────────────────
    blast_out = 'blastn.tsv'
    subprocess.run(
        f'blastn -query {Q} -db {DB}/dna-sequences.fasta '
        f'-task blastn -qcov_hsp_perc 80 '
        f'-num_threads {THREADS} -out {blast_out} '
        f'-max_target_seqs 50 -max_hsps 50 '
        f'-outfmt "6 qseqid sseqid evalue length pident nident bitscore score gaps"',
        shell=True, check=True)

    blast = pd.read_csv(blast_out, sep='\\t', header=None,
            names=['Cluster','SeqID','eval','length','pind',
                   'nind','bitscore','score','gaps'])
    blast = blast.sort_values(['bitscore','eval'], ascending=[False,False])

    # ── Load taxonomy mapping ─────────────────────────────────
    mapp = pd.read_csv(f'{DB}/taxonomy.tsv', sep='\\t')
    mapp.Taxon = mapp.Taxon.apply(
        lambda x: x.rsplit(';',1)[0] + ';' +
                  ' '.join(x.rsplit(';',1)[-1].split(' ')[:2])
                    .replace('_',' ').replace('  ','__'))
    mapping = dict(mapp[['Feature ID','Taxon']].values)

    # ── Helper: apply percent-identity masking ────────────────
    def taxonomy_thresholds(bclust):
        for ind in bclust.index:
            taxon = bclust.loc[ind, 'Taxon']
            last  = ''
            for rank, perc in THRESHOLDS.items():
                prefix = f"{rank[0].lower()}__"
                pat    = taxon.split(prefix)[-1].split(';')[0]
                if bclust.loc[ind,'pind'] >= perc:
                    last = pat
                if bclust.loc[ind,'pind'] < perc:
                    taxon = taxon.replace(prefix+pat, f"{prefix}{last} unclassified")
                    bclust.loc[ind,'Taxon'] = taxon
        return bclust

    # ── Helper: select top-hit taxonomy ──────────────────────
    def top_hit(bclust):
        if len(bclust) == 0:
            return ';'.join([f'{r}__Unclassified' for r in 'dpcofgs']), 0
        tc = bclust['Taxon'].value_counts()
        bclust = bclust.copy()
        bclust['Taxa_counts'] = bclust['Taxon'].map(tc)
        bclust.sort_values(['Taxa_counts','bitscore','pind'],
                           ascending=[False,False,False], inplace=True)
        taxon = bclust.Taxon.iloc[0]
        pind  = bclust.pind.iloc[0]
        if (len(bclust.loc[bclust.Taxa_counts == bclust.Taxa_counts.max()])
                / len(bclust) < FRAC):
            taxon = (taxon.rsplit(';',1)[0] + ';' +
                     taxon.rsplit(';',1)[-1].split(' ')[0] + ' unclassified')
        return taxon, pind

    # ── Annotate each cluster ─────────────────────────────────
    taxa = pd.DataFrame(columns=['Taxon','Perc. id.'])
    for cluster in queries:
        bclust = blast.loc[blast.Cluster == cluster].copy()
        bclust = bclust.loc[bclust.bitscore > bclust.bitscore.max() - GAP]
        bclust['Taxon'] = bclust['SeqID'].map(mapping)
        if MASK:
            bclust = taxonomy_thresholds(bclust)
        taxa.loc[cluster, ['Taxon','Perc. id.']] = top_hit(bclust)

    for cluster in queries:
        if cluster not in blast.Cluster.tolist():
            taxa.loc[cluster, 'Taxon'] = 'Unclassified'

    # ── Save Q2-style taxonomy ────────────────────────────────
    taxa.index.name = 'Feature ID'
    taxa.to_csv('Taxonomy-q2.tsv', sep='\\t')

    # ── Expand to per-rank columns ────────────────────────────
    taxa_exp = taxa.copy()
    taxa_exp.index.name = 'Cluster'
    for rank in THRESHOLDS:
        taxa_exp[rank] = taxa_exp.Taxon.apply(
            lambda x: x.split(f"{rank[0].lower()}__")[-1].split(';')[0])
    taxa_exp.drop('Taxon', axis=1, inplace=True)
    taxa_exp.to_csv('Taxonomy.tsv', sep='\\t')
    print("Taxonomy annotation complete.")
    """
}

// ── COLLAPSE_TAXA ─────────────────────────────────────────────
// Collapse abundance counts to each taxonomic rank
process COLLAPSE_TAXA {
    label 'low_cpu'

    input:
    path taxonomy   // Taxonomy.tsv
    path counts     // cluster_counts.tsv

    output:
    path "*_counts.tsv", emit: rank_counts

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    taxa   = pd.read_csv("${taxonomy}", sep='\\t', index_col=0)
    counts = pd.read_csv("${counts}",   sep='\\t', index_col=0)

    THRESHOLDS = ['Domain','Phylum','Class','Order','Family','Genus','Species']

    for rank in THRESHOLDS:
        if rank not in taxa.columns:
            continue
        coll = counts.copy()
        coll[rank] = taxa[rank]
        coll = coll.groupby(rank).sum()
        coll.to_csv(f'{rank}_counts.tsv', sep='\\t')
        print(f'Collapsed to {rank}: {len(coll)} entries')
    """
}
