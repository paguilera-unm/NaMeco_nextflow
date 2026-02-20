// ============================================================
//  modules/clustering.nf
//  Processes: KMER_COUNTER, CLUSTERING_UMAP, SHARED_CLUSTERS,
//             FILE_BY_CLUSTER
// ============================================================

// ── KMER_COUNTER ─────────────────────────────────────────────
process KMER_COUNTER {
    tag "${sample_id}"
    label 'high_cpu'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("kmers.tsv"), emit: kmers

    script:
    """
    #!/usr/bin/env python3
    import re, gzip, multiprocessing as mp
    from itertools import product
    from Bio import SeqIO

    L       = ${params.kmer}
    T       = ${params.threads}
    fasta   = "${fasta}"
    outfile = "kmers.tsv"

    # Build all possible kmers
    kmers = [''.join(c) for c in product('ACGT', repeat=L)]

    def kmer_subcounter(kmers, rec, q):
        counts = [str(len(re.findall(f'(?={mer})', str(rec.seq)))) for mer in kmers]
        q.put('\\t'.join([rec.id] + counts))

    def kmer_writer(q, outfile, kmers):
        with open(outfile, 'wt') as tab:
            tab.write('\\t'.join(['ID'] + kmers) + '\\n')
            while True:
                m = q.get()
                if m == 'kill':
                    break
                tab.write(m + '\\n')
                tab.flush()

    opener = gzip.open if fasta.endswith('.gz') else open
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(T)
    watcher = pool.apply_async(kmer_writer, (q, outfile, kmers))
    jobs = []
    with opener(fasta, 'rt') as fh:
        for rec in SeqIO.parse(fh, 'fasta'):
            job = pool.apply_async(kmer_subcounter, (kmers, rec, q))
            jobs.append(job)
    [job.get() for job in jobs]
    q.put('kill')
    pool.close()
    pool.join()
    print(f"K-mer counting done for ${sample_id}")
    """
}

// ── CLUSTERING_UMAP ──────────────────────────────────────────
// Receives ALL samples' kmer files at once → one execution
process CLUSTERING_UMAP {
    label 'high_mem'

    input:
    val   sample_ids           // list of sample names
    path  kmer_files           // all kmers.tsv files (flat, collected)

    output:
    path "*/clusters.tsv",           emit: per_sample_dir
    path "*/subsampled_ids.tsv",     emit: subsampled_ids
    path "shared_clusters_raw.tsv",  emit: shared_clusters

    script:
    def samples_list = sample_ids instanceof List ? sample_ids.join(',') : sample_ids
    """
    #!/usr/bin/env python3
    import os, re
    import numpy as np
    import pandas as pd
    import sklearn.cluster as cluster
    import umap

    os.environ["MKL_NUM_THREADS"]    = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"]    = "1"
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'

    SAMPLES    = "${samples_list}".split(',')
    T          = ${params.threads}
    EPS        = ${params.select_epsilon}
    CLUST_UQ   = ${params.cluster_size}
    RSTAT      = ${params.random_state}

    # Map sample → kmer file (files are staged flat in the work dir)
    kmer_files = "${kmer_files}".split()
    # Assign files to samples by order (Nextflow stages them in order)
    sample_to_file = dict(zip(SAMPLES, kmer_files))

    # ── Per-sample UMAP + HDBscan clustering ──────────────────
    dfs_sub = []
    for sample in SAMPLES:
        os.makedirs(sample, exist_ok=True)
        data = pd.read_csv(sample_to_file[sample], sep='\\t', index_col=0)
        size = max([CLUST_UQ, 10])

        emb    = umap.UMAP(n_jobs=T, metric='braycurtis', min_dist=0,
                           n_components=10).fit_transform(data.values)
        labels = cluster.HDBSCAN(min_cluster_size=size, n_jobs=T,
                                 cluster_selection_epsilon=EPS).fit_predict(emb)

        clusters_df = pd.DataFrame({'Feature': data.index, 'Cluster': labels})
        clusters_df = clusters_df.loc[clusters_df.Cluster >= 0]
        clusters_df.Cluster = 'Cluster_' + clusters_df.Cluster.astype(str)

        for cid in clusters_df.Cluster.unique():
            sub = clusters_df.loc[clusters_df.Cluster == cid].copy()
            if len(sub) > 100:
                sub = sub.sample(n=100, random_state=RSTAT)
            data.loc[sub.Feature.tolist(), 'FullID'] = sample + '___' + cid + '___'

        data = data[data['FullID'].notna()].copy()
        data.FullID = data.FullID + data.index.astype(str)
        data.set_index('FullID', inplace=True)
        data.to_csv(f'{sample}/subsampled_ids.tsv', sep='\\t')
        clusters_df.to_csv(f'{sample}/clusters.tsv', sep='\\t', index=False)
        dfs_sub.append(data)
        print(f"Sample {sample}: {len(clusters_df)} reads → "
              f"{len(clusters_df.Cluster.unique())} clusters")

    # ── Cross-sample UMAP + HDBscan clustering ────────────────
    dfs_sub = [df.apply(pd.to_numeric, downcast='integer') for df in dfs_sub]
    merged  = pd.concat(dfs_sub)

    emb    = umap.UMAP(n_jobs=T, metric='braycurtis', min_dist=0,
                       n_components=10).fit_transform(merged.values)
    labels = cluster.HDBSCAN(min_cluster_size=8, n_jobs=T,
                              cluster_selection_epsilon=EPS).fit_predict(emb)

    shared = pd.DataFrame({'Feature': merged.index, 'Cluster': labels})
    shared = shared.loc[shared.Cluster >= 0]
    shared.Cluster = 'Cluster_' + shared.Cluster.astype(str)
    shared.to_csv('shared_clusters_raw.tsv', sep='\\t', index=False)
    print(f"Cross-sample clustering: {len(shared)} features → "
          f"{len(shared.Cluster.unique())} shared clusters")
    """
}

// ── SHARED_CLUSTERS ───────────────────────────────────────────
process SHARED_CLUSTERS {
    label 'high_cpu'

    input:
    path shared_clusters_raw     // shared_clusters_raw.tsv
    path per_sample_dirs         // collected: */clusters.tsv + */subsampled_ids.tsv
    val  sample_ids

    output:
    path "cluster_counts.tsv",        emit: counts
    path "Clusters_subsampled",       emit: cluster_dir
    path "shared_clusters_final.tsv", emit: shared_final

    script:
    def samples_list = sample_ids instanceof List ? sample_ids.join(',') : sample_ids
    """
    #!/usr/bin/env python3
    import os, random
    import pandas as pd

    SAMPLES  = "${samples_list}".split(',')
    RSTAT    = ${params.random_state}
    SUBS     = ${params.subsample}

    shared   = pd.read_csv("${shared_clusters_raw}", sep='\\t', index_col=0)
    counts   = pd.DataFrame(columns=SAMPLES, index=shared.Cluster.unique()).astype(float).fillna(0)
    clust_dict = {c: [] for c in shared.Cluster.unique()}
    i = len(clust_dict) - 1

    for sample in SAMPLES:
        # per-sample files are staged as <sample>/clusters.tsv etc.
        unique    = pd.read_csv(f'{sample}/clusters.tsv',       sep='\\t', index_col=0)
        subunique = pd.read_csv(f'{sample}/subsampled_ids.tsv', sep='\\t', index_col=0)

        for uclust in unique.Cluster.unique():
            uniq    = unique.loc[unique.Cluster == uclust]
            subsize = len(subunique.loc[subunique.index == uclust])
            shar    = shared.loc[shared.index.str.contains(f"{sample}___{uclust}___")]

            if len(shar) > 0:
                shar  = shar.groupby('Cluster').size().reset_index(name='counts')
                shar  = shar.sort_values('counts', ascending=False).reset_index()
                if shar.loc[0, 'counts'] >= subsize / 2:
                    clust_dict[shar.loc[0, 'Cluster']] += uniq.index.tolist()
                    counts.loc[shar.loc[0, 'Cluster'], sample] += len(uniq)
                    continue
            i += 1
            counts.loc[f'Cluster_{i}', sample] = len(uniq)
            clust_dict[f'Cluster_{i}'] = uniq.index.tolist()

    counts = counts.astype(float).fillna(0)
    counts = counts.loc[counts.sum(axis=1) >= 10]
    counts.index.name = 'Cluster'
    counts.sort_index(key=lambda x: x.to_series().str[8:].astype(int), inplace=True)
    counts.to_csv('cluster_counts.tsv', sep='\\t')

    # ── Write per-cluster read ID files ──────────────────────
    random.seed(RSTAT)
    os.makedirs('Clusters_subsampled', exist_ok=True)
    clusters = counts.index.tolist()

    with open('Clusters_subsampled/Pooled.txt', 'w') as pooled:
        for k in clusters:
            v = clust_dict.get(k, [])
            pooled.write('{}\\n{}\\n'.format(k, '\\n'.join(v)))
            sub = random.sample(v, SUBS) if len(v) > SUBS else v
            with open(f'Clusters_subsampled/{k}.txt', 'w') as cf:
                cf.write('\\n'.join(sub))

    # Save final shared table for reference
    shared.to_csv('shared_clusters_final.tsv', sep='\\t')
    print(f"Pooled into {len(clusters)} shared clusters.")
    """
}

// ── FILE_BY_CLUSTER ───────────────────────────────────────────
process FILE_BY_CLUSTER {
    label 'high_cpu'

    input:
    path fastas           // all sample fasta files
    path cluster_dir      // Clusters_subsampled/ directory
    path counts           // cluster_counts.tsv

    output:
    path "*.fa",           emit: cluster_files
    path "*_consensus.fa", emit: consensus_files

    script:
    """
    #!/usr/bin/env python3
    import os, subprocess, multiprocessing as mp
    import pandas as pd
    from Bio import SeqIO

    T        = ${params.threads}
    SUBS     = ${params.subsample}
    clust_dir = "${cluster_dir}"
    counts_f  = "${counts}"

    # Pool all fastas
    df       = pd.read_csv(counts_f, sep='\\t', index_col=0)
    clusters = df.index.tolist()

    pooled_fa = "pooled.fa"
    fasta_files = "${fastas}".split()
    with open(pooled_fa, 'w') as out:
        for f in fasta_files:
            for rec in SeqIO.parse(f, 'fasta'):
                out.write(f'>{rec.id}\\n{str(rec.seq)}\\n')

    def split_and_consensus(cluster_id):
        txt  = f'{clust_dir}/{cluster_id}.txt'
        fa   = f'{cluster_id}.fa'
        cons = f'{cluster_id}_consensus.fa'

        result = subprocess.check_output(
            f"grep -f {txt} -F -A 3 {pooled_fa} | grep -v '^--\$'",
            shell=True).decode()
        with open(fa, 'w') as fh:
            fh.write(result)

        consensus = subprocess.check_output(f'spoa {fa}', shell=True).decode()
        seq = consensus.split('\\n')[1]
        with open(cons, 'w') as fh:
            fh.write(f'>{cluster_id}\\n{seq}\\n')

    with mp.Pool(T) as pool:
        pool.map(split_and_consensus, clusters)

    print(f"Created {len(clusters)} cluster fastas and consensus files")
    """
}
