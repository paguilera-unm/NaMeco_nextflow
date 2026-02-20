#!/usr/bin/env nextflow

// ============================================================
//  NaMeco — Nanopore Metabarcoding Consensus Pipeline
//  Nextflow DSL2 adaptation  (Nextflow >= 24.10)
//  Original: https://doi.org/10.1186/s12864-025-12415-x
//
//  Publishing: gestionado desde el workflow (output block)
//  publishDir en procesos está deprecado desde Nextflow 24.x
// ============================================================

log.info """
╔══════════════════════════════════════════════════════════╗
║                  N a M e C o                             ║
║    Nanopore Metabarcoding Consensus pipeline             ║
╚══════════════════════════════════════════════════════════╝
  Input dir   : ${params.inp_dir}
  Output dir  : ${params.out_dir}
  Threads     : ${params.threads}
  QC          : ${params.qc}
  DB type     : ${params.db_type}
  DB version  : ${params.db_version}
""".stripIndent()

// ============================================================
//  INCLUDE MODULES
// ============================================================
include { CHOPPER          } from './modules/quality_control.nf'
include { EXTRACT_READS    } from './modules/quality_control.nf'
include { KMER_COUNTER     } from './modules/clustering.nf'
include { CLUSTERING_UMAP  } from './modules/clustering.nf'
include { SHARED_CLUSTERS  } from './modules/clustering.nf'
include { FILE_BY_CLUSTER  } from './modules/clustering.nf'
include { READ_CORRECTION  } from './modules/read_correction.nf'
include { COLLECT_POLISHED } from './modules/read_correction.nf'
include { TAXONOMY_DB      } from './modules/taxonomy.nf'
include { TAXONOMY_ANNOTATE} from './modules/taxonomy.nf'
include { COLLAPSE_TAXA    } from './modules/taxonomy.nf'

// ============================================================
//  WORKFLOW
// ============================================================
workflow {

    main:
    // ── Input channel ──────────────────────────────────────
    ch_reads = Channel
        .fromPath("${params.inp_dir}/*.{fastq,fq,fastq.gz,fq.gz}")
        .map { file -> [ file.simpleName, file ] }

    // ── Quality Control ────────────────────────────────────
    if (params.qc) {
        CHOPPER(ch_reads)
        ch_qc = CHOPPER.out.chopped
    } else {
        log.warn "QC (Chopper) disabled — usando reads directamente"
        ch_qc = ch_reads
    }

    EXTRACT_READS(ch_qc)
    ch_fastas = EXTRACT_READS.out.fasta     // [ sample_id, fasta_file ]

    // ── Clustering ─────────────────────────────────────────
    KMER_COUNTER(ch_fastas)

    ch_all_kmers = KMER_COUNTER.out.kmers
        .map { sample, tsv -> tsv }
        .collect()

    ch_sample_ids = KMER_COUNTER.out.kmers
        .map { sample, tsv -> sample }
        .collect()

    CLUSTERING_UMAP(ch_sample_ids, ch_all_kmers)

    ch_per_sample_clusters = CLUSTERING_UMAP.out.per_sample_dir.collect()

    SHARED_CLUSTERS(
        CLUSTERING_UMAP.out.shared_clusters,
        ch_per_sample_clusters,
        ch_sample_ids
    )

    ch_all_fastas = ch_fastas
        .map { sample, fasta -> fasta }
        .collect()

    FILE_BY_CLUSTER(
        ch_all_fastas,
        SHARED_CLUSTERS.out.cluster_dir,
        SHARED_CLUSTERS.out.counts
    )

    // ── Read Correction ────────────────────────────────────
    ch_clusters = FILE_BY_CLUSTER.out.cluster_files
        .flatten()
        .filter { it.name.endsWith('_consensus.fa') }
        .map { f -> [ f.name.replaceAll('_consensus\\.fa$', ''), f ] }

    ch_cluster_reads = FILE_BY_CLUSTER.out.cluster_files
        .flatten()
        .filter { it.name =~ /^Cluster_\d+\.fa$/ && !it.name.contains('consensus') }
        .map { f -> [ f.name.replaceAll('\\.fa$', ''), f ] }

    // [ cluster_id, consensus_fa, reads_fa ]
    ch_for_correction = ch_clusters.join(ch_cluster_reads)

    READ_CORRECTION(ch_for_correction)

    COLLECT_POLISHED(
        READ_CORRECTION.out.polished.collect(),
        SHARED_CLUSTERS.out.counts
    )

    // ── Taxonomy Annotation ────────────────────────────────
    TAXONOMY_DB()

    TAXONOMY_ANNOTATE(
        COLLECT_POLISHED.out.rep_seqs,
        TAXONOMY_DB.out.db_dir,
        SHARED_CLUSTERS.out.counts
    )

    COLLAPSE_TAXA(
        TAXONOMY_ANNOTATE.out.taxonomy,
        SHARED_CLUSTERS.out.counts
    )

    // ============================================================
    //  PUBLISH SECTION — gestión de salidas centralizada aquí
    //  Reemplaza todos los publishDir dispersos en los módulos
    // ============================================================
    publish:

    // Quality Control
    CHOPPER.out.chopped                 >> 'Quality_control/Chopper'
    EXTRACT_READS.out.fasta             >> 'Quality_control/Fastas'

    // Clustering — por muestra
    KMER_COUNTER.out.kmers              >> 'Clustering'
    CLUSTERING_UMAP.out.per_sample_dir  >> 'Clustering'
    CLUSTERING_UMAP.out.subsampled_ids  >> 'Clustering'
    CLUSTERING_UMAP.out.shared_clusters >> 'Clustering'

    // Clustering — resultados globales
    SHARED_CLUSTERS.out.cluster_dir     >> 'Clustering/Clusters_subsampled'
    SHARED_CLUSTERS.out.shared_final    >> 'Clustering'
    FILE_BY_CLUSTER.out.cluster_files   >> 'Clustering/Clusters_subsampled'
    FILE_BY_CLUSTER.out.consensus_files >> 'Clustering/Clusters_subsampled'

    // Read correction
    READ_CORRECTION.out.polished        >> 'Read_correction'

    // Final outputs
    SHARED_CLUSTERS.out.counts          >> 'Final_output'
    COLLECT_POLISHED.out.rep_seqs       >> 'Final_output'
    TAXONOMY_ANNOTATE.out.taxonomy      >> 'Final_output'
    TAXONOMY_ANNOTATE.out.taxonomy_q2   >> 'Final_output'
    COLLAPSE_TAXA.out.rank_counts       >> 'Final_output'
}

// ============================================================
//  OUTPUT BLOCK — configuración global de publicación
//  Cada target corresponde a un destino del publish: de arriba
// ============================================================
output {

    // ── Quality Control ──────────────────────────────────────
    'Quality_control/Chopper' {
        mode    'copy'
    }
    'Quality_control/Fastas' {
        mode    'copy'
    }

    // ── Clustering ───────────────────────────────────────────
    'Clustering' {
        mode    'copy'
    }
    'Clustering/Clusters_subsampled' {
        mode    'copy'
    }

    // ── Read correction ──────────────────────────────────────
    'Read_correction' {
        mode    'copy'
    }

    // ── Final outputs ─────────────────────────────────────────
    'Final_output' {
        mode    'copy'
    }
}
