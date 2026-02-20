# NaMeco — Nextflow DSL2 Adaptation

Adaptación del pipeline [NaMeco](https://doi.org/10.1186/s12864-025-12415-x)
a **Nextflow DSL2** para ejecución eficiente, paralela y reproducible.

> **Versión actualizada para Nextflow ≥ 24.10**
> `publishDir` en procesos está deprecado — las salidas se gestionan centralizadamente desde `main.nf`.

---

## Estructura del proyecto

```
nameco_nextflow/
├── main.nf                  # Workflow principal + sección publish: + bloque output {}
├── nextflow.config          # Configuración, perfiles y opciones globales de publicación
├── environment.yml          # Entorno conda
└── modules/
    ├── quality_control.nf   # CHOPPER, EXTRACT_READS
    ├── clustering.nf        # KMER_COUNTER, CLUSTERING_UMAP, SHARED_CLUSTERS, FILE_BY_CLUSTER
    ├── read_correction.nf   # READ_CORRECTION (por cluster), COLLECT_POLISHED
    └── taxonomy.nf          # TAXONOMY_DB, TAXONOMY_ANNOTATE, COLLAPSE_TAXA
```

---

## Uso básico

```bash
# Instalación de Nextflow
curl -s https://get.nextflow.io | bash

# Ejecución mínima
nextflow run main.nf --inp_dir /ruta/a/reads --out_dir resultados

# Con perfil conda
nextflow run main.nf --inp_dir /ruta/a/reads -profile conda

# Con SLURM + conda
nextflow run main.nf --inp_dir /ruta/a/reads -profile slurm,conda

# Reanudar tras error
nextflow run main.nf --inp_dir /ruta/a/reads -resume

# Perfil de test
nextflow run main.nf -profile test
```

---

## Parámetros principales

| Parámetro         | Default                    | Descripción                                     |
|-------------------|----------------------------|-------------------------------------------------|
| `--inp_dir`       | *requerido*                | Carpeta con archivos `.fastq`/`.fq`/`.gz`      |
| `--out_dir`       | `NaMeco_out`               | Carpeta de salida (= `outputDir`)               |
| `--threads`       | `2`                        | CPUs por proceso                                |
| `--qc`            | `true`                     | Ejecutar Chopper (usar `--no-qc` para omitir)  |
| `--phred`         | `10`                       | Mínimo Phred para Chopper                       |
| `--min_length`    | `1200`                     | Longitud mínima de lectura                      |
| `--max_length`    | `2000`                     | Longitud máxima de lectura                      |
| `--primer_F`      | `AGAGTTTGATCMTGGCTCAG`     | Primer forward                                  |
| `--primer_R`      | `CGGTTACCTTGTTACGACTT`     | Primer reverse                                  |
| `--kmer`          | `5`                        | Longitud de k-mers                              |
| `--cluster_size`  | `10`                       | Tamaño mínimo de cluster                        |
| `--subsample`     | `200`                      | Submuestreo máximo por cluster                  |
| `--select_epsilon`| `0.1`                      | Epsilon para HDBscan                            |
| `--n_polish`      | `3`                        | Rondas de pulido con Racon                      |
| `--db_type`       | `All`                      | `All` o `SpeciesReps` (GTDB)                   |
| `--db_version`    | `226.0`                    | Versión de GTDB                                 |
| `--mask_taxa`     | `true`                     | Enmascarar rangos por % identidad               |

---

## Salidas

```
NaMeco_out/
├── Quality_control/
│   ├── Chopper/               # Lecturas filtradas por calidad
│   └── Fastas/                # Lecturas extraídas por primers
├── Clustering/
│   ├── <sample>/
│   │   ├── kmers.tsv
│   │   ├── clusters.tsv
│   │   └── subsampled_ids.tsv
│   └── Clusters_subsampled/   # FASTA + consensos por cluster
├── Read_correction/
│   └── *_polished.fa
├── Final_output/
│   ├── cluster_counts.tsv
│   ├── rep_seqs.fasta
│   ├── Taxonomy.tsv
│   ├── Taxonomy-q2.tsv
│   ├── Domain_counts.tsv
│   ├── Phylum_counts.tsv
│   └── ...
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

---

## Cita

Si usas este pipeline, cita el trabajo original:
> NaMeco: https://doi.org/10.1186/s12864-025-12415-x

Y las herramientas incluidas: Chopper, RESCRIPt, UMAP, HDBscan, SPOA, minimap2, Racon, BLAST, GTDB.

