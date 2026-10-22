# mitogenome-pipeline

[![CI](../../actions/workflows/ci.yaml/badge.svg)](../../actions/workflows/ci.yaml)

A reproducible [Snakemake](https://snakemake.readthedocs.io) pipeline that turns **public Illumina reads into a mitochondrial phylogeny**. Give it a table of species with SRA run accessions (and, optionally, published GenBank mitogenomes as extra taxa or outgroups). It downloads the reads, assembles and annotates each mitogenome, extracts the 13 protein-coding genes, aligns them, and builds a partitioned maximum-likelihood tree. Everything ends up in a one-page HTML report.

```
SRA reads ─► fastq-dump ─► Trimmomatic ─► FastQC/MultiQC
                               │
                               ▼
                MitoFinder (MEGAHIT + annotation)      GenBank mitogenomes
                               │                               │
                               └──────► 13 protein-coding genes ◄┘
                                               │
                                MAFFT ─► trimAl ─► supermatrix + partitions
                                               │
                              IQ-TREE 2 (ModelFinder per gene, UFBoot)
                                               │
                                     results/report.html
```

## Features

- **Laptop-sized:** streams only the first *N* read pairs of each run (`max_read_pairs`, default 4 million). A mitogenome is high-copy, so that's usually enough to assemble it, with no full-run download.
- **Mixed taxa:** assembled samples and published GenBank mitogenomes go in one sample sheet, so outgroups and reference species are just extra rows.
- **Fault-tolerant:** a sample where nothing assembles appears as 0/13 genes in the report instead of stopping the run. Taxa missing a gene get gaps in the supermatrix.
- **Works around the bioconda MitoFinder package:** it ships without its config file and without MiTFi, so the pipeline writes the config itself and annotates tRNAs with ARWEN. MitoFinder also breaks on paths containing spaces, so it runs in a temporary folder.
- **Gene-name normalisation:** COI/CO1/COX1, NAD4L/ND4L, `cob`/CYTB and NCBI product descriptions all map to one name.
- **Reproducible:** every step runs in its own pinned conda environment. The config and sample sheet are validated against JSON schemas. CI runs unit tests, `snakemake --lint` and a full dry run on every push.

## Quick start (Linux or WSL)

Bioconda tools run on Linux and macOS. On Windows, use WSL (Ubuntu).

```bash
# 1. Install Miniforge (once): https://github.com/conda-forge/miniforge
# 2. Create the runner environment
conda env create -f environment.yaml
conda activate mitopipe
conda config --set channel_priority strict

# 3. Put your e-mail in config/config.yaml (NCBI asks for one), then:
snakemake --sdm conda --cores 8
```

The first run installs each tool's environment, which takes a few minutes. The bundled example (`config/samples.tsv`) builds a great-ape tree from six published mitogenomes rooted on macaque, so you can check the setup before assembling anything. Open `results/report.html` when it finishes.

To test assembly too, run the second example. It adds a white-cheeked gibbon assembled from raw reads (SRR1425539, 4 million read pairs, about 800 MB). The gibbon assembles as a single 16.5 kb contig with all 13 genes, and it lands as the sister group of the great apes:

```bash
snakemake --sdm conda --cores 8 --config samples=config/samples_sra_example.tsv
```

**On WSL with the project in a Windows folder (`/mnt/c/...`):** keep the conda environments on the Linux side, where they install much faster and symlinks work. Add `--conda-prefix ~/.snakemake-conda` to the commands above.

## Your own taxa

Edit `config/samples.tsv` (tab-separated):

| sample | source | accession | species |
|---|---|---|---|
| `Species_one` | `sra` | `SRR…` (paired-end Illumina run) | Species one |
| `Outgroup` | `genbank` | `NC_…` | Outgroup species |

To use your own sequencing data, set `source` to `local`, leave `accession` empty and add two more columns, `r1` and `r2`, with the paths to the gzipped paired FASTQ files.

Then in `config/config.yaml`, set:

- `genetic_code`: 2 for vertebrates, 5 for most invertebrates (NCBI translation tables).
- `reference_accessions`: one or more published mitogenomes close to your group. MitoFinder uses them to find and name genes. Or set `reference_file` to a GenBank file on disk, such as your lab's own draft.
- `iqtree.outgroup`: the sample name the tree is rooted on.

You need at least 4 taxa with genes for a tree.

## Outputs

| Path | Contents |
|---|---|
| `results/report.html` | Tree with bootstrap support, assembly table, gene presence matrix |
| `results/qc/multiqc_report.html` | Read quality after trimming |
| `results/mitofinder/<sample>/` | Assembled mitochondrial contigs and annotated genes |
| `results/assembly_stats.tsv` | Contigs, length, GC% and genes found per taxon |
| `results/alignments/` | Per-gene protein alignments (raw and trimmed) |
| `results/phylo/` | Supermatrix, partition file, IQ-TREE tree and log |

## Development

```bash
pytest                                      # unit tests for workflow/scripts/mitolib.py
snakemake --lint --configfile .test/config.yaml
snakemake -n --configfile .test/config.yaml # dry run: builds the job graph only
snakemake --dag | dot -Tsvg > dag.svg       # picture of the workflow
```

## Tools

sra-tools, Trimmomatic, FastQC, MultiQC, MitoFinder 1.4.1 (MEGAHIT, BLAST, ARWEN), MAFFT, trimAl, IQ-TREE 2, Biopython. Please cite them if you use this pipeline.

## License

MIT
