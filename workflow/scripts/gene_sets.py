"""Regroup per-taxon protein files into one FASTA per gene (headers = taxon)."""

import sys
from pathlib import Path


def main(snakemake):
    sys.path.insert(0, snakemake.scriptdir)
    from mitolib import genes_by_taxon, read_fasta, write_fasta

    sys.stderr = open(snakemake.log[0], "w")
    per_gene = {g: [] for g in snakemake.params.genes}
    for taxon, path in zip(snakemake.params.taxa, snakemake.input):
        found = genes_by_taxon(read_fasta(path))
        print(f"{taxon}: {len(found)} genes", file=sys.stderr)
        for gene, seq in found.items():
            per_gene[gene].append((taxon, seq))
    outdir = Path(snakemake.params.outdir)
    for gene, records in per_gene.items():
        write_fasta(outdir / f"{gene}.faa", records)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
