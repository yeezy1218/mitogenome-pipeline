"""Per-taxon table: contigs, length, GC, genes recovered."""

import sys

import pandas as pd


def main(snakemake):
    sys.path.insert(0, snakemake.scriptdir)
    from mitolib import genes_by_taxon, read_fasta, seq_stats

    sys.stderr = open(snakemake.log[0], "w")
    p = snakemake.params
    rows = []
    for taxon, source, species, contigs, genes in zip(
        p.taxa, p.sources, p.species, snakemake.input.contigs, snakemake.input.genes
    ):
        found = genes_by_taxon(read_fasta(genes))
        row = {"taxon": taxon, "species": species or taxon.replace("_", " "), "source": source}
        row.update(seq_stats(read_fasta(contigs)))
        row["genes_found"] = len(found)
        row.update({g: g in found for g in p.genes})
        rows.append(row)
    pd.DataFrame(rows).to_csv(snakemake.output[0], sep="\t", index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
