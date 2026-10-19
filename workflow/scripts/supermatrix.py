"""Concatenate trimmed gene alignments into a supermatrix plus a partition file."""

import sys
from pathlib import Path


def main(snakemake):
    sys.path.insert(0, snakemake.scriptdir)
    from mitolib import build_supermatrix, partitions_nexus, read_fasta, write_fasta

    sys.stderr = open(snakemake.log[0], "w")
    alignments = {}
    for gene, path in zip(snakemake.params.genes, snakemake.input):
        # trimAl may append extra words to headers; the taxon is the first one.
        alignments[gene] = {h.split()[0]: s.upper() for h, s in read_fasta(path)}
    taxa = list(snakemake.params.taxa)
    matrix, parts = build_supermatrix(alignments, taxa)
    if not parts:
        raise SystemExit("No gene alignments to concatenate: did any assembly find genes?")
    present = [t for t in taxa if matrix[t].strip("-")]
    dropped = sorted(set(taxa) - set(present))
    if dropped:
        print(f"Dropping taxa with no genes: {', '.join(dropped)}", file=sys.stderr)
    if len(present) < 4:
        raise SystemExit(f"Need at least 4 taxa with genes for a tree, have {len(present)}: {present}")
    write_fasta(snakemake.output.matrix, [(t, matrix[t]) for t in present])
    Path(snakemake.output.partitions).write_text(partitions_nexus(parts))
    print(f"{len(present)} taxa x {parts[-1][2]} sites in {len(parts)} genes", file=sys.stderr)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
