"""Protein-coding genes and genome sequence of a published GenBank mitogenome."""

import sys

from Bio import SeqIO


def main(snakemake):
    sys.path.insert(0, snakemake.scriptdir)
    from mitolib import clean_protein, gene_from_product, normalize_gene, write_fasta

    sys.stderr = open(snakemake.log[0], "w")
    taxon = snakemake.wildcards.sample
    genes, genome = {}, []
    for rec in SeqIO.parse(snakemake.input[0], "genbank"):
        genome.append((rec.id, str(rec.seq)))
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            gene = normalize_gene(q.get("gene", [None])[0]) or gene_from_product(q.get("product", [None])[0])
            if gene is None:
                continue
            if "translation" in q:
                protein = q["translation"][0]
            else:
                table = int(q.get("transl_table", [snakemake.params.code])[0])
                protein = str(feat.extract(rec.seq).translate(table=table, to_stop=True))
            protein = clean_protein(protein)
            if len(protein) > len(genes.get(gene, "")):
                genes[gene] = protein
    write_fasta(snakemake.output.contigs, genome)
    write_fasta(snakemake.output.aa, [(f"{taxon}@{g}", s) for g, s in sorted(genes.items())])
    print(f"{taxon}: {len(genes)} protein-coding genes: {', '.join(sorted(genes))}", file=sys.stderr)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
