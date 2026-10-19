def sample_genes(wildcards=None):
    """Protein FASTA of every taxon, assembled or downloaded."""
    return [f"results/mitofinder/{s}/{s}_genes_AA.fasta" for s in SRA] + [
        f"results/genbank/{s}/{s}_genes_AA.fasta" for s in GENBANK
    ]
