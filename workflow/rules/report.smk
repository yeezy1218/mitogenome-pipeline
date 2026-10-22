rule report:
    """Single self-contained HTML page: tree, assembly table, gene matrix."""
    input:
        stats="results/assembly_stats.tsv",
        tree="results/phylo/mito.treefile",
        summary="results/phylo/mito.iqtree",
        partitions="results/phylo/partitions.nex",
    output:
        "results/report.html",
    params:
        genes=GENES,
        outgroup=config["iqtree"]["outgroup"],
        has_qc=bool(ASSEMBLED),
    log:
        "logs/report.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/report.py"
