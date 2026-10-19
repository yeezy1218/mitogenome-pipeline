rule gene_sets:
    """Regroup sequences from per-taxon files into one FASTA per gene."""
    input:
        sample_genes,
    output:
        expand("results/genes/{gene}.faa", gene=GENES),
    params:
        taxa=TAXA,
        genes=GENES,
        outdir=lambda wc, output: os.path.dirname(output[0]),
    log:
        "logs/gene_sets.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/gene_sets.py"


rule align:
    input:
        "results/genes/{gene}.faa",
    output:
        "results/alignments/{gene}.aln.faa",
    log:
        "logs/align/{gene}.log",
    threads: 2
    conda:
        "../envs/phylo.yaml"
    shell:
        """
        if [ "$(grep -c '^>' {input})" -ge 2 ]; then
            mafft --auto --thread {threads} {input} > {output} 2> {log}
        else
            cp {input} {output}; echo "fewer than 2 sequences, not aligned" > {log}
        fi
        """


rule trim_alignment:
    """Drop poorly aligned columns (trimAl's automated1 heuristic)."""
    input:
        "results/alignments/{gene}.aln.faa",
    output:
        "results/alignments/{gene}.trim.faa",
    log:
        "logs/trim/{gene}.log",
    conda:
        "../envs/phylo.yaml"
    shell:
        """
        if [ "$(grep -c '^>' {input})" -ge 2 ]; then
            trimal -in {input} -out {output} -automated1 > {log} 2>&1 || : > {output}
        else
            cp {input} {output}; echo "fewer than 2 sequences, not trimmed" > {log}
        fi
        """


rule supermatrix:
    """Concatenate gene alignments; missing genes become gaps. Writes the
    gene boundaries as a NEXUS partition file for IQ-TREE."""
    input:
        expand("results/alignments/{gene}.trim.faa", gene=GENES),
    output:
        matrix="results/phylo/supermatrix.faa",
        partitions="results/phylo/partitions.nex",
    params:
        taxa=TAXA,
        genes=GENES,
    log:
        "logs/supermatrix.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/supermatrix.py"


rule iqtree:
    """Maximum-likelihood tree: best model per gene (ModelFinder), ultrafast
    bootstrap support. Needs at least 4 taxa."""
    input:
        matrix="results/phylo/supermatrix.faa",
        partitions="results/phylo/partitions.nex",
    output:
        tree="results/phylo/mito.treefile",
        summary="results/phylo/mito.iqtree",
    params:
        prefix=lambda wc, output: output.tree.removesuffix(".treefile"),
        bootstrap=config["iqtree"]["bootstrap"],
        outgroup=(f"-o {config['iqtree']['outgroup']}" if config["iqtree"]["outgroup"] else ""),
    log:
        "logs/iqtree.log",
    threads: config["threads"]["iqtree"]
    conda:
        "../envs/phylo.yaml"
    shell:
        """
        iqtree2 -s {input.matrix} -p {input.partitions} -m MFP -B {params.bootstrap} \
            {params.outgroup} -T {threads} --prefix {params.prefix} -redo > {log} 2>&1
        """
