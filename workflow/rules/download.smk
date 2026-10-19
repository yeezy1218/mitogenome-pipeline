rule download_reads:
    """Stream paired reads straight from SRA. With max_read_pairs > 0 only the
    first N spots are fetched, which is plenty for a high-copy mitogenome and
    keeps downloads laptop-sized."""
    output:
        r1="results/reads/raw/{sample}_1.fastq.gz",
        r2="results/reads/raw/{sample}_2.fastq.gz",
    params:
        acc=lambda wc: SAMPLES.loc[wc.sample, "accession"],
        limit=lambda wc: (
            f"-X {config['max_read_pairs']}" if int(config["max_read_pairs"]) > 0 else ""
        ),
        outdir=lambda wc: f"results/reads/raw/tmp_{wc.sample}",
    log:
        "logs/download/{sample}.log",
    conda:
        "../envs/sra.yaml"
    retries: 2
    shell:
        """
        (
        rm -rf {params.outdir} && mkdir -p {params.outdir}
        fastq-dump --split-3 --gzip {params.limit} -O {params.outdir} {params.acc}
        mv {params.outdir}/{params.acc}_1.fastq.gz {output.r1}
        mv {params.outdir}/{params.acc}_2.fastq.gz {output.r2}
        rm -rf {params.outdir}
        ) > {log} 2>&1
        """


rule fetch_genbank:
    """Download a published mitogenome (GenBank flat file) from NCBI."""
    output:
        "results/genbank/{sample}/{sample}.gb",
    params:
        accessions=lambda wc: [SAMPLES.loc[wc.sample, "accession"]],
        email=config["entrez_email"],
    log:
        "logs/genbank/{sample}.log",
    conda:
        "../envs/python.yaml"
    retries: 2
    script:
        "../scripts/fetch_genbank.py"


rule fetch_reference:
    """Reference mitogenome(s) MitoFinder uses to find and annotate genes."""
    output:
        "results/reference/reference.gb",
    params:
        accessions=config["reference_accessions"],
        email=config["entrez_email"],
    log:
        "logs/reference.log",
    conda:
        "../envs/python.yaml"
    retries: 2
    script:
        "../scripts/fetch_genbank.py"
