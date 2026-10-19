rule trimmomatic:
    input:
        r1="results/reads/raw/{sample}_1.fastq.gz",
        r2="results/reads/raw/{sample}_2.fastq.gz",
    output:
        r1="results/reads/trimmed/{sample}_1P.fastq.gz",
        r2="results/reads/trimmed/{sample}_2P.fastq.gz",
        u1=temp("results/reads/trimmed/{sample}_1U.fastq.gz"),
        u2=temp("results/reads/trimmed/{sample}_2U.fastq.gz"),
    params:
        adapters=config["trimmomatic"]["adapters"],
        settings=config["trimmomatic"]["settings"],
    log:
        "logs/trimmomatic/{sample}.log",
    threads: 4
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        """
        ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/{params.adapters}"
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} {output.u1} {output.r2} {output.u2} \
            ILLUMINACLIP:"$ADAPTERS":2:30:10 {params.settings} \
            2> {log}
        """


rule fastqc:
    input:
        "results/reads/trimmed/{sample}_{read}P.fastq.gz",
    output:
        html="results/qc/fastqc/{sample}_{read}P_fastqc.html",
        zip="results/qc/fastqc/{sample}_{read}P_fastqc.zip",
    wildcard_constraints:
        read="1|2",
    log:
        "logs/fastqc/{sample}_{read}.log",
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        "fastqc --quiet --threads {threads} -o results/qc/fastqc {input} > {log} 2>&1"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_{read}P_fastqc.zip", sample=SRA, read=["1", "2"]),
        expand("logs/trimmomatic/{sample}.log", sample=SRA),
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc --force --outdir results/qc --filename multiqc_report.html {input} > {log} 2>&1"
