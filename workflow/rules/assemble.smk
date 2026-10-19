rule mitofinder:
    """Assemble with MEGAHIT and annotate the mitogenome against the reference.
    A sample where nothing is found gets empty outputs (the report shows 0 of 13
    genes) instead of stopping the whole run."""
    input:
        r1="results/reads/trimmed/{sample}_1P.fastq.gz",
        r2="results/reads/trimmed/{sample}_2P.fastq.gz",
        ref="results/reference/reference.gb",
    output:
        contigs="results/mitofinder/{sample}/{sample}_mtDNA.fasta",
        nt="results/mitofinder/{sample}/{sample}_genes_NT.fasta",
        aa="results/mitofinder/{sample}/{sample}_genes_AA.fasta",
    params:
        code=config["genetic_code"],
        assembler=config["assembler"],
        mem=config["mitofinder_memory_gb"],
    log:
        "logs/mitofinder/{sample}.log",
    threads: config["threads"]["mitofinder"]
    conda:
        "../envs/mitofinder.yaml"
    shell:
        """
        R1=$(realpath {input.r1}); R2=$(realpath {input.r2}); REF=$(realpath {input.ref})
        LOG=$(realpath -m {log}); OUT=$(realpath -m "$(dirname {output.nt})")
        mkdir -p "$OUT/run" && cd "$OUT/run"
        mitofinder --override -j {wildcards.sample} -1 "$R1" -2 "$R2" -r "$REF" \
            -o {params.code} --{params.assembler} -p {threads} -m {params.mem} > "$LOG" 2>&1

        FINAL=$(ls -d {wildcards.sample}/{wildcards.sample}_MitoFinder_*_Final_Results 2>/dev/null | head -n 1)
        copy_or_empty() {{ if [ -n "$FINAL" ] && [ -s "$FINAL/$1" ]; then cp "$FINAL/$1" "$2"; else : > "$2"; fi; }}
        copy_or_empty {wildcards.sample}_final_genes_NT.fasta "$OUT/{wildcards.sample}_genes_NT.fasta"
        copy_or_empty {wildcards.sample}_final_genes_AA.fasta "$OUT/{wildcards.sample}_genes_AA.fasta"
        if [ -n "$FINAL" ] && ls "$FINAL/{wildcards.sample}_mtDNA_contig"*.fasta > /dev/null 2>&1; then
            cat "$FINAL/{wildcards.sample}_mtDNA_contig"*.fasta > "$OUT/{wildcards.sample}_mtDNA.fasta"
        else
            : > "$OUT/{wildcards.sample}_mtDNA.fasta"
            echo "WARNING: MitoFinder found no mitochondrial contig for {wildcards.sample}" >> "$LOG"
        fi
        """


rule genbank_genes:
    """Protein-coding genes of a published mitogenome, from its CDS features."""
    input:
        "results/genbank/{sample}/{sample}.gb",
    output:
        contigs="results/genbank/{sample}/{sample}_mtDNA.fasta",
        aa="results/genbank/{sample}/{sample}_genes_AA.fasta",
    params:
        code=config["genetic_code"],
    log:
        "logs/genbank/{sample}_genes.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/genbank_genes.py"


rule assembly_stats:
    """One row per taxon: contigs, length, GC and which genes were recovered."""
    input:
        contigs=[f"results/mitofinder/{s}/{s}_mtDNA.fasta" for s in SRA]
        + [f"results/genbank/{s}/{s}_mtDNA.fasta" for s in GENBANK],
        genes=sample_genes,
    output:
        "results/assembly_stats.tsv",
    params:
        taxa=TAXA,
        sources=[SAMPLES.loc[s, "source"] for s in TAXA],
        species=[SAMPLES.loc[s, "species"] for s in TAXA],
        genes=GENES,
    log:
        "logs/assembly_stats.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assembly_stats.py"
