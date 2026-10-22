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

        # MitoFinder builds unquoted shell commands from absolute paths, so it must
        # run somewhere without spaces; a Linux temp dir is also faster under WSL.
        WORK=$(mktemp -d "${{TMPDIR:-/tmp}}/mitofinder_{wildcards.sample}_XXXXXX")
        trap 'rm -rf "$WORK"' EXIT
        cd "$WORK"
        cp "$R1" reads_1.fastq.gz; cp "$R2" reads_2.fastq.gz; cp "$REF" reference.gb

        # The bioconda build ships without its installer's marker file and config,
        # and without MiTFi: point it at the tools in this env and use ARWEN for tRNAs.
        # Its helper scripts ignore -c and always read bin/Mitofinder.config.
        [ -e "$CONDA_PREFIX/bin/install.sh.ok" ] || touch "$CONDA_PREFIX/bin/install.sh.ok"
        for tool in megahit blast idba metaspades arwen trnascan mitfi; do
            echo "${{tool}}folder = $CONDA_PREFIX/bin/"
        done > "$CONDA_PREFIX/bin/Mitofinder.config"

        mitofinder --override -t arwen -j {wildcards.sample} \
            -1 reads_1.fastq.gz -2 reads_2.fastq.gz -r reference.gb \
            -o {params.code} --{params.assembler} -p {threads} -m {params.mem} > "$LOG" 2>&1
        # MitoFinder's own error logs live in the temp dir; keep them with ours.
        find {wildcards.sample} -name '*error*.log' -size +0 -exec sh -c 'echo "== $1"; cat "$1"' _ {{}} \\; >> "$LOG"

        FINAL=$(ls -d {wildcards.sample}/{wildcards.sample}_MitoFinder_*_Final_Results 2>/dev/null | head -n 1)
        if [ -n "$FINAL" ]; then rm -rf "$OUT/final_results"; cp -r "$FINAL" "$OUT/final_results"; fi
        if [ -n "$FINAL" ] && [ ! -s "$FINAL/{wildcards.sample}_final_genes_AA.fasta" ]; then
            echo "WARNING: {wildcards.sample}: a contig was assembled but no genes were annotated (see above)" | tee -a "$LOG" >&2
        fi
        copy_or_empty() {{ if [ -n "$FINAL" ] && [ -s "$FINAL/$1" ]; then cp "$FINAL/$1" "$2"; else : > "$2"; fi; }}
        copy_or_empty {wildcards.sample}_final_genes_NT.fasta "$OUT/{wildcards.sample}_genes_NT.fasta"
        copy_or_empty {wildcards.sample}_final_genes_AA.fasta "$OUT/{wildcards.sample}_genes_AA.fasta"
        # Contig files only (_mtDNA_contig.fasta or _mtDNA_contig_N.fasta), not the
        # per-gene _mtDNA_contig_genes_*.fasta files beside them.
        CONTIGS=$([ -n "$FINAL" ] && find "$FINAL" -maxdepth 1 -regextype posix-extended \
            -regex ".*/{wildcards.sample}_mtDNA_contig(_[0-9]+)?\\.fasta" | sort || true)
        if [ -n "$CONTIGS" ]; then
            echo "$CONTIGS" | xargs -d '\\n' cat > "$OUT/{wildcards.sample}_mtDNA.fasta"
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
        contigs=[f"results/mitofinder/{s}/{s}_mtDNA.fasta" for s in ASSEMBLED]
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
