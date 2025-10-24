configfile: "config/qcspades.yaml"

rule all:
    input:
        expand("output/quast/{sample}/report.html", sample=config["sample_id"])

# ------------------------------------------------------------------
# Rule 1: Run Quality Control (FastQC)
# ------------------------------------------------------------------

rule run_fastqc:
    # INPUT: Lambda is used to get raw reads from config.yaml file. 
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    # OUTPUT: FastQC creates both reports and fastqc files. 
    output:
        "output/fastqc/{sample}_1_fastqc.html",
        "output/fastqc/{sample}_2_fastqc.html",
        "output/fastqc/{sample}_1_fastqc.zip",
        "output/fastqc/{sample}_2_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    threads: 4
    message: "Running fastQC for {wildcards.sample} on {input.r1} and {input.r2}"
    shell:
        "fastqc --threads {threads} -o output/fastqc/ {input.r1} {input.r2}"


# -----------------------------------------------------
# Rule 2: Run Assembly (SPAdes)
# -----------------------------------------------------
rule run_spades:
    # INPUT:
    # 1. Dependency on the FastQC HTML/ZIP reports to ensure QC ran first.
    # 2. Assembly reads (r1, r2) using a lambda to get the FASTQ files from the config.
    input:
        # Dependency only: Ensures FastQC is done before assembly
        qc_dep = ["output/fastqc/{sample}_1_fastqc.html", "output/fastqc/{sample}_2_fastqc.html"],
        # Actual FASTQ input for SPAdes command (from config)
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        scaffolds = "output/spades/{sample}/scaffolds.fasta"
    log:
        "logs/spades/{sample}.log"
    params:
        outdir = "output/spades/{sample}"
    threads: 4
    message: "Running SPAdes assembly for {wildcards.sample}"
    shell:
        "spades.py "
        "-k 21,33,55,77,99,127 "
        "-t {threads} "
        "--pe1-1 {input.r1} "  # Correctly uses the FASTQ file from the config lambda
        "--pe1-2 {input.r2} "  # Correctly uses the FASTQ file from the config lambda
        "--careful "
        "-o {params.outdir}"

#-----------------------------------------------------
# Rule 3: Run Assembly Assessment (QUAST)
# ----------------------------------------------------

rule run_quast:
    input:
        scaffolds = rules.run_spades.output.scaffolds
    output:
        "output/quast/{sample}/report.html"
    params:
        # Replaces os.path.join with a lambda function, as requested
        outdir = lambda wildcards: f"output/quast/{wildcards.sample}"
    threads: 4
    message: "Running QUAST for {wildcards.sample}"
    shell:
        "quast.py "
        "-o {params.outdir} "
        "-t {threads} "
        "{input.scaffolds}"
