configfile: "config/annotate.yaml"

SAMPLES = config["samples"]

rule all:
    input:
        expand(config["barrnap_dir"] + "/rrna.gff", sample=SAMPLES),
        expand(config["prodigal_dir"] + "/prodigal.gff", sample=SAMPLES),

# ------------------------------------------------------------------
# Rule 1: Run Prodigal (Predict Genes)
# ------------------------------------------------------------------

rule run_prodigal:
    input:
        scaffolds="output/spades/{sample}/scaffolds.fasta"

    output:
        gff="output/prodigal/{sample}/prodigal.gff",
        translated_proteins="output/prodigal/{sample}/prodigal_translated.out",
        genes="output/prodigal/{sample}/prodigal_geneseqs.out",
        stats="output/prodigal/{sample}/prodigal_table.out",

    params:
        out_dir="output/prodigal/{sample}",
        base_name="prodigal",

    log:
        "logs/prodigal/prodigal_{sample}.log"

    threads: 4

    message: "Running Prodigal for {input}"

    shell:
        """
        mkdir -p {params.out_dir}
        prodigal -i {input} \
        -o {params.out_dir}/{params.base_name}.out -s {output.stats} -f gff -d {output.genes} \
        -a {output.translated_proteins} > {output.gff} 2>> {log}
        """




# ------------------------------------------------------------------
# Rule 3: Run Barrnap (Predict rRNA's)
# ------------------------------------------------------------------

rule run_barrnap:
    input:
        scaffolds="output/spades/{sample}/scaffolds.fasta"

    output:
        gff="output/barrnap/{sample}/rrna.gff"

    params:
        out_dir="output/barrnap/{sample}"

    log:
        "logs/barrnap/barrnap_{sample}.log"

    threads: 4

    message: "Running Barrnap for {input}"

    shell:
        """
        mkdir -p {params.out_dir}
        barrnap --threads {threads} --kingdom bac {input} 1> {output.gff} 2> {log}
