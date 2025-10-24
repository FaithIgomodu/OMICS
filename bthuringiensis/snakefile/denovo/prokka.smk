configfile: "config/prokka.yaml"

SAMPLES = config["samples"]

rule all:
    input:
        expand(config["prokka_dir"] + "/{sample}/{sample}.gff", sample=SAMPLES)



# ------------------------------------------------------------------
# Rule 4: Run Prokka (Annotate Predicted Genes)
# ------------------------------------------------------------------

rule run_prokka:
    input:
        assembly=config["assembly_dir"] + "/scaffolds.fasta"

    output:
        gff=config["prokka_dir"] + "/{sample}/{sample}.gff"
    params:
        out_dir=config["prokka_dir"] + "/{sample}",
        prefix="{sample}"

    log:
        "logs/prokka/{sample}_prokka.log"

    threads: 4

    message: "Running Prokka for assembly: {input.assembly}"

    shell:
        """
        prokka --outdir {params.out_dir} --prefix {params.prefix} --centre X \
        --force --cpus {threads} {input.assembly} 1> {log} 2>&1
        """

