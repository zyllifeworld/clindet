rule sequenza_bam2seqz:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        gc=lambda wildcards: get_config_value(config,['singularity', 'sequenza', wildcards.genome_version, 'gc'],default=""),
    output:
        seqz="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}.seqz.gz",
    threads: 8
    singularity:
        flexible_container_img(config,['singularity','sequenza','sif'],image_url = config['singularity']['sequenza']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/cnv/{sample}.sequenza_bam2seqz.benchmark.txt"
    shell:
        """
        sequenza-utils bam2seqz \
        --normal {input.NC} \
        --tumor {input.Tum} \
        --fasta {input.ref} -gc {input.gc} \
        --output {output.seqz}
        """

rule sequenza_seqz_binning:
    input:
        seqz="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}.seqz.gz"
    output:
        bin_seqz="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}_bin50.seqz.gz",
    singularity:
        flexible_container_img(config,['singularity','sequenza','sif'],image_url = config['singularity']['sequenza']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/cnv/{sample}.sequenza_seqz_binning.benchmark.txt"
    shell:
        """
        sequenza-utils seqz_binning -w 50 --seqz {input.seqz} -o {output.bin_seqz}
        """


rule sequenza_call:
    input:
        bin_seqz="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}_bin50.seqz.gz",
    output:
        segment="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}_segments.txt",
    params:
        wd="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}",
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    benchmark:
        "{project}/{genome_version}/results/benchmarks/cnv/{sample}.sequenza_call.benchmark.txt"
    script:
        "../../../../scripts/sequenza.R"