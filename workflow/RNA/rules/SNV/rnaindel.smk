rule rnaindel:
    input:
        bam="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/rnaindel.vcf.gz"
    params:
        data_dir=config['singularity']['rnaindel'][genome_version]['data_dir']
    singularity: config['singularity']['rnaindel']['sif']
    threads:8
    shell:
        """
        rnaindel PredictIndels -i {input.bam} \
        -o {output.vcf} \
        -r {input.ref} -d {params.data_dir} -p {threads}
        """
