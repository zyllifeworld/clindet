rule SV_unp_delly:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/unpaired/{sample}-NC.bam",
    output:
        sv="{project}/{genome_version}/results/sv/unpaired/{sample}/SV_delly_{sample}.vcf",
    params:
        ref=config['resources'][genome_version]['REFFA'],
    singularity: config['singularity']['delly']['sif']
    shell:
        """
        delly call -g {params.ref} {input.Tum} {input.NC} > {output}
        """

rule delly_unp_filter:
    input:
        vcf="{project}/{genome_version}/results/sv/unpaired/{sample}/SV_delly_{sample}.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/unpaired/{sample}/SV_delly_{sample}_filter.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
    shell:
        """
        bcftools filter 'FILTER="PASS"'  {input.vcf} > {output.vcf}
        """



