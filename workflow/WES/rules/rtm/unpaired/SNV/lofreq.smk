rule lofreq_somatic_unpaired:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions=get_sample_bed,
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam"
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/Lofreq/{sample}.lofreq.vcf",
        # indel="{project}/{genome_version}/results/vcf/paired/{sample}/lofreq/out_somatic_final.indels.vcf.gz"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/unpaired/{sample}/Lofreq/out_",
        dbsnp="",
    threads: 10
    singularity: config['singularity']['lofreq']['sif']
    shell:
        """
        lofreq call-parallel --pp-threads {threads} -f {params.ref} -l {input.regions} -o {output.vcf} {input.Tum} -s
        """

rule unpair_lofreq_filter:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/Lofreq/{sample}.lofreq.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/lofreq.vcf"
        #bcftools filter -e 'QUAL<20 | INFO/DP[0] < 10'  {input.vcf} > {output.vcf}
    shell:
        """
         cp {input.vcf}  {output.vcf}
        """