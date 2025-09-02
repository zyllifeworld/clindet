rule lofreq_call_up:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        Tum="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/unpaired/{sample}/lofreq.vcf"
    threads: 10
    singularity: config['singularity']['lofreq']['sif']
    shell:
        """
        lofreq call-parallel --pp-threads {threads} -f {input.ref} -o {output.vcf} {input.Tum}
        """

rule lofreq_norm_filter:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/lofreq.vcf",
    output:
        vcf="{project}/{genome_version}/results/vcf_nf/unpaired/{sample}/lofreq.vcf"
        #bcftools filter -e 'QUAL<20 | INFO/DP[0] < 10'  {input.vcf} > {output.vcf}
    shell:
        """
         cp {input.vcf}  {output.vcf}
        """
