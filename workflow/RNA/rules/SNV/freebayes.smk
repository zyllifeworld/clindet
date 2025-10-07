rule unpaired_freebayes:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        bam="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        bam_bai="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam.bai",
        regions=config['resources'][genome_version]['WES_BED'],
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/freebayes.vcf",
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 10
    shell:
        """
        freebayes -f {input.ref} -t {input.regions} --min-coverage 10 -C 3 --pooled-continuous {input.bam} > {output.vcf}
        """


rule norm_filter_freebayes:
    input:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/freebayes.vcf",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/vcf_nf/unpaired/{sample}/freebayes.vcf",
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
    shell:
        """
        bcftools norm -f {input.ref} -m +both -Ov {input.vcf} | bcftools filter -e 'QUAL<20 | FORMAT/DP[0] < 10' > {output.vcf}
        """