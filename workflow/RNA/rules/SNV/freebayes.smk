rule unpaired_freebayes:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        samples="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        indexes="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam.bai",
    output:
        "{project}/{genome_version}/results/mut/vcf/{sample}/freebayes.vcf",
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 10
    wrapper:
        "v1.7.0/bio/freebayes"


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