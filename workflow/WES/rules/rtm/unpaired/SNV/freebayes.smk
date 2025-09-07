rule unpaired_freebayes:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        bam="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        regions=get_sample_bed
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/freebayes.vcf"
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
        AF_THR=0.01
    threads: 10
    conda:'snake'
    shell:
        """
        freebayes -f {input.ref} -t {input.regions} --min-coverage 10 -C 3 --pooled-continuous {input.bam} > {output.vcf}
        """