rule lofreq_somatic:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions=get_sample_bed,
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/lofreq/out_tumor_stringent.snvs.vcf.gz",
        # indel="{project}/{genome_version}/results/vcf/paired/{sample}/lofreq/out_somatic_final.indels.vcf.gz"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/paired/{sample}/lofreq/out_",
        dbsnp="",
    threads: 10
    singularity: config['singularity']['lofreq']['sif']
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.lofreq.benchmark.txt"
    shell:
        """
        lofreq somatic -n {input.NC} -t {input.Tum} -f {params.ref} -l {input.regions} --threads {threads} -o {params.out_prefix}
        """

rule lofreq_norm_filter:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/lofreq/out_somatic_final.snvs.vcf.gz",
        indel="{project}/{genome_version}/results/vcf/paired/{sample}/lofreq/out_somatic_final.indels.vcf.gz"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/lofreq.vcf"
    conda: config['conda']['clindet_main']
    shell:
        """
        bcftools concat -a {input.snp} {input.indel} | bcftools filter -e 'QUAL<20 | INFO/DP[0] < 10' -o {output.vcf}
        """