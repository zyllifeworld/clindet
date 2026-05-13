# lancet2 only support bed file with three cols, so i will generate this file first
rule lancet2_somatic_call:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions=get_sample_bed,
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        vcf_gz="{project}/{genome_version}/results/vcf/paired/{sample}/lancet2/lancet2.vcf.gz",
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/lancet2/lancet2.vcf"
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','lancet2','sif'],image_url = config['singularity']['lancet2']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.lancet2.benchmark.txt"
    shell:
        """
        [ ! -f {input.regions}.lancet2 ] && cut -f 1-3 {input.regions} > {input.regions}.lancet2
        Lancet2 pipeline --reference {input.reference} -T {threads} \
        --tumor {input.Tum} --normal {input.NC} \
        --bed-file {input.regions}.lancet2 --out-vcfgz {output.vcf_gz}
        zcat {output.vcf_gz} > {output.vcf}
        """