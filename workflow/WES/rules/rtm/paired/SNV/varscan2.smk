## varscan2 call #####
### based on https://pubmed.ncbi.nlm.nih.gov/25553206/
rule varscan2_mpileup:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        normal="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        tumor="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        regions=get_sample_bed,
    output:
        normal=temp("{project}/{genome_version}/results/recal/paired/{sample}-NC.mp"),
        tumor=temp("{project}/{genome_version}/results/recal/paired/{sample}-T.mp")
    threads: 2
    conda: config['conda']['clindet_main']
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.varscan_mpileup.benchmark.txt"
    shell:
        """
        samtools mpileup -q 1 -Q 1 -f {input.ref} -l {input.regions} {input.normal} > {output.normal}
        samtools mpileup -q 1 -Q 1 -f {input.ref} -l {input.regions} {input.tumor} > {output.tumor}
        """

rule varscan2_call:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        normal="{project}/{genome_version}/results/recal/paired/{sample}-NC.mp",
        tumor="{project}/{genome_version}/results/recal/paired/{sample}-T.mp",
        regions=get_sample_bed,
    output:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.vcf",
        indel="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.vcf"
    log:
        "{project}/{genome_version}/logs/varscan2/paired/{sample}.log"
    conda: config['conda']['clindet_main']
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.varscan_somatic.benchmark.txt"
    shell:
        """
        varscan somatic {input.normal} {input.tumor} --output-snp {output.snp} --output-indel {output.indel} --output-vcf 1 --strand-filter 1
        """

rule varscan2_processSomatic:
    input:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.vcf",
        indel="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.vcf"
    output:
        snp_som="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.Somatic.vcf",
        snp_som_hc="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.Somatic.hc.vcf",
        snp_loh="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.LOH.vcf",
        snp_loh_hc="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.LOH.hc.vcf",
        snp_germ="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.Germline.vcf",
        snp_germ_hc="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.Germline.hc.vcf",
        indel_som="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.Somatic.vcf",
        indel_som_hc="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.Somatic.hc.vcf",
        indel_loh="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.LOH.vcf",
        indel_loh_hc="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.LOH.hc.vcf",
        indel_germ="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.Germline.vcf",
        indel_germ_hc="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.Germline.hc.vcf",
    threads: 1
    conda: config['conda']['clindet_main']
    shell:
        """
        varscan processSomatic {input.snp}
        varscan processSomatic {input.indel}
        """

rule varscan2_som_filter:
    input:
        indel="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.vcf",
        snp_som_hc="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.Somatic.hc.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.Somatic.hc.filter.vcf"
    conda: config['conda']['clindet_main']
    shell:
        """
        varscan somaticFilter {input.som_hc} --indel-file {input.indel} --output-file {output.vcf}
        """

rule varscan2_merge_somatic:
    input:
        snp_som_hc_filter="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.Somatic.hc.filter.vcf",
        indel_som_hc="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.Somatic.hc.vcf",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.vcf"
    threads: 1
    params:
        caller='varscan2'
    conda: config['conda']['clindet_main']
    shell:
        """
        bgzip {input.snp_som_hc_filter} -k -o {input.snp_som_hc_filter}.gz
        tabix {input.snp_som_hc_filter}.gz
        bgzip {input.indel_som_hc} -k -o {input.indel_som_hc}.gz
        tabix {input.indel_som_hc}.gz
        bcftools concat -a  {input.snp_som_hc_filter}.gz  {input.indel_som_hc}.gz | bcftools view -o {output.vcf}
        """
