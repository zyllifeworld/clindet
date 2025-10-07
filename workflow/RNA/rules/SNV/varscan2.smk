## varscan2 call #####
varscan2_conda = 'snake'
rule varscan2_mpileup_unpaired:
    input:
        regions=config['resources'][genome_version]['WES_BED'],
        tumor="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        tumor="{project}/{genome_version}/results/mut/dedup/{sample}-T.mp"
    threads: 2
    conda:varscan2_conda
    shell:
        """
        samtools mpileup -q 1 -Q 1 -f {input.ref} -l {input.regions} {input.tumor} > {output.tumor}
        """

rule varscan2_call_unpaired_snp:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        tumor="{project}/{genome_version}/results/mut/dedup/{sample}-T.mp",
        regions=config['resources'][genome_version]['WES_BED'],
    output:
        snp="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/varscan2.snp.vcf"
    conda:varscan2_conda
    shell:
        """
        varscan mpileup2snp  {input.tumor} --min-coverage 10 --min-var-freq 0.20 --p-value 0.05 --output-vcf 1 > {output.snp}
        """


rule varscan2_call_unpaired_indel:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        tumor="{project}/{genome_version}/results/mut/dedup/{sample}-T.mp",
        regions=config['resources'][genome_version]['WES_BED'],
    output:
        indel="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/varscan2.indel.vcf"
    conda:varscan2_conda
    shell:
        """
        varscan mpileup2indel {input.tumor} --min-coverage 10 --min-var-freq 0.10 --p-value 0.10 --output-vcf 1 > {output.indel}
        """

rule varscan2_filter_snp:
    input:
        snp="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/varscan2.snp.vcf",
        indel="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/varscan2.indel.vcf"
    output:
        snp="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/varscan2.snp.filter.vcf"
    conda:varscan2_conda
    shell:
       """
        varscan filter {input.snp} --indel-file {input.indel} --output-file {output.snp}
       """

rule varscan2_filter_indel:
    input:
        indel="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/varscan2.indel.vcf"
    output:
        indel="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/varscan2.indel.filter.vcf"
    conda:varscan2_conda
    shell:
       """
        varscan filter {input.indel} --min-reads2 4 --min-var-freq 0.15 --p-value 0.05 --output-file {output.indel}
       """

rule varscan2_merge_unpaired:
    input:
        snp="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/varscan2.snp.filter.vcf",
        indel="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/varscan2.indel.filter.vcf"
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/varscan2.vcf",
        name_change="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/sample_name.txt"
    threads: 1
    params:
        name_change="{project}/{genome_version}/results/mut/vcf/{sample}/varscan/sample_name.txt"
    shell:
        """
        echo -e "Sample1\t{wildcards.sample}" > {output.name_change}
        bgzip {input.snp} -k -o {input.snp}.gz
        tabix {input.snp}.gz
        bgzip {input.indel} -k -o {input.indel}.gz
        tabix {input.indel}.gz
        bcftools concat -a  {input.snp}.gz  {input.indel}.gz | bcftools view | bcftools reheader -s {output.name_change}  -o {output.vcf}
        """
