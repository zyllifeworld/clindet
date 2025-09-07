## varscan2 call #####
varscan2_conda = 'snake'
rule varscan2_mpileup_unpaired:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        tumor="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        regions=get_sample_bed,
    output:
        tumor="{project}/{genome_version}/results/recal/unpaired/{sample}-T.mp"
    threads: 2
    conda:varscan2_conda
    shell:
        """
        samtools mpileup -q 1 -Q 1 -f {input.ref} -l {input.regions} {input.tumor} > {output.tumor}
        """

rule varscan2_call_unpaired_snp:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        tumor="{project}/{genome_version}/results/recal/unpaired/{sample}-T.mp",
        regions=get_sample_bed,
    output:
        snp="{project}/{genome_version}/results/vcf/unpaired/{sample}/varscan/varscan2.snp.vcf"
    conda:varscan2_conda
    shell:
        """
        varscan mpileup2snp   {input.tumor} --output-vcf 1 > {output.snp}
        """


rule varscan2_call_unpaired_indel:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        tumor="{project}/{genome_version}/results/recal/unpaired/{sample}-T.mp",
        regions=get_sample_bed,
    output:
        indel="{project}/{genome_version}/results/vcf/unpaired/{sample}/varscan/varscan2.indel.vcf"
    conda:varscan2_conda
    shell:
        """
        varscan mpileup2indel {input.tumor} --output-vcf 1 > {output.indel}
        """

rule varscan2_merge_unpaired:
    input:
        snp="{project}/{genome_version}/results/vcf/unpaired/{sample}/varscan/varscan2.snp.vcf",
        indel="{project}/{genome_version}/results/vcf/unpaired/{sample}/varscan/varscan2.indel.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/varscan2.vcf",
        name_change="{project}/{genome_version}/results/vcf/unpaired/{sample}/varscan/sample_name.txt"
    threads: 1
    params:
        name_change="{project}/{genome_version}/results/vcf/unpaired/{sample}/varscan/sample_name.txt"
    shell:
        """
        echo -e "Sample1\t{wildcards.sample}_T" > {output.name_change}
        bgzip {input.snp} -k -o {input.snp}.gz
        tabix {input.snp}.gz
        bgzip {input.indel} -k -o {input.indel}.gz
        tabix {input.indel}.gz
        bcftools concat -a  {input.snp}.gz  {input.indel}.gz | bcftools view | bcftools reheader -s {output.name_change}  -o {output.vcf}
        """

