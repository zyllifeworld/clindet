rule cnv_facets:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        png="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}.cnv.png"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}",
        common_vcf=config['resources'][genome_version]['common_vcf'],
    threads: 8
    conda:config['softwares']['facets']['conda']
    shell:
        """
        cnv_facets.R -t {input.Tum} -n {input.NC} -vcf {params.common_vcf} -o {params.wd}
        """
# cnv_facets.R -t MMWES/hg19/results/recal/paired/MM-099-T.bam -n MMWES/hg19/results/recal/paired/MM-099-NC.bam -vcf /public/ClinicalExam/lj_sih/resource/genome/human/hg19/00-common_all_hg19.vcf.gz -o MMWES/hg19/results/cnv/paired/facets/MM-099/MM-099