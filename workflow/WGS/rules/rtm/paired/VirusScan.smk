# virusbreakend \
#   -r /public/home/lijf/env/genome/broad/hg19/ucsc.hg19.fasta \
#   -j /usr/local/share/gridss-2.13.2-2/gridss.jar \
#   -o /public/ClinicalExam/lj_sih/projects/project_pipeline/WGS/MMWGSPE300_RJ/hg19/results/virus/virusbreaked/MM-215/MM-215.vcf \
#   --db /public/ClinicalExam/lj_sih/resource/genome/human/virusbreakenddb_20210401 \
#   /public/ClinicalExam/lj_sih/projects/project_pipeline/WGS/MMWGSPE300_RJ/hg19/results/recal/paired/MM-215-T.bam

rule virusbreakend:
    input:
      ref=config['resources'][genome_version]['REFFA'],
      bam="{project}/{genome_version}/results/recal/paired/{sample}-T.bam"
    output:
      vcf="{project}/{genome_version}/results/virus/virusbreaked/{sample}/{sample}.vcf"
    params:
      db="/public/ClinicalExam/lj_sih/resource/genome/human/virusbreakenddb_20210401",
      wd="{project}/{genome_version}/results/virus/virusbreaked/{sample}"
    threads: 10
    singularity: config['singularity']['gridss']['sif']
    shell:
      """
        virusbreakend \
        -r {input.ref} -t {threads} --minreads 30 \
        -j /usr/local/share/gridss-2.13.2-2/gridss.jar \
        -o {output.vcf}  -w {params.wd}\
        --db {params.db} --minviralcoverage 0 \
        {input.bam}
      """