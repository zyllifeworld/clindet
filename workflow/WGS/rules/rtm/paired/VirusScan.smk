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