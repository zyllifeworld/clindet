
rule mutect2_call:
    input:
        Tum="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        ref=config['softwares']['rsem']['ref'][genome_version],
        pon=config['resources'][genome_version]['WES_PON'],
        germ_vcf=config['resources'][genome_version]['MUTECT2_germline_vcf'],
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/Mutect2.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
        bed=config['resources'][genome_version]['WES_BED'],
        ref='/public/ClinicalExam/lj_sih/resource/genome/human/hg38/Homo_sapiens_assembly38.fasta'
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        Mutect2 -R {params.ref} \
        --native-pair-hmm-threads {threads} \
        -I {input.Tum} \
        -O {output.vcf} \
        -pon {input.pon} \
        --germline-resource {input.germ_vcf} \
        --intervals {params.bed}
        """

rule M2_filter_unpaired:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/vcf/{sample}/Mutect2.vcf",
    output:
        vcf="{project}/{genome_version}/results/vcf/{sample}/Mutect2_filter.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
    threads: 1
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        FilterMutectCalls \
        -R {input.ref} \
        -V {input.vcf} \
        --stats {input.vcf}.stats \
        -O {output.vcf}
        """

#  export _JAVA_OPTIONS=-Djava.io.tmpdir=/public/ClinicalExam/lj_sih/projects/project_pipeline/WES/_tmp && /public/ClinicalExam/lj_sih/softwares/gatk-4.2.6.1/gatk          Mutect2 -R /public/ClinicalExam/lj_sih/resource/genome/human/hg38/Homo_sapiens_assembly38.fasta         --native-pair-hmm-threads 1         -I MMRNA/hg38/results/dedup/MM-066.split.bam         -O MMRNA/hg38/results/vcf/MM-066/Mutect2.vcf         -pon /public/ClinicalExam/lj_sih/projects/project_clindet/reference/hg38/1000g_pon.hg38.vcf.gz         --germline-resource /public/ClinicalExam/lj_sih/projects/project_clindet/reference/hg38/af-only-gnomad.hg38.vcf.gz         --intervals /public/ClinicalExam/lj_sih/resource/genome/human/hg38/hg38.bed