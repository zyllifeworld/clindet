rule M2_ST:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        vcf=config['resources'][genome_version]['MUTECT2_germline_vcf'],
        ref=config['resources'][genome_version]['REFFA']
    output:
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table"
    params:
        bed=get_sample_bed,
        temp_directory=config['params']['java']['temp_directory']
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk \
        GetPileupSummaries -R {input.ref}  \
        -I {input.Tum} \
        -V  {input.vcf} \
        -L {params.bed} \
        -O {output.table}
        """

rule M2_SNC:
    input:
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        vcf=config['resources'][genome_version]['MUTECT2_germline_vcf']
    output:
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table"
    params:
        temp_directory=config['params']['java']['temp_directory'],
        bed=get_sample_bed
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk \
        GetPileupSummaries -R {input.ref}  \
        -I {input.NC} \
        -V {input.vcf} \
        -L {params.bed} \
        -O {output.table}
        """


rule M2_contam:
    input:
        T="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table",
        NC="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        seg="{project}/{genome_version}/results/recal/contam/{sample}/{sample}_segments.table",
        ctam="{project}/{genome_version}/results/recal/contam/{sample}/{sample}_calculatecontamination.table"
    params:
        temp_directory=config['params']['java']['temp_directory'],
        bed=get_sample_bed
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk \
        CalculateContamination \
        -I {input.T} \
        -matched {input.NC} \
        -tumor-segmentation {output.seg} \
        -O {output.ctam} 
        """


rule mutect2:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed,
        germ_vcf=config['resources'][genome_version]['MUTECT2_germline_vcf'],
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2_raw.vcf"
    params:
        temp_directory=config['params']['java']['temp_directory'],
        # pon=config['resources'][genome_version]['WES_PON'],
        pon=lambda wildcards: get_config_value(config,['resources', wildcards.genome_version, 'WES_PON'], params='-pon', default=""),
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk \
        Mutect2 -R {input.ref} \
        --native-pair-hmm-threads {threads} \
        -I {input.Tum} \
        -I {input.NC} \
        -O {output.vcf} \
        -normal {wildcards.sample}_NC \
        {params.pon} \
        --germline-resource {input.germ_vcf} \
        --intervals {input.bed}
        """

rule M2_filter:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed,
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2_raw.vcf",
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table",
        seg="{project}/{genome_version}/results/recal/contam/{sample}/{sample}_segments.table",
        ctam="{project}/{genome_version}/results/recal/contam/{sample}/{sample}_calculatecontamination.table"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2.vcf"
    params:
        temp_directory=config['params']['java']['temp_directory']
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk \
        FilterMutectCalls \
        -R {input.ref} \
        -V {input.vcf} \
        --contamination-table {input.ctam} \
        --stats {input.vcf}.stats \
        --tumor-segmentation {input.seg} \
        -O {output.vcf}
        """
