if redup:
    rule mark_duplicates:
        input:
            bam="{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam"
        output:
            bam=temp("{project}/{genome_version}/results/mut/dedup/{sample}.rmdep.bam"),
        params:
            temp_directory=config['params']['java']['temp_directory'],
        threads:10
        singularity:
            flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
        shell:
            """
            export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && \
            gatk MarkDuplicates --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
            -I {input} \
            -O {output.bam} \
            -M {params.metrics}
            """
else:
    rule link_bam:
        input:
            bam="{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam",
            bai="{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam.bai"
        output:
            bam="{project}/{genome_version}/results/mut/dedup/{sample}.rmdep.bam",
            bai="{project}/{genome_version}/results/mut/dedup/{sample}.rmdep.bai",
        threads:1
        shell:
            """
            ln -s $(realpath {input.bam}) {output.bam}
            ln -s $(realpath {input.bai}) {output.bai}
            """

rule SplitNCigarReads:
    input:
        bam="{project}/{genome_version}/results/mut/dedup/{sample}.rmdep.bam",
        ref=config['resources'][genome_version]['REFFA']
    output:
        bam="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam"
    params:
        temp_directory=config['params']['java']['temp_directory']
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    shell:
        """ 
            export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk SplitNCigarReads \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.bam}
        """