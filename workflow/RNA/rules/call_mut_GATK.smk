rule mark_duplicates_score:
    input:
       bam="{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam"
    #    bai="{project}/{genome_version}/results/dedup/{sample}.addRG.bam.bai"
    output:
        score="{project}/{genome_version}/results/mut/dedup/{sample}.score.gz",
    params:
        temp_directory=config['params']['java']['temp_directory'],
        metrics="{project}/{genome_version}/results/qc/dedup/{sample}.metrics.txt"
    threads:10
    shell:
        """
        {config[softwares][gatk4][MarkDuplicates][call]} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
        -I {input} \
        -O {output.bam} \
        -M {params.metrics}
        """

rule mark_duplicates:
    input:
        bam="{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam",
    output:
        bam=temp("{project}/{genome_version}/results/mut/dedup/{sample}.rmdep.bam"),
    params:
        temp_directory=config['params']['java']['temp_directory'],
    threads:10
    shell:
        """
        {config[softwares][gatk4][MarkDuplicates][call]} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
        -I {input} \
        -O {output.bam} \
        -M {params.metrics}
        """

rule SplitNCigarReads:
    input:
        bam="{project}/{genome_version}/results/mut/dedup/{sample}.rmdep.bam",
        # bai="{project}/{genome_version}/results/dedup/{sample}.rmdep.bam.bai",
        #bai="{project}/{genome_version}/{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam.bai",
        rsem_ref=config['resources'][genome_version]['REFFA']
    output:
        bam="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam"
    params:
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """ 
            export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {config[softwares][gatk4][call]} SplitNCigarReads \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.bam}
        """

## will not recalibrate base qualities
# rule recalibrate_base_qualities:
#     input:
#         bam="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
#         #bai="{project}/{genome_version}/{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam.bai",
#         ref=config['resources'][genome_version]['REFFA'],
#         known_1=config['resources']['varanno'][genome_version]['KNOWN_SITES1'],
#         known_2=config['resources']['varanno'][genome_version]['KNOWN_SITES2']
#     output:
#         recal_table="{project}/{genome_version}/results/mut/recal/{sample}.grp",
#     params:
#         temp_directory=config['params']['java']['temp_directory']
#     shell:
#         """ 
#             {config[softwares][gatk4][BaseRecalibrator][call]} --use-original-qualities -R {input.ref} \
#             -I {input.bam} \
#             -O {output.recal_table} \
#             --known-sites {input.known_1} \
#             --known-sites {input.known_2}
#         """


# rule apply_base_quality_recalibration:
#     input:
#         bam="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
#         #bai="{project}/{genome_version}/{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam.bai",
#         ref=config['resources'][genome_version]['REFFA'],
#         recal_table="{project}/{genome_version}/results/mut/recal/{sample}.grp",
#     output:
#         bam="{project}/{genome_version}/results/recal/{sample}.bam",
#         bam_bai="{project}/{genome_version}/results/recal/{sample}.bam.bai"
#     conda:
#         config['softwares']['samtools']['conda']
#     params:
#         temp_directory=config['params']['java']['temp_directory']
#     shell:
#         """{config[softwares][gatk4][ApplyBQSR][call]} --add-output-sam-program-record -use-original-qualities \
#             --bqsr-recal-file {input.recal_table} \
#             -R {input.ref} \
#             -I {input.bam} \
#             -O {output.bam}  
#             samtools index {output.bam}
#         """
