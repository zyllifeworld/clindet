rule map_reads:
    input:
        R1="{project}/{genome_version}/results/trimmed/{sample}-{group}_R1.fastq.gz",
        R2="{project}/{genome_version}/results/trimmed/{sample}-{group}_R2.fastq.gz"
    output:
        "{project}/{genome_version}/results/mapped/{sample_type}/{sample}-{group}.sorted.bam",
    params:
        ref=config['resources'][genome_version]['REFFA'],
        rg=r"@RG\tID:{sample}_{group}\tPL:Illumia\tSM:{sample}_{group}\tLB:DNA-Seq"
    threads: 30 
    conda:
        config['softwares']['samtools']['conda']
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mapping/{sample}.mapping.benchmark.txt"
    shell:
        """ {config[softwares][bwa][mem][call]} -t 30 -MR '{params.rg}' \
        {params.ref} \
        {input.R1} {input.R2} | samtools fixmate -O bam - - | \
        samtools sort -@ 30 -O bam -o {output}
        """

rule mark_duplicates:
    input:
        "{project}/{genome_version}/results/mapped/{sample_type}/{sample}-{group}.sorted.bam",
    output:
        bam=temp("{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam"),
        metrics="{project}/{genome_version}/results/qc/dedup/{sample_type}/{sample}-{group}.metrics.txt"
    params:
        temp_directory=config['params']['java']['temp_directory']
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mapping/{sample}.markdup.benchmark.txt"
    shell:
        """
        {config[softwares][gatk4][MarkDuplicates][call]} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
        -I {input} \
        -O {output.bam} \
        -M {output.metrics}
        """

recal = True 
### for faster run, may consider not run applyBQSR, but i will keep this step in WES, you can customize as you own.
## For Noveseq data, don't do this step, meaningless!
## Aslo this step will not significantly import downstream analysis see: https://www.biostars.org/p/9605712/ and anywhere else.

if recal:
    ## if reacal, let dedup bam as temp file to save space
    rule recalibrate_base_qualities:
        input:
            bam="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam",
            #bai="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam.bai",
            ref=config['resources'][genome_version]['REFFA'],
            known_1=config['resources']['varanno'][genome_version]['KNOWN_SITES1'],
            known_2=config['resources']['varanno'][genome_version]['KNOWN_SITES2']
        output:
            recal_table="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.grp",
        params:
            temp_directory=config['params']['java']['temp_directory']
        benchmark:
            "{project}/{genome_version}/results/benchmarks/mapping/{sample}.recal.benchmark.txt"
        shell:
            """ 
                {config[softwares][gatk4][BaseRecalibrator][call]} --use-original-qualities -R {input.ref} \
                -I {input.bam} \
                -O {output.recal_table} \
                --known-sites {input.known_1} \
                --known-sites {input.known_2}
            """

    rule apply_base_quality_recalibration:
        input:
            bam="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam",
            #bai="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam.bai",
            ref=config['resources'][genome_version]['REFFA'],
            recal_table="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.grp",
        output:
            bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
            bai="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam.bai",
        params:
            temp_directory=config['params']['java']['temp_directory']
        benchmark:
            "{project}/{genome_version}/results/benchmarks/mapping/{sample}.applybqsr.benchmark.txt"
        shell:
            """{config[softwares][gatk4][ApplyBQSR][call]} --add-output-sam-program-record -use-original-qualities \
                --bqsr-recal-file {input.recal_table} \
                -R {input.ref} \
                -I {input.bam} \
                -O {output.bam}  
                samtools index {output.bam}
            """
    else:
        rule recal_link:
            input:
                bam="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam",
                bai="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bai"
            output:
                bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
                bai="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam.bai"
            params:
                temp_directory=config['params']['java']['temp_directory']
            shell:
                """
                ln -s $(realpath {input.bam}) {output.bam}
                ln -s $(realpath {input.bai}) {output.bai}
                """


rule bed_to_interval_list:
    input:
        bed=get_sample_bed,
        dict=config['resources'][genome_version]['REFFA_DICT'],
    output:
        it="{project}/{genome_version}/picard/{sample}.interval_list"
    params:
        extra="--SORT true",  # sort output interval list before writing
        picard=config['softwares']['picard']['call'],
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """
            export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && \
            java -jar {params.picard} BedToIntervalList --INPUT {input.bed} \
            --SEQUENCE_DICTIONARY {input.dict} {params.extra} \
            --OUTPUT {output.it}
        """


rule picard_collect_wes:
    input:
        bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
        ref=config['resources'][genome_version]['REFFA'],
        reference=config['resources'][genome_version]['REFFA'],
        bait_intervals="{project}/{genome_version}/picard/{sample}.interval_list",
        target_intervals="{project}/{genome_version}/picard/{sample}.interval_list",
    output:
        hs="{project}/{genome_version}/results/stats/{sample_type}/wes_metrics/{sample}-{group}.txt"
    resources:
        mem_mb=1024
    params:
        gatk4=config['softwares']['gatk4']['call'],
        picard=config['softwares']['picard']['call'],
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && \
        java -jar {params.picard} CollectHsMetrics \
        --INPUT {input.bam} \
        --OUTPUT {output.hs} \
        --REFERENCE_SEQUENCE {input.reference} \
        --BAIT_INTERVALS {input.bait_intervals} \
        --TARGET_INTERVALS {input.target_intervals}
        """
