rule STAR_1_pass:
    input:
        # R1="{project}/{genome_version}/results/trimmed/{sample}_R1.fastq.gz",
        # R2="{project}/{genome_version}/results/trimmed/{sample}_R2.fastq.gz",
        unpack(get_rna_fastq)
    output:
        pass1="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_pass1.log",
        sj_filt='{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/SJ.filt'
    params:
        out_dir="{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/", # "/"" mustq in the config string
        ref=config['resources'][genome_version]['REFFA'],
        star_index=config['softwares']['star']['index'][genome_version],
        rg=r"ID:{sample} PL:ILLUMINA.NovaSeq LB:RNA-Seq SM:{sample}"
    threads: 10
    conda:
        config['softwares']['star']['conda']
    shell:
        """ 
        STAR --genomeDir {params.star_index} --runThreadN={threads} \
            --outSAMtype None --outFileNamePrefix {params.out_dir} \
            --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat
        awk '$1~/chr[1-2XY]/ && $6==0 && $5>0 && $7>0' {params.out_dir}SJ.out.tab > {output.sj_filt}
        touch {output.pass1}
        """



rule STAR_arriba_map:
    input:
        unpack(get_rna_fastq),
        sj_filt='{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/SJ.filt'
    output:
        bam=temp("{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam"),
        stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log",
        um_fq1="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_unmapped_R1.fq",
        um_fq2="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_unmapped_R2.fq"
    params:
        out_dir="{project}/{genome_version}/results/mapped/STAR/{sample}/", # "/"" must in the config string
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
        star_index=config['softwares']['star']['index'][genome_version],
        rg=r"ID:{sample} PL:ILLUMINA.NovaSeq LB:RNA-Seq SM:{sample}"
    threads: 10
    conda:
        config['softwares']['star']['conda']
    shell:
        """ 
        STAR --genomeDir {params.star_index} --runThreadN={threads} \
            --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.out_dir} \
            --outReadsUnmapped Fastx \
            --outSAMattrRGline '{params.rg}' --sjdbFileChrStartEnd {input.sj_filt}\
            --sjdbGTFfile {params.gtf} \
            --outFilterMultimapNmax 50 \
            --peOverlapNbasesMin 10 \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimSegmentMin 10 \
            --chimOutType WithinBAM HardClip \
            --chimJunctionOverhangMin 10 \
            --chimScoreDropMax 30 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50 \
            --readFilesIn {input.R1} {input.R2} --readFilesCommand gunzip -c

        mv {params.out_dir}/Aligned.sortedByCoord.out.bam {output.bam}
        mv {params.out_dir}/Unmapped.out.mate1 {output.um_fq1}
        mv {params.out_dir}/Unmapped.out.mate2 {output.um_fq2}

        touch {output.stamp}
        """

# maybe sort by samtools,STAR always face RAM error
rule STAR_mut_map:
    input:
        unpack(get_rna_fastq),
    output:
        bam="{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam",
        unsort_bam="{project}/{genome_version}/results/mut/STAR/{sample}/Aligned.out.bam",
    params:
        out_dir="{project}/{genome_version}/results/mut/STAR/{sample}/", # "/"" must in the config string
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
        star_index=config['softwares']['star']['index'][genome_version],
        rg="ID:{sample} -r PL:ILLUMINA.NovaSeq -r LB:RNA-Seq -r SM:{sample}",
        sort_mem_per_thread='1G'
    threads: 10
    conda:
        config['softwares']['star']['conda']
    shell:
        """ 
        STAR --genomeDir {params.star_index} --runThreadN={threads} \
            --outSAMtype BAM Unsorted --outFileNamePrefix {params.out_dir} \
            --outReadsUnmapped Fastx \
            --sjdbGTFfile {params.gtf} \
            --alignIntronMax 1000000 \
            --alignIntronMin 20 \
            --alignMatesGapMax 1000000 \
            --alignSJDBoverhangMin 1 \
            --alignSJoverhangMin 8 \
            --alignSoftClipAtReferenceEnds Yes \
            --chimJunctionOverhangMin 15 \
            --chimMainSegmentMultNmax 1 \
            --chimOutJunctionFormat 1 \
            --chimSegmentMin 15 \
            --limitSjdbInsertNsj 1200000 \
            --outFilterIntronMotifs None \
            --outFilterMatchNminOverLread 0.33 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.1 \
            --outFilterMultimapNmax 20 \
            --outFilterScoreMinOverLread 0.33 \
            --twopassMode Basic \
            --outSAMmapqUnique 60 \
            --readFilesIn {input.R1} {input.R2} --readFilesCommand gunzip -c
        samtools addreplacerg -@ -{threads} -o {params.out_dir}/addRG.bam {params.out_dir}/Aligned.out.bam
        samtools sort -n -@ {threads} -m 1G -o {output.bam} {params.out_dir}/addRG.bam
        samtools index {output.bam}

        rm {params.out_dir}/addRG.bam
        """


rule cal_exp_RSEM:
    input:
        unpack(get_rna_fastq),
    params:
        rsem_ref=config['softwares']['rsem']['index'][genome_version],
        result_prefix="{project}/{genome_version}/results/summary/RSEM"
    threads: 10
    conda:
        config['softwares']['rsem']['conda']
    output:
        genes="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.genes.results",
        isoforms="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.isoforms.results",
        bam=temp("{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.bam"),
        tx_bam=temp("{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.transcript.bam"),
        tx_sort_bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.transcript.sorted.bam"
        # bam_sort=temp("{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.sorted.bam"),
    shell:
        """
            rsem-calculate-expression  --paired-end \
            --star -p {threads} --star-output-genome-bam \
            --star-gzipped-read-file  --sort-bam-by-coordinate \
            {input.R1} {input.R2} \
            {params.rsem_ref} {params.result_prefix}/{wildcards.sample}/{wildcards.sample}
        """

rule RSEM_sort_genome:
    input:
        bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.bam"
    threads: 10
    conda:
        config['softwares']['rsem']['conda']
    output:
        bam_sort=temp("{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.sorted.bam"),
    shell:
        """
        samtools sort -@ {threads} -m 1G -o {output.bam} {input.bam}
        samtools index {output.bam}
        """

# rule RSEM_bam2bigwig:
#     input:
#         bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.sort.bam"
#     output:
#         bw="{project}/{genome_version}/results/RSEM/bigwig/{sample}/{sample}.bw"
#     conda:'deeptools'
#     threads:10
#     shell:
#         """
#         bamCoverage -b {input.bam} -o {output.bw} -p {threads} --normalizeUsing RPKM
#         """


###### kallisto quanto
rule kallisto:
    input:
        unpack(get_rna_fastq),
    params:
        index=config['softwares']['kallisto']['index'][genome_version],
        result_prefix="{project}/{genome_version}/results/summary/kallisto/{sample}"
    threads: 10
    conda:
        config['softwares']['kallisto']['conda']
    output:
        tsv="{project}/{genome_version}/results/summary/kallisto/{sample}/abundance.tsv",
    shell:
        """
        kallisto quan -i {params.index} -o {params.result_prefix} {input.R1} {input.R2} -t {threads}
        """
### samlom 
rule salmon:
    input:
        unpack(get_rna_fastq),
    params:
        index=config['softwares']['salmon']['index'][genome_version],
        result_prefix="{project}/{genome_version}/results/summary/salmon/{sample}"
    threads: 10
    conda:
        config['softwares']['salmon']['conda']
    output:
        tsv="{project}/{genome_version}/results/summary/salmon/{sample}/quant.sf",
    shell:
        """
        salmon quant -i {params.index} -l A -o {params.result_prefix} -1 {input.R1} -2 {input.R2} -p {threads}
        """