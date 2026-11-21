rule STAR_1_pass:
    input:
        # R1="{project}/{genome_version}/results/trimmed/{sample}_R1.fastq.gz",
        # R2="{project}/{genome_version}/results/trimmed/{sample}_R2.fastq.gz",
        unpack(get_rna_fastq)
    output:
        pass1="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_pass1.log",
        sj_filt='{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/SJ.filt'
    params:
        out_dir="{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/", # "/"" must in the config string
        ref=config['resources'][genome_version]['REFFA'],
        star_index=config['softwares_params'][genome_version]['star']['index'],
        rg=r"ID:{sample} PL:ILLUMINA.NovaSeq LB:RNA-Seq SM:{sample}",
        read_files_in=get_star_read_files_in
    threads: 10
    benchmark:
        "{project}/{genome_version}/results/benchmarks/fusion/{sample}.star_arriba_map_1.benchmark.txt"
    conda:
        flexible_conda_env(config,['conda','rna'],env_yaml = 'envs/rsem.yaml')
    shell:
        """ 
        if [ ! -s {params.star_index}/SA ] || \
        [ ! -s {params.star_index}/Genome ] || \
        [ ! -s {params.star_index}/SAindex ]; then
            echo "[STAR index] STAR index not found or incomplete. Building STAR index..."
            mkdir -p {params.star_index}
            STAR \
                --runThreadN {threads} \
                --runMode genomeGenerate \
                --genomeDir {params.star_index} \
                --genomeFastaFiles {params.ref}
        else
            echo "[STAR index] Existing STAR index found. Skip building."
        fi
        STAR --genomeDir {params.star_index} --runThreadN={threads} \
            --outSAMtype None --outFileNamePrefix {params.out_dir} \
            --readFilesIn {params.read_files_in} --readFilesCommand zcat
        awk '$1~/(chr)?[1-2XY]/ && $6==0 && $5>0 && $7>0' {params.out_dir}SJ.out.tab > {output.sj_filt}
        touch {output.pass1}
        """



rule STAR_arriba_map:
    input:
        unpack(get_rna_fastq),
        sj_filt="{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/SJ.filt"
    output:
        bam=temp("{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam"),
        bai=temp("{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam.bai"),
        stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log",
        um_fq1="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_unmapped_R1.fq",
        um_fq2="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_unmapped_R2.fq",
    params:
        out_dir="{project}/{genome_version}/results/mapped/STAR/{sample}/",
        ref=lambda wildcards: config["resources"][wildcards.genome_version]["REFFA"],
        gtf=lambda wildcards: config["resources"][wildcards.genome_version]["GTF"],
        star_index=lambda wildcards: config["softwares_params"][wildcards.genome_version]["star"]["index"],
        rg=r"ID:{sample} PL:ILLUMINA.NovaSeq LB:RNA-Seq SM:{sample}",
        read_files_in=get_star_read_files_in,
        layout=get_rna_layout
    threads: 10
    benchmark:
        "{project}/{genome_version}/results/benchmarks/fusion/{sample}.star_arriba_map.benchmark.txt"
    conda:
        flexible_conda_env(config, ["conda", "rna"], env_yaml="envs/rsem.yaml")
    shell:
        """
        set -euo pipefail

        mkdir -p {params.out_dir}
        mkdir -p $(dirname {output.stamp})

        echo "[STAR_arriba_map] Sample: {wildcards.sample}" > {output.stamp}
        echo "[STAR_arriba_map] Layout: {params.layout}" >> {output.stamp}
        echo "[STAR_arriba_map] Reads: {params.read_files_in}" >> {output.stamp}

        STAR \
            --genomeDir {params.star_index} \
            --runThreadN {threads} \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.out_dir} \
            --outReadsUnmapped Fastx \
            --sjdbFileChrStartEnd {input.sj_filt} \
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
            --readFilesIn {params.read_files_in} \
            --readFilesCommand gunzip -c \
            >> {output.stamp} 2>&1

        samtools sort \
            -@ {threads} \
            -m 1G \
            -o {output.bam} \
            {params.out_dir}Aligned.out.bam \
            >> {output.stamp} 2>&1

        samtools index -@ {threads} {output.bam} >> {output.stamp} 2>&1

        if [ -s {params.out_dir}Unmapped.out.mate1 ]; then
            mv {params.out_dir}Unmapped.out.mate1 {output.um_fq1}
        else
            touch {output.um_fq1}
        fi

        if [ "{params.layout}" = "PE" ]; then
            if [ -s {params.out_dir}Unmapped.out.mate2 ]; then
                mv {params.out_dir}Unmapped.out.mate2 {output.um_fq2}
            else
                touch {output.um_fq2}
            fi
        else
            # Single-end RNA-seq has no mate2; create empty placeholder for downstream compatibility.
            touch {output.um_fq2}
        fi

        echo "[STAR_arriba_map] Finished." >> {output.stamp}
        """

# maybe sort by samtools,STAR always face RAM error
rule STAR_mut_map:
    input:
        unpack(get_rna_fastq)
    output:
        bam="{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam",
        unsort_bam="{project}/{genome_version}/results/mut/STAR/{sample}/Aligned.out.bam",
        bai="{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam.bai"
    params:
        out_dir="{project}/{genome_version}/results/mut/STAR/{sample}/",
        ref=config["resources"][genome_version]["REFFA"],
        gtf=config["resources"][genome_version]["GTF"],
        star_index=lambda wildcards: config["softwares_params"][genome_version]["star"]["index"],
        rg="'ID:{sample}' -r 'PL:ILLUMINA.NovaSeq' -r 'LB:RNA-Seq' -r 'SM:{sample}'",
        sort_mem_per_thread="1G",
        two_pass="--twopassMode Basic",
        read_files_in=get_star_read_files_in,
        layout=get_rna_layout
    threads: 10
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.star_mut_map.benchmark.txt"
    conda:
        flexible_conda_env(config, ["conda", "rna"], env_yaml="envs/rsem.yaml")
    shell:
        """
        if [ ! -s {params.star_index}/SA ] || \
           [ ! -s {params.star_index}/Genome ] || \
           [ ! -s {params.star_index}/SAindex ]; then
            echo "[STAR index] STAR index not found or incomplete. Building STAR index..."
            mkdir -p {params.star_index}
            STAR \
                --runThreadN {threads} \
                --runMode genomeGenerate \
                --genomeDir {params.star_index} \
                --genomeFastaFiles {params.ref}
        else
            echo "[STAR index] Existing STAR index found. Skip building."
        fi

        echo "[STAR_mut_map] Sample: {wildcards.sample}"
        echo "[STAR_mut_map] Layout: {params.layout}"
        echo "[STAR_mut_map] Reads: {params.read_files_in}"

        STAR \
            --genomeDir {params.star_index} \
            --runThreadN {threads} \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.out_dir} \
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
            --outSAMmapqUnique 60 \
            --readFilesIn {params.read_files_in} \
            --readFilesCommand gunzip -c

        samtools addreplacerg \
            -@ {threads} \
            -r {params.rg} \
            -o {params.out_dir}/addRG.bam \
            {params.out_dir}/Aligned.out.bam

        samtools sort \
            -@ {threads} \
            -m {params.sort_mem_per_thread} \
            -o {output.bam} \
            {params.out_dir}/addRG.bam

        samtools index -@ {threads} {output.bam}

        rm -f {params.out_dir}/addRG.bam
        """

rule cal_exp_RSEM:
    input:
        unpack(get_rna_fastq)
    output:
        genes="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.genes.results",
        isoforms="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.isoforms.results",
        bam=temp("{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.bam"),
        tx_bam=temp("{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.transcript.bam"),
        tx_sort_bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.transcript.sorted.bam",
    params:
        rsem_ref=lambda wildcards: config["softwares_params"][wildcards.genome_version]["rsem"]["index"],
        result_prefix="{project}/{genome_version}/results/summary/RSEM",
        paired_flag=get_rsem_paired_flag,
        reads=get_rsem_reads,
        layout=get_rna_layout
    threads: 10
    conda:
        flexible_conda_env(config, ["conda", "rna"], env_yaml="envs/rsem.yaml")
    benchmark:
        "{project}/{genome_version}/results/benchmarks/exp/{sample}.rsem.benchmark.txt"
    shell:
        """
        echo "[RSEM] Sample: {wildcards.sample}"
        echo "[RSEM] Layout: {params.layout}"
        echo "[RSEM] Reads: {params.reads}"

        rsem-calculate-expression \
            {params.paired_flag} \
            --star \
            -p {threads} \
            --star-output-genome-bam \
            --star-gzipped-read-file \
            --sort-bam-by-coordinate \
            {params.reads} \
            {params.rsem_ref} \
            {params.result_prefix}/{wildcards.sample}/{wildcards.sample}
        """

rule RSEM_sort_genome:
    input:
        bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.bam"
    threads: 10
    conda:
        flexible_conda_env(config,['conda','rna'],env_yaml = 'envs/rsem.yaml')
    output:
        bam_sort=temp("{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.sorted.bam"),
    shell:
        """
        samtools sort -@ {threads} -m 1G -o {output.bam_sort} {input.bam}
        samtools index {output.bam_sort}
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
        unpack(get_rna_fastq)
    output:
        tsv="{project}/{genome_version}/results/summary/kallisto/{sample}/abundance.tsv"
    params:
        index=config["softwares_params"][genome_version]["kallisto"]["index"],
        result_prefix="{project}/{genome_version}/results/summary/kallisto/{sample}",
        reads=get_kallisto_reads
    threads: 10
    conda:
        flexible_conda_env(config, ["conda", "rna"], env_yaml="envs/rsem.yaml")
    benchmark:
        "{project}/{genome_version}/results/benchmarks/exp/{sample}.kallisto.benchmark.txt"
    shell:
        """
        kallisto quant \
            -i {params.index} \
            -o {params.result_prefix} \
            -t {threads} \
            {params.reads}
        """
### samlom 
rule salmon:
    input:
        unpack(get_rna_fastq)
    output:
        tsv="{project}/{genome_version}/results/summary/salmon/{sample}/quant.sf"
    params:
        index=config["softwares_params"][genome_version]["salmon"]["index"],
        result_prefix="{project}/{genome_version}/results/summary/salmon/{sample}",
        reads=get_salmon_reads
    threads: 10
    conda:
        flexible_conda_env(config, ["conda", "rna"], env_yaml="envs/rsem.yaml")
    benchmark:
        "{project}/{genome_version}/results/benchmarks/exp/{sample}.salmon.benchmark.txt"
    shell:
        """
        salmon quant \
            -i {params.index} \
            -l A \
            -o {params.result_prefix} \
            -p {threads} \
            {params.reads}
        """