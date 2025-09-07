# maybe sort by samtools,STAR always face RAM error
rule STAR_isofox_map:
    input:
        unpack(get_rna_fastq),
    output:
        bam="{project}/{genome_version}/results/summary/isofox/{sample}/{sample}.sorted.bam",
        unsort_bam="{project}/{genome_version}/results/summary/isofox/{sample}/Aligned.out.bam",
    params:
        out_dir="{project}/{genome_version}/results/summary/isofox/{sample}/", 
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
        star_index=config['softwares']['star']['index'][genome_version],
        rg="'ID:{sample}' -r 'PL:ILLUMINA.NovaSeq' -r 'LB:RNA-Seq' -r 'SM:{sample}'",
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
            --outSAMunmapped Within \
            --outBAMcompression 0 --outSAMattributes All --outFilterMultimapNmax 10 \
            --outFilterMismatchNmax 3 limitOutSJcollapsed 3000000 -chimSegmentMin 10 \
            --chimOutType WithinBAM SoftClip \
            --chimJunctionOverhangMin 10 --chimSegmentReadGapMax 3 --chimScoreMin 1 \
            --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 --outFilterScoreMinOverLread 0.33 \
            --outFilterMatchNminOverLread 0.33 --outFilterMatchNmin 35 \
            --alignSplicedMateMapLminOverLmate 0.33 --alignSplicedMateMapLmin 35 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --readFilesIn {input.R1} {input.R2} --readFilesCommand gunzip -c
        samtools addreplacerg -@ {threads} -r {params.rg} -o {params.out_dir}/addRG.bam {params.out_dir}/Aligned.out.bam
        samtools sort -@ {threads} -m 1G -o {output.bam} {params.out_dir}/addRG.bam
        samtools index {output.bam}
        rm {params.out_dir}/addRG.bam
        """


rule isofox_call:
    input:
        bam="{project}/{genome_version}/results/summary/isofox/{sample}/{sample}.sorted.bam",
        ref_genome=config['resources'][genome_version]['REFFA']
    output:
        fusion="{project}/{genome_version}/results/summary/isofox/{sample}/{sample}.isf.fusions.csv",
        gene="{project}/{genome_version}/results/summary/isofox/{sample}/{sample}.isf.gene_data.csv"
    params:
        out_dir="{project}/{genome_version}/results/summary/isofox/{sample}/",
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['purple']['ensembl_data_dir'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['sage']['ref_genome_version'],
        # xms=lambda wildcards, input: max(1.5 * input.size_mb, 300),
        # xmx=lambda wildcards, input: max(1.5 * input.size_mb, 300),
    threads: 10
    resources:
        # Allocate memory based on input file size, with a minimum of 300 MB
        mem_mb=lambda wildcards, input: max(1.5 * input.size_mb, 300),
    conda:
        config['softwares']['star']['conda']
    shell:
        """ 
        isofox -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        -sample {wildcards.sample} \
        -functions "TRANSCRIPT_COUNTS;ALT_SPLICE_JUNCTIONS;FUSIONS" \
        -ref_genome_version {params.ref_genome_version} \
        -bam_file {input.bam} \
        -ref_genome {input.ref_genome} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -output_dir {params.out_dir} \
        -threads {threads}
        """