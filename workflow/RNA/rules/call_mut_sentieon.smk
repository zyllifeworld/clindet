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
        {config[softwares][sentieon][call]}  driver -t {threads} \
        -i {input.bam} --algo LocusCollector --fun score_info {output.score}
        """

rule mark_duplicates:
    input:
        bam="{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam",
        score="{project}/{genome_version}/results/mut/dedup/{sample}.score.gz",
    output:
        bam=temp("{project}/{genome_version}/results/mut/dedup/{sample}.rmdep.bam"),
        # bai=temp("{project}/{genome_version}/results/dedup/{sample}.rmdep.bam.bai"),
        txt=temp("{project}/{genome_version}/results/mut/dedup/{sample}_metric.txt"),
    params:
        temp_directory=config['params']['java']['temp_directory'],
    threads:10
    shell:
        """
        {config[softwares][sentieon][call]}  driver -t {threads} \
        -i {input.bam} --algo Dedup --score_info {input.score} \
        --metrics {output.txt} {output.bam}
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
    threads:10
    shell:
        """ 
        {config[softwares][sentieon][call]}  driver -t {threads} -r {input.rsem_ref} \
        -i {input.bam} --algo RNASplitReadsAtJunction --reassign_mapq 255:60 {output.bam}
        """

rule sentieon_call:
    input:
        bam="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        #bai="{project}/{genome_version}/{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam.bai",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        vcf_raw="{project}/{genome_version}/results/mut/vcf/{sample}/sentieon_raw.vcf",
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/sentieon.vcf",
    params:
        temp_directory=config['params']['java']['temp_directory'],
        bed=config['resources'][genome_version]['WES_BED'],
    threads:10
    shell:
        """ 
        {config[softwares][sentieon][call]}  driver -t {threads} -r {input.ref} \
        -i {input.bam} --interval {params.bed} --algo Haplotyper --trim_soft_clip  \
        --call_conf 20 --emit_conf 20 {output.vcf_raw}
        bcftools filter -i 'FORMAT/DP[0] > 10' {output.vcf_raw} > {output.vcf}
        """


rule sentieon_anno_rna_edit:
    input:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/sentieon.vcf",
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/sentieon_anno_rnaedit.vcf",
    params:
        temp_directory=config['params']['java']['temp_directory'],
        bed=config['resources'][genome_version]['WES_BED'],
        rna_vcf=config['resources'][genome_version]['RNA_EDIT_VCF'],
    threads:10
    conda:'clindet'
    shell:
        """ 
        bgzip -k {input.vcf}
        tabix {input.vcf}.gz
        bcftools annotate \
            --output-type v \
            --annotations {params.rna_vcf} \
            --columns "INFO/RNAEDIT" \
            --output {output.vcf} \
            {input.vcf}.gz
        """
