### for amber, only work for genome without chr prefix, it will not be add to the pipeline in hg19 version test
### purple
### amber first
rule paired_amber:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        indexes="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bai",
        bed=get_sample_bed
    output:
        pcf="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber/{sample}.amber.baf.pcf",
        output_dir=directory("{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber")
    params:
        output_dir="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber",
        loci=config['singularity']['hmftools'][genome_version]['amber']['loci'],
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 10
    conda:config['singularity']['hmftools']['conda']
    shell:
        """
        amber  -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -reference {wildcards.sample}_NC -reference_bam {input.NC} \
        -target_regions_bed {input.bed} \
        -output_dir {params.output_dir} \
        -threads {threads} \
        -loci {params.loci}
        """

rule paired_cobalt:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        indexes="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bai",
        bed=get_sample_bed
    output:
        output_dir=directory("{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt")
    params:
        output_dir="{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt",
        tumor_only_diploid_bed=config['singularity']['hmftools'][genome_version]['cobalt']['tumor_only_diploid_bed'],
        gc_profile=config['singularity']['hmftools'][genome_version]['cobalt']['gc_profile'],
    threads: 10
    # singularity:config['singularity']['hmftools']['sif']
    resources:
        mem_mb=lambda wildcards, input: max(0.4 * input.size_files_mb[0], 1000) 
    # singularity:config['singularity']['hmftools']['sif']
    conda:config['singularity']['hmftools']['conda']
    shell:
        """
        cobalt -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -reference {wildcards.sample}_NC -reference_bam {input.NC} \
        -output_dir {params.output_dir} \
        -threads {threads} \
        -pcf_gamma 50 \
        -gc_profile {params.gc_profile}
        """

rule paired_purple:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        indexes="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bai",
        amber="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber",
        cobalt="{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt",
        sage_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.pave.vcf.gz",
        ref_genome=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        purple_som_snv_vcf   = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.somatic.vcf.gz',
        purple_som_cnv       = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.cnv.somatic.tsv',
        purple_som_gene_cnv  = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.cnv.gene.tsv',
        purple_circos_png    = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.circos.png',
        purple_input_png     = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.input.png',
        purple_cn_png        = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.copynumber.png',
        purple_ma_png        = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.map.png',
        purple_purity_png    = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.purity.range.png',
        purple_segment_png   = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.segment.png',
        purple_clonality_png = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.somatic.clonality.png',
        purple_ploidy_png    = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.somatic.png',
        qc="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc",
        pp='{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.purity.tsv',
        output_dir=directory("{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple")
    params:
        output_dir="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple",
        tumor_only_diploid_bed=config['singularity']['hmftools'][genome_version]['purple']['tumor_only_diploid_bed'],
        gc_profile=config['singularity']['hmftools'][genome_version]['purple']['gc_profile'],
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['purple']['ensembl_data_dir'],
        somatic_hotspots=config['singularity']['hmftools'][genome_version]['purple']['somatic_hotspots'],
        driver_gene_panel=config['singularity']['hmftools'][genome_version]['purple']['driver_gene_panel'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['purple']['ref_genome_version'],
    threads: 10
    # singularity:config['singularity']['hmftools']['sif']
    conda:config['singularity']['hmftools']['conda']
    resources:
        mem_mb=lambda wildcards, input: max(0.45 * input.size_files_mb[0], 1000) 
    shell:
        """
        purple -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        -tumor {wildcards.sample} \
        -reference {wildcards.sample}_NC \
        -amber {input.amber} \
        -cobalt {input.cobalt} \
        -target_regions_bed {input.bed} \
        -gc_profile {params.gc_profile} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -threads {threads} \
        -somatic_vcf {input.sage_vcf} \
        -somatic_hotspots {params.somatic_hotspots} \
        -driver_gene_panel {params.driver_gene_panel} \
        -circos $(readlink -f $(which circos)) \
        -output_dir {output.output_dir}
        """
