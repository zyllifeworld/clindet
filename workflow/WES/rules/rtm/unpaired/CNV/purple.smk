### purple
### amber first
### for tumor only mode, because bam sample add _T prefix, so need add a _T suffix to run
rule unpaired_amber:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        # indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        bed=get_sample_bed
    output:
        pcf="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber/{sample}.amber.baf.pcf",
        output_dir=directory("{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber")
    params:
        output_dir="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber",
        loci=config['singularity']['hmftools'][genome_version]['amber']['loci'],
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 10
    # singularity:config['singularity']['hmftools']['sif']
    conda:config['singularity']['hmftools']['conda']
    shell:
        """
        amber  -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -target_regions_bed {input.bed} \
        -output_dir {params.output_dir} \
        -threads {threads} \
        -loci {params.loci}
        """

rule unpaired_cobalt:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        # indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        bed=get_sample_bed
    output:
        pcf="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/cobalt/{sample}.cobalt.ratio.pcf",
        output_dir=directory("{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/cobalt")
    params:
        output_dir="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/cobalt",
        tumor_only_diploid_bed=config['singularity']['hmftools'][genome_version]['cobalt']['tumor_only_diploid_bed'],
        gc_profile=config['singularity']['hmftools'][genome_version]['cobalt']['gc_profile'],
    threads: 10
    resources:
        mem_mb=lambda wildcards, input: max(0.4 * input.size_files_mb[0], 1000) 
    # singularity:config['singularity']['hmftools']['sif']
    conda:config['singularity']['hmftools']['conda']
    shell:
        """
        cobalt -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -output_dir {params.output_dir} \
        -threads {threads} \
        -pcf_gamma 50 \
        -tumor_only_diploid_bed {params.tumor_only_diploid_bed} \
        -gc_profile {params.gc_profile}
        """

rule unpaired_purple:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        pcf="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/cobalt/{sample}.cobalt.ratio.pcf",
        # indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        amber="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber",
        cobalt="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/cobalt",
        sage_vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/sage/{sample}.sage.pave.vcf.gz",
        ref_genome=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        qc="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/purple/{sample}.purple.qc",
        output_dir=directory("{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/purple")
    params:
        output_dir="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber",
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