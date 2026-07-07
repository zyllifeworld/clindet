rule paired_amber:
    input:
        # Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        # indexes="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bai",
    output:
        pcf="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber/{sample}.amber.baf.pcf",
        output_dir=directory("{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber")
    params:
        output_dir="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber",
        loci=config['softwares_params'][genome_version]['hmftools']['amber']['loci'],
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 10
    # singularity:config['singularity']['hmftools']['sif']
    conda: flexible_conda_env(config,['conda','hmftools'],env_yaml = 'envs/hmftools.yaml')
    shell:
        """
        amber  -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -reference {wildcards.sample}_NC -reference_bam {input.NC} \
        -output_dir {params.output_dir} \
        -threads {threads} \
        -loci {params.loci}
        """

rule paired_cobalt:
    input:
        # Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        # indexes="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bai",
    output:
        pcf="{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt/{sample}.cobalt.ratio.pcf",
        output_dir=directory("{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt")
    params:
        output_dir="{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt",
        gc_profile=config['softwares_params'][genome_version]['hmftools']['cobalt']['gc_profile'],
    threads: 10
    resources:
        mem_mb=lambda wildcards, input: max(0.4 * input.size_files_mb[0], 1000) 
    # singularity:config['singularity']['hmftools']['sif']
    conda: flexible_conda_env(config,['conda','hmftools'],env_yaml = 'envs/hmftools.yaml')
    shell:
        """
        cobalt -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -reference {wildcards.sample}_NC -reference_bam {input.NC} \
        -output_dir {params.output_dir} \
        -threads {threads} \
        -gc_profile {params.gc_profile}
        """

purple_run_with_sv = True
if purple_run_with_sv:
    rule paired_purple:
        input:
            # Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
            # Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
            # NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
            pcf="{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt/{sample}.cobalt.ratio.pcf",
            # indexes="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bai",
            amber="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber",
            cobalt="{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt",
            sage_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.pave.vcf.gz",#use pave annotate vcf
            sv_vcf=purple_sv_vcf,
            ref_genome=config['resources'][genome_version]['REFFA'],
        output:
            qc="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc",
            output_dir=directory("{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple"),
            pp="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.purity.tsv"
        params:
            output_dir="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber",
            tumor_only_diploid_bed=config['softwares_params'][genome_version]['hmftools']['purple']['tumor_only_diploid_bed'],
            gc_profile=config['softwares_params'][genome_version]['hmftools']['purple']['gc_profile'],
            ensembl_data_dir=config['softwares_params'][genome_version]['hmftools']['purple']['ensembl_data_dir'],
            somatic_hotspots=config['softwares_params'][genome_version]['hmftools']['purple']['somatic_hotspots'],
            driver_gene_panel=config['softwares_params'][genome_version]['hmftools']['purple']['driver_gene_panel'],
            ref_genome_version=config['softwares_params'][genome_version]['hmftools']['purple']['ref_genome_version'],
            sv_vcf=get_purple_sv_vcf
        threads: 10
        resources:
            mem_mb=lambda wildcards, input: max(0.45 * input.size_files_mb[0], 1000) 
        # singularity:config['singularity']['hmftools']['sif']
        conda: flexible_conda_env(config,['conda','hmftools'],env_yaml = 'envs/hmftools.yaml')
        shell:
            """
            purple -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
            -reference {wildcards.sample}_NC \
            -tumor {wildcards.sample} \
            -amber {input.amber} \
            -cobalt {input.cobalt} \
            -gc_profile {params.gc_profile} \
            -ref_genome_version {params.ref_genome_version} \
            -ref_genome {input.ref_genome} \
            -ensembl_data_dir {params.ensembl_data_dir} \
            -somatic_vcf {input.sage_vcf} \
            {params.sv_vcf} \
            -threads {threads} \
            -somatic_hotspots {params.somatic_hotspots} \
            -driver_gene_panel {params.driver_gene_panel} \
            -circos /opt/circos-0.69-2/bin/circos \
            -output_dir {output.output_dir}
            """
else:
    rule paired_purple:
        input:
            # Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
            # Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
            # NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
            pcf="{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt/{sample}.cobalt.ratio.pcf",
            # indexes="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bai",
            amber="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber",
            cobalt="{project}/{genome_version}/results/cnv/paired/purple/{sample}/cobalt",
            sage_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}.sage.vcf.gz",
            ref_genome=config['resources'][genome_version]['REFFA'],
        output:
            qc="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc",
            output_dir=directory("{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple"),
            pp="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.purity.tsv"
        params:
            output_dir="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber",
            tumor_only_diploid_bed=config['softwares_params'][genome_version]['hmftools']['purple']['tumor_only_diploid_bed'],
            gc_profile=config['softwares_params'][genome_version]['hmftools']['purple']['gc_profile'],
            ensembl_data_dir=config['softwares_params'][genome_version]['hmftools']['purple']['ensembl_data_dir'],
            somatic_hotspots=config['softwares_params'][genome_version]['hmftools']['purple']['somatic_hotspots'],
            driver_gene_panel=config['softwares_params'][genome_version]['hmftools']['purple']['driver_gene_panel'],
            ref_genome_version=config['softwares_params'][genome_version]['hmftools']['purple']['ref_genome_version'],
            sv_vcf=get_purple_sv_vcf
        threads: 10
        # singularity:config['singularity']['hmftools']['sif']
        conda: flexible_conda_env(config,['conda','hmftools'],env_yaml = 'envs/hmftools.yaml')
        resources:
            mem_mb=lambda wildcards, input: max(0.45 * input.size_files_mb[0], 1000) 
        shell:
            """
            purple -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
            -reference {wildcards.sample}_NC \
            -tumor {wildcards.sample} \
            -amber {input.amber} \
            -cobalt {input.cobalt} \
            -gc_profile {params.gc_profile} \
            -ref_genome_version {params.ref_genome_version} \
            -ref_genome {input.ref_genome} \
            -ensembl_data_dir {params.ensembl_data_dir} \
            -somatic_vcf {input.sage_vcf} \
            {params.sv_vcf} \
            -threads {threads} \
            -somatic_hotspots {params.somatic_hotspots} \
            -driver_gene_panel {params.driver_gene_panel} \
            -circos $(readlink -f $(which circos)) \
            -output_dir {output.output_dir}
            """

