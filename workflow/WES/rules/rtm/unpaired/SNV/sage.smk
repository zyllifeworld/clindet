rule unpaired_sage:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        # indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        ref_genome=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/sage/{sample}.sage.vcf.gz",
    params:
        output_dir="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber",
        loci=config['singularity']['hmftools'][genome_version]['amber']['loci'],
        high_confidence_bed=config['singularity']['hmftools'][genome_version]['sage']['high_confidence_bed'],
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['sage']['ensembl_data_dir'],
        coverage_bed=config['singularity']['hmftools'][genome_version]['sage']['coverage_bed'],
        hotspots=config['singularity']['hmftools'][genome_version]['sage']['hotspots'],
        panel_bed=config['singularity']['hmftools'][genome_version]['sage']['panel_bed'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['sage']['ref_genome_version'],
    threads: 10
    # singularity:config['singularity']['hmftools']['sif']
    conda:config['singularity']['hmftools']['conda']
    shell:
        """
        sage \
        -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -threads {threads} \
        -skip_msi_jitter \
        -coverage_bed {params.coverage_bed} \
        -hotspots  {params.hotspots} \
        -panel_bed {params.panel_bed} \
        -high_confidence_bed {params.high_confidence_bed} \
        -output_vcf {output.vcf} -write_bqr_plot
        """

rule unpaired_sage_filter_pass:
    input:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/sage/{sample}.sage.vcf.gz",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/sage.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    threads: 1
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf} 
        """

rule unpaired_pave_anno_sage:
    input:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/sage/{sample}.sage.vcf.gz",
        ref_genome=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/sage/{sample}.sage.pave.vcf.gz",
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    conda:config['singularity']['hmftools']['conda']
    # singularity:"https://depot.galaxyproject.org/singularity/hmftools-pave:1.7.1--hdfd78af_0" # snakemake not work well for conda & singularity
    params:
        wd="{project}/{genome_version}/results/vcf/unpaired/{sample}/sage",
        driver_gene_panel=config['singularity']['hmftools'][genome_version]['purple']['driver_gene_panel'],
        high_confidence_bed=config['singularity']['hmftools'][genome_version]['sage']['high_confidence_bed'],
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['sage']['ensembl_data_dir'],
        coverage_bed=config['singularity']['hmftools'][genome_version]['sage']['coverage_bed'],
        hotspots=config['singularity']['hmftools'][genome_version]['sage']['hotspots'],
        panel_bed=config['singularity']['hmftools'][genome_version]['sage']['panel_bed'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['sage']['ref_genome_version'],
    threads: 8
    resources:
        mem_mb=lambda wildcards, input: max(100 * input.size_files_mb[0], 1000) # 100 times vcf file size
    shell:
        """
        pave -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        -sample {wildcards.sample} \
        -vcf_file {input.vcf} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -driver_gene_panel {params.driver_gene_panel} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -output_dir {params.wd} \
        -threads {threads}
        """

rule unpaired_pave_filter_pass:
    input:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/sage/{sample}.sage.pave.vcf.gz",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/pave.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    threads: 1
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf} 
        """