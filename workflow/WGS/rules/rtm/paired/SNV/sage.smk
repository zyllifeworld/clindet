rule paired_sage:
    input:
        # Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        # indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        ref_genome=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.vcf.gz",
    params:
        high_confidence_bed=config['singularity']['hmftools'][genome_version]['sage']['high_confidence_bed'],
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['sage']['ensembl_data_dir'],
        coverage_bed=config['singularity']['hmftools'][genome_version]['sage']['coverage_bed'],
        hotspots=config['singularity']['hmftools'][genome_version]['sage']['hotspots'],
        panel_bed=config['singularity']['hmftools'][genome_version]['sage']['panel_bed'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['sage']['ref_genome_version'],
    threads: 30
    conda:config['singularity']['hmftools']['conda']
    # singularity:config['singularity']['hmftools']['sif']
    shell:
        """
        sage \
        -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -reference {wildcards.sample}_NC -reference_bam {input.NC} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -skip_msi_jitter \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -threads {threads} \
        -coverage_bed {params.coverage_bed} \
        -hotspots  {params.hotspots} \
        -panel_bed {params.panel_bed} \
        -high_confidence_bed {params.high_confidence_bed} \
        -output_vcf {output.vcf} -write_bqr_plot
        """

rule pave_anno_sage:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.vcf.gz",
        ref_genome=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.pave.vcf.gz",
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    conda:config['singularity']['hmftools']['conda']
    # singularity:"https://depot.galaxyproject.org/singularity/hmftools-pave:1.7.1--hdfd78af_0"
    params:
        wd="{project}/{genome_version}/results/vcf/paired/{sample}/sage",
        driver_gene_panel=config['singularity']['hmftools'][genome_version]['purple']['driver_gene_panel'],
        high_confidence_bed=config['singularity']['hmftools'][genome_version]['sage']['high_confidence_bed'],
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['sage']['ensembl_data_dir'],
        coverage_bed=config['singularity']['hmftools'][genome_version]['sage']['coverage_bed'],
        hotspots=config['singularity']['hmftools'][genome_version]['sage']['hotspots'],
        panel_bed=config['singularity']['hmftools'][genome_version]['sage']['panel_bed'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['sage']['ref_genome_version'],
    threads: 8
    shell:
        """
        export _JAVA_OPTIONS="-Xmx30g" && pave \
        -sample {wildcards.sample} \
        -vcf_file {input.vcf} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -driver_gene_panel {params.driver_gene_panel} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -output_dir {params.wd} \
        -threads {threads}
        """


rule sage_filter_pass:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.vcf.gz",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage.vcf"
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    params:
        ref=config['resources'][genome_version]['REFFA']
    threads: 1
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf} 
        """


