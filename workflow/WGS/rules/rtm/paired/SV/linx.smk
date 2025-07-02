rule paired_linx:
    input:
        qc="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc",
        purple="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple"
    output:
        svs="{project}/{genome_version}/results/sv/paired/linx/{sample}/{sample}.linx.svs.tsv",
        output_dir=directory("{project}/{genome_version}/results/sv/paired/linx/{sample}"),
    params:
        vcf="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.sv.vcf.gz",
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['linx']['ensembl_data_dir'],
        driver_gene_panel=config['singularity']['hmftools'][genome_version]['linx']['driver_gene_panel'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['linx']['ref_genome_version'],
        known_fusion_file=config['singularity']['hmftools'][genome_version]['linx']['known_fusion_file'],
        sv_vcf=get_purple_sv_vcf
    threads: 10
    singularity:config['singularity']['hmftools']['sif']
    shell:
        """
        linx \
        -sample {wildcards.sample} \
        -ref_genome_version {params.ref_genome_version} \
        -sv_vcf {params.vcf} \
        -purple_dir {input.purple} \
        -threads {threads} \
        -output_dir {output.output_dir} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -known_fusion_file {params.known_fusion_file} \
        -driver_gene_panel {params.driver_gene_panel} \
        -log_debug -write_vis_data -write_all_vis_fusions
        """


rule report_linx:
    input:
        # output_dir=directory("{project}/{genome_version}/results/sv/paired/linx/{sample}"),
        svs="{project}/{genome_version}/results/sv/paired/linx/{sample}/{sample}.linx.svs.tsv",
    output:
        stamp="{project}/{genome_version}/results/sv/paired/linx/{sample}/plot.stamp"
    params:
        output_dir=directory("{project}/{genome_version}/results/sv/paired/linx/{sample}"),
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['linx']['ensembl_data_dir'],
        driver_gene_panel=config['singularity']['hmftools'][genome_version]['linx']['driver_gene_panel'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['linx']['ref_genome_version'],
        known_fusion_file=config['singularity']['hmftools'][genome_version]['linx']['known_fusion_file'],
        sv_vcf=get_purple_sv_vcf
    threads: 10
    singularity:config['singularity']['hmftools']['sif']
    shell:
        """
        java -cp /usr/bin/linx_v1.25.jar com.hartwig.hmftools.linx.visualiser.SvVisualiser \
        -sample {wildcards.sample} \
        -vis_file_dir {params.output_dir} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -circos /opt/circos-0.69-2/bin/circos \
        -plot_reportable \
        -debug 
        touch {output.stamp}
        """