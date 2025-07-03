## case report by hmftools orange model
rule bamMetrics_tumor:
    input: 
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        # indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        ref_genome=config['resources'][genome_version]['REFFA']
    output:
        output_dir=directory("{project}/{genome_version}/results/qc/bamtools/{sample}"),
    threads: 16
    shell:
        """
        bamtools \
        -output_id {wildcards.sample} \
        -sample {wildcards.sample} \
        -bam_file {input.Tum} \
        -ref_genome {input.ref_genome} \
        -output_dir {output.output_dir} \
        -threads {threads}
        """

rule orange:
    input:
        t_metrics="{project}/{genome_version}/results/qc/bamtools/{sample}".
        linx="{project}/{genome_version}/results/sv/paired/linx/{sample}",
        purple="{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple",    
    output:
        output_dir=directory("{project}/{genome_version}/results/report/{sample}"),
    conda:config['singularity']['hmftools']['conda']
    singularity:"https://depot.galaxyproject.org/singularity/hmftools-pave:1.7.1--hdfd78af_0"
    params:
        sage="{project}/{genome_version}/results/vcf/paired/{sample}/sage",
        ref_genome_version=config['singularity']['hmftools'][genome_version]['sage']['ref_genome_version'],
        doid_json=config['singularity']['hmftools'][genome_version]['disease_ontology'],
        cohort_mapping_tsv=config['singularity']['hmftools'][genome_version]['cohort_mapping'],
        cohort_percentiles_tsv=config['singularity']['hmftools'][genome_version]['cohort_percentiles'],
        driver_gene_panel=config['singularity']['hmftools'][genome_version]['driver_gene_panel'],
        known_fusion_file=config['singularity']['hmftools'][genome_version]['known_fusion_data'],
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['ensembl_data_resources'],
        signatures_etiology_tsv=config['singularity']['hmftools'][genome_version]['sigs_etiology']
    threads: 8
    shell:
        """
        orange \
        -experiment_type "PANEL" \
        -tumor_sample_id {wildcards.sample} \
        -primary_tumor_doids "" \
        -ref_genome_version {params.ref_genome_version} \
        -output_dir {outpur.output_dir} \
        -doid_json  {params.doid_json} \
        -cohort_mapping_tsv  {params.cohort_mapping_tsv} \
        -cohort_percentiles_tsv  {params.cohort_percentiles_tsv} \
        -driver_gene_panel {params.driver_gene_panel} \
        -known_fusion_file {params.known_fusion_file} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -tumor_metrics_dir {input.t_metrics} \
        -sage_dir {params.sage} \
        -purple_dir {input.purple} \
        -purple_plot_dir {input.putple}/plot \
        -linx_dir {input.linx} \
        -linx_plot_dir {input.linx}/plot  \
        -signatures_etiology_tsv {params.signatures_etiology_tsv}
        """
