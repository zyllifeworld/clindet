import os
rule run_cancer_report:
    input:
        # rmd_file="workflow/WES/rules/report/cancer_report.Rmd",
        # af_keygenes="workflow/WES/rules/report/somatic_panel.tsv",
        af_keygenes=Path(str(workflow.current_basedir) + '/./somatic_panel.tsv'),
        rmd_file=Path(str(workflow.current_basedir) + "/./main.Rmd"),
        ## somatic mut MAF ########
        somatic_mut_maf = "{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf" if 'call_mut' in stages else [],
        ## CNV results
        ascat_rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata",
        purple_res_dir='{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple',
        # somatic_snv_summary  = rules.somatic_snv_summary.output.json,
        # somatic_snv_vcf      = '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz',
        # somatic_sv_tsv       = lambda wc: rules.prep_sv_tsv.output[0]
        #                        if (batch_by_name[wc.batch].sv_vcf and 'structural' in stages) else [],
        # somatic_sv_vcf       = lambda wc: '{batch}/structural/{batch}-manta.vcf.gz'
        #                        if (batch_by_name[wc.batch].sv_vcf and 'structural' in stages) else [],
        #purple result
        purple_som_snv_vcf   = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.somatic.vcf.gz',
        purple_som_cnv       = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.cnv.somatic.tsv',
        purple_som_gene_cnv  = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.cnv.gene.tsv',
        purple_purity        = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.purity.tsv',
        purple_qc            = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc',
        purple_circos_png    = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.circos.png',
        purple_input_png     = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.input.png',
        purple_cn_png        = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.copynumber.png',
        purple_ma_png        = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.map.png',
        purple_purity_png    = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.purity.range.png',
        purple_segment_png   = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.segment.png',
        purple_clonality_png = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.somatic.clonality.png',
        purple_ploidy_png    = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/plot/{sample}.somatic.png',
        ## patient and sample info
        # patient_csv: '~/Documents/Project/project_clindet/source/clindet_report/report/test_data/patient_info.csv'
        ## Assay QC
        ngs_bit_sample_geneder = "{project}/{genome_version}/results/qc/conpair/qc/ngs_bit/{sample}-gender.tsv",
        conpair_contamination = "{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T_contamination.txt",
        conpair_concordance = "{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T_concordance.txt",
        ngs_bit_mapqc_tumor = "{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}_T_mapqc.txt",
        ngs_bit_mapqc_normal = "{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}_NC_mapqc.txt"
    params:
        tumor_name = lambda wc: wc.sample,
        genome_version = lambda wc: wc.genome_version,
        batch_name = lambda wc: wc.project,
        report_rmd  = lambda wc, input: abspath(input.rmd_file),
        af_keygenes  = lambda wc, input: abspath(input.af_keygenes),
        output_file         = lambda wc, output: join(os.getcwd(), output['report_html']),
        ## Assay QC
        ngs_bit_sample_geneder = lambda wc, input: os.path.basename(abspath(input.ngs_bit_sample_geneder)),
        conpair_contamination =  lambda wc, input: os.path.basename(abspath(input.conpair_contamination)),
        conpair_concordance = lambda wc, input: os.path.basename(abspath(input.conpair_concordance)),
        ngs_bit_mapqc_tumor = lambda wc, input: os.path.basename(abspath(input.ngs_bit_mapqc_tumor)),
        ngs_bit_mapqc_normal = lambda wc, input: os.path.basename(abspath(input.ngs_bit_mapqc_normal)),
        ## somatic mut MAF ########
        somatic_mut_maf = lambda wc, input: (
            os.path.basename(abspath(input.somatic_mut_maf))
            if getattr(input, "somatic_mut_maf", None)
            else None
        ),
        ### CNV results dir
        ascat_res_dir=lambda wc, input: os.path.basename(abspath(input.ascat_rdata)),
        purple_res_dir=lambda wc, input: abspath(input.purple_res_dir),
        ## purple results
        purple_som_snv_vcf  = lambda wc, input: abspath(input.purple_som_snv_vcf),
        purple_som_cnv      = lambda wc, input: abspath(input.purple_som_cnv),
        purple_som_gene_cnv = lambda wc, input: abspath(input.purple_som_gene_cnv),
        # purple_germ_cnv     = lambda wc, input: abspath(input.purple_germ_cnv),
        purple_purity       = lambda wc, input: abspath(input.purple_purity),
        purple_qc           = lambda wc, input: abspath(input.purple_qc),
        ## conda list
        # conda_list          = lambda wc, input: abspath(input.conda_list),
        # img_dir_abs         = lambda wc, output: abspath(output.img_dir),
    output:
        report_html = '{project}/{genome_version}/results/report/{sample}/{sample}_cancer_report.html',
        # img_dir = directory('{project}/{genome_version}/results/reprot/{sample}/img'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000
    message:
        lambda wildcards, params: (
            f"Parameters for {wildcards.sample}\n"
            + format_params(params)
        )
    conda:'cancer_report'
    script:
        '../../scripts/rmd.R'