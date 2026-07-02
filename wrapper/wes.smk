import pandas as pd
sample_sheet_file = config['project']['sample_sheet']
samples_info = pd.read_csv(sample_sheet_file,index_col='Sample_name')

unpaired_samples = samples_info.loc[pd.isna(samples_info['Normal_R1_file_path'])].index.tolist()
paired_samples = samples_info.loc[~pd.isna(samples_info['Normal_R1_file_path'])].index.tolist()

## split config files 
configfile: "workflow/config/conf/genomes.yaml"
configfile: "workflow/config/conf/singularity.yaml"
configfile: "workflow/config/conf/softwares.yaml"
configfile: "workflow/config/conf/softwares_params.yaml"

project = config['project']['output_dir']
genome_version = config['project']['genome_version']
recal = config['project']['recal_BQSR']
recal_config = config['resources']['varanno'].get(genome_version, False)
if recal_config:
    recal = False
else:
    recal = recal

groups = ['NC','T']
## somatic mutation calling softwares
caller_list = config.get("run_params", {}).get("somatic_caller_list", [])
stages = config.get("run_params", {}).get("stages", [])
# germline mutation calling softwares
germ_caller_list = config.get("run_params", {}).get("germ_caller_list", [])
# somatic CNV calling softwares
somatic_cnv_list = config.get("run_params", {}).get("somatic_cnv_list", [])
# somatic SV calling softwares
somatic_sv_list = config.get("run_params", {}).get("somatic_sv_list", [])

### tumor only call
tumor_only_caller = config.get("run_params", {}).get("tumor_only_caller", [])
tumor_only_cnv_caller = config.get("run_params", {}).get("tumor_only_cnv_caller", [])

recall_pon =  False
custome_pon_db = True
pre_pon_db = True
recall_pon_pindel =  False
purple_sv = config.get("run_params", {}).get("purple_sv", 'gridss')


## paired sample list
paired_res_list = [
    ##### for QC report ######
    # rules.conpair_contamination.output           if 'conpair'          in stages else None,
    '{project}/{genome_version}/logs/paired/conpair/{sample}.done' if 'conpair'          in stages else None,

    ##### for SNV/INDEL calling #####
    "{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf" if 'call_mut' in stages else None,
    ## MAF need download VEP cache
    "{project}/{genome_version}/results/vcf_norm/paired/{sample}/merge/{sample}.vcf" if 'call_mut_vcf' in stages else None,

    # somatic_cnv_list = ['purple','ASCAT','facets','sequenza','freec','dryclean']
    ##### for CNV result ##### There is a bug for snakemake rules namelist when include *smk for 3-4 levels
    # rules.paired_purple.output.qc  if 'purple' in somatic_cnv_list else None, # purple call
    "{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc"  if 'purple' in somatic_cnv_list else None, # purple call
    # rules.CNA_ASCAT.output.rdata   if 'ASCAT'  in somatic_cnv_list else None, # ASCAT call
    "{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata"   if 'ASCAT'  in somatic_cnv_list else None, # ASCAT call

    "{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}_purity.png"         if 'facets' in somatic_cnv_list else None, # facets call
    "{project}/{genome_version}/results/cnv/paired/freec/{sample}/{sample}_config_freec.ini" if 'freec' in somatic_cnv_list else None, # Control-FREEC call
    "{project}/{genome_version}/results/cnv/paired/exomedepth/{sample}/{sample}_exomedepth.tsv"  if 'exomedepth' in somatic_cnv_list else None, # sequenza call
    # "{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}_segments.txt"  if 'sequenza' in somatic_cnv_list else None, # sequenza call
    #### Case report #####
    '{project}/{genome_version}/results/report/{sample}/{sample}_cancer_report.html' if 'case_report' in stages else None,
    "{project}/{genome_version}/results/report/{sample}/moalmanac/{sample}.report.html" if 'case_report' in stages else None,
    #### Multiple QC report #####
    '{project}/{genome_version}/results/multiqc_report.html' if 'multiqc' in stages else None,
    
]
paired_res_list = list(filter(None, paired_res_list))


## unpaired sample list
unpaired_res_list = [
    ##### only mapping #####
    "{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
    ##### for SNV/INDEL calling #####
    "{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}.maf"  if 'call_mut' in stages else None,
    "{project}/{genome_version}/results/vcf/unpaired/{sample}/merge/{sample}.vcf" if 'call_mut_vcf' in stages else None,
    ##### for CNV result ##### There is a bug for snakemake rules namelist when include *smk for 3-4 levels
    # rules.paired_purple.output.qc  if 'purple' in somatic_cnv_list else None, # purple call
    "{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/purple/{sample}.purple.qc"    if 'purple' in tumor_only_cnv_caller else None,
    ### if you want call CNV from use tumor-only WES data, take you own risk
    "{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}-T.bam_ratio.txt.png"  if 'freec' in tumor_only_cnv_caller else None,

]
unpaired_res_list = list(filter(None, unpaired_res_list))

##### Modules #####
rule wes_pipeline:
    input:
        ## paired sample
        expand(paired_res_list,
        sample = paired_samples,
        project = project,
        genome_version = genome_version,
        group = groups,
        caller = caller_list),
        #### unpaired sample
        expand(unpaired_res_list,
        project = project,
        genome_version = genome_version,
        sample = unpaired_samples,
        caller = caller_list),
        ##### multiqc report ########
        f'{project}/{genome_version}/results/multiqc/filelist.txt' if 'report' in stages else [],
        f'{project}/{genome_version}/results/multiqc_report.html' if 'report' in stages else [],


include: '../workflow/WES/Snakefile'