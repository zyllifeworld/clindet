import pandas as pd

samples_info = pd.read_csv('/public/ClinicalExam/lj_sih/projects/project_clindet/resources/test_pipeline/test_fq/wes_test.csv',index_col='Sample_name')

unpaired_samples = samples_info.loc[pd.isna(samples_info['Normal_R1_file_path'])].index.tolist()
paired_samples = samples_info.loc[~pd.isna(samples_info['Normal_R1_file_path'])].index.tolist()


configfile: "/AbsoPath/of/clindet/folder/config/config.yaml"
project = samples_info["Project"].unique().tolist()[0]
genome_version = 'b37'

# stages you want run. conpair check  sample swap, case_report generate HTML case report, multiqc for QC report 
stages = ['conpair','case_report','multiqc']
#

## germline mutation calling softwares
germ_caller_list = ['strelkamanta','caveman']
## somatic mutation calling softwares
caller_list = ['HaplotypeCaller','strelkasomaticmanta','cgppindel_filter','caveman','muse','sage','Mutect2_filter','lofreq','vardict','verscan2','deepvariant']

## somatic CNV calling softwares
somatic_cnv_list = ['purple','ASCAT','facets','sequenza','freec','exomedepth']
# somatic SV calling softwares
# somatic_sv_list = ['BRASS','delly','gridss','svaba','Manta']


# tumor only mode
tumor_only_caller = ['sage']
tumor_only_cnv_caller = ['freec']

recall_pon =  False
custome_pon_db = True
recall_pon_pindel =  False

paired_res_list = [
    ##### for QC report ######
    # rules.conpair_contamination.output           if 'conpair'          in stages else None,
    '{project}/{genome_version}/logs/paired/conpair/{sample}.done' if 'conpair'          in stages else None,

    ##### for SNV/INDEL calling #####
    "{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf",

    ##### for CNV result ##### There is a bug for snakemake rules namelist when include *smk for 3-4 levels
    # rules.paired_purple.output.qc  if 'purple' in somatic_cnv_list else None, # purple call
    "{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc"  if 'purple' in somatic_cnv_list else None, # purple call
    # rules.CNA_ASCAT.output.rdata   if 'ASCAT'  in somatic_cnv_list else None, # ASCAT call
    "{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata"   if 'ASCAT'  in somatic_cnv_list else None, # ASCAT call
    # rules.facets.output.qc         if 'facets' in somatic_cnv_list else None, # facets call
    # rules.facets.output.qc         if 'freec' in somatic_cnv_list else None, # freec call
    "{project}/{genome_version}/results/cnv/paired/freec/{sample}/{sample}_config_freec.ini" if 'freec' in somatic_cnv_list else None, # Control-FREEC call
    # rules.CNA_exomedepth.output.tsv       if 'exomedepth' in somatic_cnv_list else None, # sequenza call
    "{project}/{genome_version}/results/cnv/paired/exomedepth/{sample}/{sample}_exomedepth.tsv"  if 'exomedepth' in somatic_cnv_list else None, # sequenza call
    # rules.sequenza_call.output.segment       if 'sequenza' in somatic_cnv_list else None, # sequenza call
    "{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}_segments.txt"  if 'sequenza' in somatic_cnv_list else None, # sequenza call
    
    ##### for SV result #####
    # somatic_sv_list = ['BRASS','delly','gridss','igcaller','linx','svaba','Manta']
    # BRASS call
    "{project}/{genome_version}/results/sv/paired/BRASS/{sample}/{sample}_brass.log"  if 'BRASS' in somatic_sv_list else None, # purple call
    # DELLY call
    "{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf"   if 'delly'  in somatic_sv_list else None, # ASCAT call
    # gridss call
    "{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf.bgz" if 'gridss' in somatic_sv_list else None, # Control-FREEC call
    # linx call
    "{project}/{genome_version}/results/sv/paired/linx/{sample}/{sample}.linx.svs.tsv"  if 'linx' in somatic_sv_list else None, # sequenza call
    # svaba call
    "{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf"  if 'svaba' in somatic_sv_list else None, # sequenza call
    # igcaller call
    "{project}/{genome_version}/results/sv/paired/igcaller/{sample}/{sample}-T_IgCaller/{sample}-T_output_filtered.tsv"  if 'igcaller' in somatic_sv_list else None, # sequenza call
    # Manta call
    "{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/somaticSV.vcf.gz"  if 'Manta' in somatic_sv_list else None, # sequenza call

    #### Case report #####
    '{project}/{genome_version}/results/report/{sample}/{sample}_cancer_report.html' if 'case_report' in stages else None,
    #### Multiple QC report #####
    '{project}/{genome_version}/results/multiqc_report.html' if 'multiqc' in stages else None,

]
paired_res_list = list(filter(None, paired_res_list))


## unpaired sample list
unpaired_res_list = [
    ##### for SNV/INDEL calling #####
    "{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}.maf",
    ##### for CNV result ##### There is a bug for snakemake rules namelist when include *smk for 3-4 levels
    # rules.paired_purple.output.qc  if 'purple' in somatic_cnv_list else None, # purple call
    "{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/purple/{sample}.purple.qc"    if 'purple' in tumor_only_cnv_caller else None,
    ### if you want call CNV from use tumor-only WES data, take you own risk
    "{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}-T.bam_ratio.txt.png"  if 'freec' in tumor_only_cnv_caller else None,

]
unpaired_res_list = list(filter(None, unpaired_res_list))

##### Modules #####
rule all:
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


include: '/public/ClinicalExam/lj_sih/projects/project_clindet/workflow/WES/Snakefile'


