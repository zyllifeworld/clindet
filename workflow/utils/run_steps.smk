# def generate_output_file_list(stages, somatic_cnv_list, somatic_sv_list, tumor_only_cnv_caller):
#     """
#     Generate output file list based on analysis type and tools configuration
    
#     Args:
#         stages (list): List of analysis stages
#         somatic_cnv_list (list): List of CNV callers for somatic analysis
#         somatic_sv_list (list): List of SV callers for somatic analysis
#         tumor_only_cnv_caller (list): List of CNV callers for tumor-only analysis
#     Returns:
#         list: List of output file paths
#     """
#     ## setup tools to run， based on genome version
#     paired_res_list = [
#         ##### for QC report ######
#         # rules.conpair_contamination.output           if 'conpair'          in stages else None,
#         f'{project}/{genome_version}/logs/paired/conpair/{sample}.done' if 'conpair'          in stages else None,

#         ##### for SNV/INDEL calling #####
#         f"{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf",

#         ##### for CNV result ##### There is a bug for snakemake rules namelist when include *smk for 3-4 levels
#         # rules.paired_purple.output.qc  if 'purple' in somatic_cnv_list else None, # purple call
#         f"{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc"  if 'purple' in somatic_cnv_list else None, # purple call
#         # rules.CNA_ASCAT.output.rdata   if 'ASCAT'  in somatic_cnv_list else None, # ASCAT call
#         f"{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata"   if 'ASCAT'  in somatic_cnv_list else None, # ASCAT call
#         # rules.facets.output.qc         if 'facets' in somatic_cnv_list else None, # facets call
#         # rules.facets.output.qc         if 'facets' in somatic_cnv_list else None, # facets call
#         f"{project}/{genome_version}/results/cnv/paired/freec/{sample}/{sample}_config_freec.ini" if 'freec' in somatic_cnv_list else None, # Control-FREEC call
#         # rules.CNA_exomedepth.output.tsv       if 'exomedepth' in somatic_cnv_list else None, # sequenza call
#         f"{project}/{genome_version}/results/cnv/paired/exomedepth/{sample}/{sample}_exomedepth.tsv"  if 'exomedepth' in somatic_cnv_list else None, # sequenza call
#         # rules.sequenza_call.output.segment       if 'sequenza' in somatic_cnv_list else None, # sequenza call
#         f"{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}_segments.txt"  if 'sequenza' in somatic_cnv_list else None, # sequenza call
        
#         ##### for SV result #####
#         # somatic_sv_list = ['BRASS','delly','gridss','igcaller','linx','svaba','Manta']
#         # BRASS call
#         f"{project}/{genome_version}/results/sv/paired/BRASS/{sample}/{sample}_brass.log"  if 'BRASS' in somatic_sv_list else None, # purple call
#         # DELLY call
#         f"{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf"   if 'delly'  in somatic_sv_list else None, # ASCAT call
#         # gridss call
#         f"{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf.bgz" if 'gridss' in somatic_sv_list else None, # Control-FREEC call
#         # linx call
#         f"{project}/{genome_version}/results/sv/paired/linx/{sample}/{sample}.linx.svs.tsv"  if 'linx' in somatic_sv_list else None, # sequenza call
#         # svaba call
#         f"{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf"  if 'svaba' in somatic_sv_list else None, # sequenza call
#         # igcaller call
#         f"{project}/{genome_version}/results/sv/paired/igcaller/{sample}/{sample}-T_IgCaller/{sample}-T_output_filtered.tsv"  if 'igcaller' in somatic_sv_list else None, # sequenza call
#         # Manta call
#         f"{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/somaticSV.vcf.gz"  if 'Manta' in somatic_sv_list else None, # sequenza call

#         #### Case report #####
#         f'{project}/{genome_version}/results/report/{sample}/{sample}_cancer_report.html' if 'case_report' in stages else None,
#         #### Multiple QC report #####
#         f'{project}/{genome_version}/results/multiqc_report.html' if 'multiqc' in stages else None,

#     ]
#     paired_res_list = list(filter(None, paired_res_list))


#     ## unpaired sample list
#     unpaired_res_list = [
#         ##### for SNV/INDEL calling #####
#         f"{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}.maf",
#         ##### for CNV result ##### There is a bug for snakemake rules namelist when include *smk for 3-4 levels
#         # rules.paired_purple.output.qc  if 'purple' in somatic_cnv_list else None, # purple call
#         f"{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/purple/{sample}.purple.qc"    if 'purple' in tumor_only_cnv_caller else None,
#         ### if you want call CNV from use tumor-only WES data, take you own risk
#         f"{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}-T.bam_ratio.txt.png"  if 'freec' in tumor_only_cnv_caller else None,

#     ]
#     unpaired_res_list = list(filter(None, unpaired_res_list))

#     return paired_res_list,unpaired_res_list

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
    # rules.facets.output.qc         if 'facets' in somatic_cnv_list else None, # facets call
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
