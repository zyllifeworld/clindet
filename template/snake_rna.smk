import pandas as pd
samples_info = pd.read_csv('./pipe_RNA.csv',index_col='Sample_name')
SE_samples = samples_info.loc[pd.isna(samples_info['R2_file_path'])].index.tolist()
PE_samples = samples_info.loc[~pd.isna(samples_info['R1_file_path'])].index.tolist()

configfile: "/public/ClinicalExam/lj_sih/projects/project_clindet/build_log/config.yaml"

rna_stages = ['salmon','kallisto','RSEM','arriba','TRUST4']
rna_caller_list = ['sentieon','Mutect2_filter','vardict','varscan2','lofreq','freebayes']
project = samples_info["Project"].unique().tolist()[0]
genome_version = 'b37'

rna_res_list = [
    ##### for isoform expression RSEM ######
    "{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.genes.results" if 'RSEM'      in rna_stages else None,
    ##### kallisto
    "{project}/{genome_version}/results/summary/kallisto/{sample}/abundance.tsv"      if 'kallisto'  in rna_stages else None,
    ##### salmon
    "{project}/{genome_version}/results/summary/salmon/{sample}/quant.sf"             if 'salmon'    in rna_stages else None,
    ##### for Immu analysis #####
    "{project}/{genome_version}/results/IG/TRUST4/{sample}_report.tsv"                if 'TRUST4'    in rna_stages else None,
    ##### for fusion gene detection #####
    "{project}/{genome_version}/results/fusion/{sample}_arriba_fusion.tsv"            if 'arriba'    in rna_stages else None,
    ##### for isofox immu analysis #####
    "{project}/{genome_version}/results/summary/isofox/{sample}/{sample}.sorted.bam"  if 'isofox'    in rna_stages else None,

    #### mutation section #####
    "{project}/{genome_version}/results/mut/maf/{sample}/merge/{sample}.maf"
]
rna_res_list = list(filter(None, rna_res_list))
rule all:
    input:
        ## paired sample
        expand(rna_res_list,
        sample = PE_samples,
        project = project,
        genome_version = genome_version
        )
        
##### Modules #####
include: "/public/ClinicalExam/lj_sih/projects/project_clindet/workflow/RNA/Snakefile"
