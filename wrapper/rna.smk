import pandas as pd
sample_sheet_file = config['project']['sample_sheet']
samples_info = pd.read_csv(sample_sheet_file,index_col='Sample_name')
SE_samples = samples_info.loc[pd.isna(samples_info['R2_file_path'])].index.tolist()
PE_samples = samples_info.loc[~pd.isna(samples_info['R1_file_path'])].index.tolist()

## split config files 
configfile: "workflow/config/conf/genomes.yaml"
configfile: "workflow/config/conf/singularity.yaml"
configfile: "workflow/config/conf/softwares.yaml"
configfile: "workflow/config/conf/softwares_params.yaml"

stages = config['run_params']['stages']
rna_caller_list = config['run_params']['rna_caller_list']
project = config['project']['output_dir']

genome_version = config['project']['genome_version']

rna_res_list = [
    ##### for isoform expression RSEM ######
    "{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.genes.results" if 'RSEM'      in stages else None,
    ##### kallisto
    "{project}/{genome_version}/results/summary/kallisto/{sample}/abundance.tsv"      if 'kallisto'  in stages else None,
    ##### salmon
    "{project}/{genome_version}/results/summary/salmon/{sample}/quant.sf"             if 'salmon'    in stages else None,
    ##### for Immu analysis #####
    "{project}/{genome_version}/results/IG/TRUST4/{sample}_report.tsv"                if 'TRUST4'    in stages else None,
    ##### for fusion gene detection #####
    "{project}/{genome_version}/results/fusion/{sample}_arriba_fusion.tsv"            if 'arriba'    in stages else None,
    ##### for TRUST4 immu analysis #####
    "{project}/{genome_version}/results/summary/isofox/{sample}/{sample}.sorted.bam"  if 'isofox'    in stages else None,
    # "{project}/{genome_version}/results/mut/STAR/{sample}/{sample}.sorted.bam",
    #### Case report #####

    # #### mutation section #####
    "{project}/{genome_version}/results/mut/maf/{sample}/merge/{sample}.maf" if 'call_mut' in stages else None,
]
rna_res_list = list(filter(None, rna_res_list))
rule rna_pipeline:
    input:
        ## paired sample
        expand(rna_res_list,
        sample = PE_samples,
        project = project,
        genome_version = genome_version
        ),
        
##### Modules #####
include: "../workflow/RNA/Snakefile"

