import pandas as pd
sample_sheet_file = config['project']['sample_sheet']
samples_info = pd.read_csv(sample_sheet_file,index_col='Sample_name')
SE_samples = samples_info.loc[pd.isna(samples_info['R2_file_path'])].index.tolist()
PE_samples = samples_info.loc[~pd.isna(samples_info['R1_file_path'])].index.tolist()

## split config files 
configfile: "/public/ClinicalExam/lj_sih/projects/project_clindet/workflow/config/conf/genomes.yaml"
configfile: "/public/ClinicalExam/lj_sih/projects/project_clindet/workflow/config/conf/singularity.yaml"
configfile: "/public/ClinicalExam/lj_sih/projects/project_clindet/workflow/config/conf/softwares.yaml"

groups = ['NC','T']
stages = ['salmon','isofox','RSEM','arriba']#]'RSEM','arriba','salmon','TRUST4','salmon','kallisto','isofox']
rna_caller_list = ['sentieon','Mutect2_filter','vardict','varscan2','lofreq']#,'strelka']# vardict,varscan2,
project = config['project']['output_dir']
# genome_version = 'hg38'
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
    "{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
    "{project}/{genome_version}/results/summary/isofox/{sample}/{sample}.isf.fusions.csv",
    #### Case report #####

    #### mutation section #####
    ## mutect2 ##
    "{project}/{genome_version}/results/mut/vcf/{sample}/Mutect2_filter.vcf",
    ## vardict ##
    "{project}/{genome_version}/results/mut/vcf/{sample}/vardict.vcf",
    ## freebayes ##
    # "{project}/{genome_version}/results/mut/vcf/{sample}/freebayes.vcf",
    ## Strelka ##
    # "{project}/{genome_version}/results/mut/vcf/{sample}/strelka.vcf",
    ## Lofreq ##
    "{project}/{genome_version}/results/mut/vcf/{sample}/lofreq.vcf",
    ## sentieon ##
    "{project}/{genome_version}/results/mut/vcf/{sample}/sentieon.vcf",
    "{project}/{genome_version}/results/mut/vcf/{sample}/varscan2.vcf",
    ## ##
    "{project}/{genome_version}/results/mut/maf/{sample}/merge/{sample}.maf"
]
rna_res_list = list(filter(None, rna_res_list))
rule all:
    input:
        ## paired sample
        expand(rna_res_list,
        # sample = PE_samples,
        sample = ['CD1','COLO829'],
        project = project,
        genome_version = genome_version
        ),
        
##### Modules #####
include: "../workflow/RNA/Snakefile"
include: "../workflow/RNA/rules/isofox.smk"
