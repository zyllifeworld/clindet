import pandas as pd
samples_info = pd.read_csv('./pipe_WGS.csv',index_col='Sample_name')

unpaired_samples = samples_info.loc[pd.isna(samples_info['Normal_R1_file_path'])].index.tolist()
paired_samples = samples_info.loc[~pd.isna(samples_info['Normal_R1_file_path'])].index.tolist()

## split config files 
configfile: "/public/ClinicalExam/lj_sih/projects/project_clindet/workflow/config/conf/genomes.yaml"
configfile: "/public/ClinicalExam/lj_sih/projects/project_clindet/workflow/config/conf/singularity.yaml"
configfile: "/public/ClinicalExam/lj_sih/projects/project_clindet/workflow/config/conf/softwares.yaml"

project = config['project']['output_dir']
genome_version = config['project']['genome_version']

import os
# checking if the directory demo_folder 
# exist or not.
if not os.path.exists("logs/slurm"):
    os.makedirs("logs/slurm")

pre_pon_db = True


groups = ['NC','T']
germ_caller_list = ['caveman','deepvariant']#'vardict_germline',
somatic_caller_list = ['strelkasomaticmanta','muse','cgppindel_filter','sentieon','deepvariant','sage']#,'Mutect2_filter','cgppindel_filter'] # ,'Mutect2_filter','vardict_filter','vardict_filter','cgppindel_filter','vardict','strelkasomatic'pindel freebayes too slow 'UnifiedGenoTyper'
somatic_cnv_list = ['purple','ascat']
purple_sv = 'svaba'

recall_pon =  False
recall_pon_pindel =  False
recal = False
rule all:
    input:
        ## paired sample
        expand([
            "{project}/{genome_version}/results/recal/paired/{sample}-{group}.bam",
            # "{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber/{sample}.amber.baf.pcf",
            # "{project}/{genome_version}/results/cnv/paired/frag/tumor/{sample}/cov.rds",
            # "{project}/{genome_version}/results/cnv/paired/frag/normal/{sample}/cov.rds",
            # "{project}/{genome_version}/results/vcf/paired/{sample}/{sample}.sage.vcf.gz",
            # "{project}/{genome_version}/results/sv/paired/jabba/COL0829/jabba.rds",
            # "{project}/{genome_version}/results/sv/paired/het/COL0829/sites.txt",
            # "{project}/{genome_version}/results/cnv/Battenberg/COL0829/COL0829_fitcnv.txt",
            # "{project}/{genome_version}/results/virus/virusbreaked/COL0829/COL0829.vcf",
            # "{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf",
            # "{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf.bgz",
            # "{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc",
            # "{project}/{genome_version}/results/sv/paired/linx/{sample}/plot.stamp",
            # "{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf",
            # "{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_germ.vcf",
            # # "{project}/{genome_version}/results/sv/paired/esvee/{sample}/SV_gridss_{sample}.vcf",
            # "{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf",
            "{project}/{genome_version}/results/cnv/paired/dryclean/{sample}/drycleaned.cov.rds",
            "{project}/{genome_version}/results/cnv/paired/ascat/COL0829/COL0829_ASCAT.rdata",
            "{project}/{genome_version}/results/cnv/paired/ascat/HG008/HG008_ASCAT.rdata",
            "{project}/{genome_version}/results/sv/paired/BRASS/COL0829/COL0829_brass.log",
            "{project}/{genome_version}/results/sv/paired/svaba/COL0829/COL0829.svaba.somatic.indel.maf",
            "{project}/{genome_version}/results/sv/paired/svaba/COL0829/COL0829.svaba.somatic.sv.maf",
            # "{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf",
            # "{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}_filter.vcf",
            # "{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}_mutationtime.pdf",
            # "{project}/{genome_version}/results/sv/paired/gridss/{sample}/SV_gridss_{sample}.vcf"
        ],
        project = project,
        genome_version = genome_version,
        sample = ['COL0829','HG008'],
        # sample = ['MM-094','MM-099'],
	    group = groups,
        caller = somatic_caller_list),
        # 'logs/paired/Mutect2_PoNDB_MM.log',
        # "{project}/{genome_version}/results/vcf/paired/PoN/MM_pon.vcf.gz",
        # '{project}/{genome_version}/logs/paired/Mutect2_PoNVCF_MM.log',
        # 'analysis/normalPanel/pindel_MM.gff3.gz',
        ## unpaired sample
        expand([
            # "{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
            # "{project}/{genome_version}/results/stats/unpaired/wgs_metrics/{sample}-{group}.txt"
            # "{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}.maf",
            # "{project}/{genome_version}/results/maf/unpaired/{sample}/{caller}.vcf.maf",
        ],
        project = project,
        group = ['T'],
        genome_version = genome_version,
        sample = unpaired_samples,
        caller = somatic_caller_list)

include: "workflow/WGS/Snakefile"
include: "workflow/WGS/rules/mutationTime.smk"
include: "workflow/WGS/rules/rtm/paired/SV/jabba.smk"
include: "workflow/WGS/rules/rtm/paired/SV/linx.smk"
include: "workflow/WGS/rules/rtm/paired/SV/esvee.smk"
include: "workflow/WGS/rules/rtm/paired/VirusScan.smk"
include: '/public/ClinicalExam/lj_sih/projects/project_pipeline/WGS/BALL_WGS/anno_sv.smk'