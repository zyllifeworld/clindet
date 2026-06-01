rule unpair_freec_config:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ini_template=config['singularity']['freec']['config_wes']
    output:
        config="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}_config_freec.ini",
    params:
        bed=get_sample_bed,
        chrFiles=lambda wildcards: get_config_value(config,['singularity', 'freec', wildcards.genome_version, 'chrFiles'],default=""),
        chrLenFile=lambda wildcards: get_config_value(config,['singularity', 'freec', wildcards.genome_version, 'chrLenFile'],default=""),
        snp_file=lambda wildcards: get_config_value(config,['singularity', 'freec', wildcards.genome_version, 'snp_file'],default=""),
        forceGCcontentNormalization=0,# 1 for WES,0 for WGS
        maxThreads='30',
        outputDir="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}",
        sambamba=lambda wildcards: get_config_value(config,['singularity', 'freec', wildcards.genome_version, 'sambamba'],default=""),
        ref=config['resources'][genome_version]['REFFA'],
    threads: 10
    conda:flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    script:
        "../../../../scripts/freec/config_freec.py"

rule freec_call_unpaired:
    input:
        config="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}_config_freec.ini",
    output:
        ratio="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}-T.bam_ratio.txt",
        info="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}-T.bam_info.txt",
        cnv="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}-T.bam_CNVs",
    params:
    threads: 30
    singularity:
        flexible_container_img(config,['singularity','freec','sif'],image_url = config['singularity']['svaba']['repo'])
    shell:
        """
        freec -conf {input.config}
        """

rule plot_freec_unpaired:
    input:
        ratio="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}-T.bam_ratio.txt",
        cnv="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}-T.bam_CNVs",
    output:
        png="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}-T.bam_ratio.txt.png"
    singularity:
        flexible_container_img(config,['singularity','freec','sif'],image_url = config['singularity']['freec']['repo'])
    shell:
        """
        cat /usr/local/bin/makeGraph2.0.R | R --slave --args {input.ratio} 
        """
