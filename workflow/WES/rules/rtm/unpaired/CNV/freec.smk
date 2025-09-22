rule unpair_freec_config:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ini_template=config['singularity']['freec']['config_wes']
    output:
        config="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}/{sample}_config_freec.ini",
    params:
        bed=get_sample_bed,
        chrFiles=config['singularity']['freec'][genome_version]['chrFiles'],
        chrLenFile=config['singularity']['freec'][genome_version]['chrLenFile'],
        snp_file=config['singularity']['freec'][genome_version]['snp_file'],
        forceGCcontentNormalization=1,#for WES
        maxThreads='30',
        outputDir="{project}/{genome_version}/results/cnv/unpaired/freec/{sample}",
        sambamba=config['singularity']['freec'][genome_version]['sambamba'],
        ref=config['resources'][genome_version]['REFFA'],
    threads: 10
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
    singularity:config['singularity']['freec']['sif']
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
    singularity:config['singularity']['freec']['sif']
    shell:
        """
        cat /usr/local/bin/makeGraph2.0.R | R --slave --args {input.ratio} 
        """
