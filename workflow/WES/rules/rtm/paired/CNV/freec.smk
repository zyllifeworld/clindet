rule freec_config:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ini_template="workflow/WES/scripts/freec/config_exome.ini"
        #ascat_pre_dir= directory(config['resources'][genome_version]['ASCAT_WES_PreDIR'])
    output:
        config="{project}/{genome_version}/results/cnv/paired/freec/{sample}/{sample}_config_freec.ini",
    params:
        bed=get_sample_bed,
        chrFiles=config['singularity']['freec'][genome_version]['chrFiles'],
        chrLenFile=config['singularity']['freec'][genome_version]['chrLenFile'],
        forceGCcontentNormalization=1,#for WES
        maxThreads=30,
        outputDir="{project}/{genome_version}/results/cnv/paired/freec/{sample}",
    threads: 1
    script:
        "../../../../scripts/freec/config_freec.py"

rule freec_call_paired:
    input:
        config="{project}/{genome_version}/results/cnv/paired/freec/{sample}/{sample}_config_freec.ini",
    output:
        rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
    params:
    threads: 30
    singularity:config['singularity']['hmftools']['sif']
    shell:
        """
        freec -conf {input.config}
        """