rule CNA_ASCAT:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata"
    params:
        # ASCAT should change config file because allelCounter need chr prefix in hg19 version
        wd="{project}/{genome_version}/results/cnv/paired/ascat/{sample}",
        lociprefix=config['softwares']['ascat'][genome_version]['loci_1000'],
        allelesprefix=config['softwares']['ascat'][genome_version]['alleles_1000'],
        GCcontentfile=config['softwares']['ascat'][genome_version]['GCcontentfile'],
        replictimingfile=config['softwares']['ascat'][genome_version]['replictimingfile'],
    threads: 8
    conda: config['conda']['clindet_main']
    script:
        "../../../../scripts/ASCAT.R"

rule ASCAT_GISTIC:
    input:
        cnv_rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata"
    output:
        seg="{project}/{genome_version}/results/cnv/paired/GISTIC2_seg/{sample}/{sample}.seg"
    conda: config['conda']['clindet_main']
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../../scripts/ASCAT_to_GISTIC2.R"