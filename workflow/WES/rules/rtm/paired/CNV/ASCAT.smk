rule CNA_ASCAT:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata"
    params:
        # ASCAT should change config file because allelCounter need chr prefix in hg19 version
        wd="{project}/{genome_version}/results/cnv/paired/ascat/{sample}",
        lociprefix=config['softwares_params'][genome_version]['ascat']['loci_1000'],
        allelesprefix=config['softwares_params'][genome_version]['ascat']['alleles_1000'],
        GCcontentfile=config['softwares_params'][genome_version]['ascat']['GCcontentfile'],
        replictimingfile=config['softwares_params'][genome_version]['ascat']['replictimingfile'],
    threads: 8
    benchmark:
        "{project}/{genome_version}/results/benchmarks/cnv/{sample}.ascat.benchmark.txt"
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    script:
        "../../../../scripts/ASCAT.R"

rule ASCAT_EXTRACT_PURITYPLOIDY:
    input:
        rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata",
    output:
        tsv="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_purity.ploidy.tsv",
    params:
        # sample_index= lambda wildcards: wildcards.sample
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    script:
        "../../../../scripts/ascat_pp.R"


rule ASCAT_GISTIC:
    input:
        cnv_rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata"
    output:
        seg="{project}/{genome_version}/results/cnv/paired/GISTIC2_seg/{sample}/{sample}.seg"
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../../scripts/ASCAT_to_GISTIC2.R"