rule unpair_CNA_ASCAT_sc:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
    output:
        rdata="{project}/{genome_version}/results/cnv/ASCATsc/unpaired/{sample}/{sample}_ASCATsc.rdata"
    params:
        wd="{project}/{genome_version}/results/cnv/ASCATsc/unpaired/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../scripts/ASCATsc.R"
