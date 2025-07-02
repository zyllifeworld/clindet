include:"CNV/ASCAT.smk"
### PoN of factesCH
include:"CNV/facets.smk"
include:"CNV/purple.smk"
include:"CNV/Battenberg.smk"
##### sequenza section
include:"CNV/sequenza.smk"

rule SM_check:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        tsv="{project}/{genome_version}/results/stats/paired/sample_check/{sample}/{sample}_check.tsv",
    params:
        wd="{project}/{genome_version}/results/stats/paired/sample_check/{sample}",
        rdata="{project}/{genome_version}/results/stats/paired/sample_check/{sample}/{sample}_check.rdata",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../scripts/sample_check.R"