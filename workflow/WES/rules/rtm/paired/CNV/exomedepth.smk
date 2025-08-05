rule CNA_exomedepth:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        # bed=get_sample_bed,
        bed='/public/ClinicalExam/lj_sih/resource/WES/HyperExomeV2_hg19/HyperExomeV2_capture_targets.hg19.bed',
        ref=config['resources'][genome_version]['REFFA']
    output:
        rds="{project}/{genome_version}/results/cnv/paired/exomedepth/{sample}/{sample}_exomedepth.rds",
        tsv="{project}/{genome_version}/results/cnv/paired/exomedepth/{sample}/{sample}_exomedepth.tsv"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/exomedepth/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../../scripts/ExomeDepth.R"