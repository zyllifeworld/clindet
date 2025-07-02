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
        "../../../scripts/ExomeDepth.R"

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

rule CNA_ABSOLUTE_GISTIC:
    input:
        cnv_rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata",
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}_filter.maf"
    output:
        rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
        seg="{project}/{genome_version}/results/cnv/GISTIC2/{sample}/{sample}.seg",
        ndt_seg="{project}/{genome_version}/results/clone/PhylogicNDT/{sample}/{sample}.seg.txt",
        ndt_maf="{project}/{genome_version}/results/clone/PhylogicNDT/{sample}/{sample}.maf",
        absolute_dir=directory("{project}/{genome_version}/results/cnv/ABSOLUTE/{sample}"),
        # absolute_pdf="{project}/{genome_version}/results/cnv/ABSOLUTE/{sample}/DoAbsolute.called.ABSOLUTE.plots.pdf"
    # conda:'snakemake'
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../scripts/ABSOLUTE.R"



rule CNA_Battenberg:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        # wd=directory("{project}/{genome_version}/results/cnv/Battenberg/{sample}"),
        sclo="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_copynumber.txt"
    params:
        # wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        wd="{project}/{genome_version}/results/cnv/Battenberg/{sample}",
        script=Path(str(workflow.current_basedir) + '/../../../scripts/subclone.R'),
        # mkdir -p {params.wd}
        # gender=,
        sample_index = lambda wildcards: wildcards.sample
    threads: 8
    shell:
        """
        Rscript {params.script} -t {wildcards.sample}  -n {wildcards.sample}_NC \
        --tb {input.Tum} \
        --nb {input.NC} --sex Male \
        -o {params.wd}
        """

include:'CNV/ASCAT.smk'
include:'CNV/freec.smk'
include:'CNV/cnv_facets.smk'
include:'CNV/sequenza.smk'
include:'CNV/purple.smk'