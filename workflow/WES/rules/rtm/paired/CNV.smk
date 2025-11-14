# I will keep SM_chek,CNA_ABSOLUTE_GISTIC,CNA_Battenberg rule for future development
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
    conda: config['conda']['clindet_main']
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

ascat_config = config['softwares']['ascat_wgs'].get(genome_version, False)
if ascat_config:
    include:"CNV/ASCAT.smk"
### PoN of factesCH
include:"CNV/facets.smk"
### purple only work for b37 and hg38
if genome_version in ['b37','hg38']:
    include:"CNV/purple.smk"

##### freec section
include:'CNV/freec.smk'
##### exomedepth section
include:'CNV/exomedepth.smk'
##### cnv_facets section
include:'CNV/facets_suite.smk'
##### sequenza section
include:'CNV/sequenza.smk'
