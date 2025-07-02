rule CNA_Battenberg:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        # wd=directory("{project}/{genome_version}/results/cnv/Battenberg/{sample}"),
        sclo="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_copynumber.txt",
        pp="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_purity_ploidy.txt"
    params:
        # wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        wd="{project}/{genome_version}/results/cnv/Battenberg/{sample}",
        script=Path(str(workflow.current_basedir) + '/../../../../scripts/subclone.R'),
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

rule CNA_Battenberg_v2:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        # wd=directory("{project}/{genome_version}/results/cnv/Battenberg/{sample}"),
        sclo="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_subclone.txt",
        fitcnv="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_fitcnv.txt"
    params:
        # wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        wd="{project}/{genome_version}/results/cnv/Battenberg/{sample}",
        script=Path(str(workflow.current_basedir) + '/../../../../scripts/battenberg_wgs.R'),
        # mkdir -p {params.wd}
        # gender=,
        sample_index = lambda wildcards: wildcards.sample
    threads: 8
    shell:
        """
        Rscript --vanilla {params.script} \
        -a "paired" \
        -t {wildcards.sample}  -n {wildcards.sample}_NC \
        --tb {input.Tum} \
        --nb {input.NC} --sex Male \
        -o {params.wd} \
        --cpu {threads} \
        -g {wildcards.genome_version} \
        --fit_csv {output.fitcnv}
        """


rule CNA_Battenberg_combine:
    input:
        sclo="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_copynumber.txt"
    output:
        # wd=directory("{project}/{genome_version}/results/cnv/Battenberg/{sample}"),
        seg="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_average_seg.txt"
    params:
        # wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        cnv_cutoff=2000000,# 2Mb
        sample_name = lambda wildcards: wildcards.sample
    script:
        "../../../../scripts/Battenberg_merge.R"

rule CNA_Battenberg_ABSOLUTE_GISTIC:
    input:
        cnv_rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
        batt_cnv="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_average_seg.txt",
        batt_purity="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_purity_ploidy.txt",
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}_filter.maf"
    output:
        # rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
        seg="{project}/{genome_version}/results/cnv/GISTIC2/{sample}/{sample}.seg",
        ndt_seg="{project}/{genome_version}/results/clone/PhylogicNDT/{sample}/{sample}.seg.txt",
        ndt_maf="{project}/{genome_version}/results/clone/PhylogicNDT/{sample}/{sample}.maf",
        absolute_dir=directory("{project}/{genome_version}/results/cnv/ABSOLUTE/{sample}"),
        # absolute_pdf="{project}/{genome_version}/results/cnv/ABSOLUTE/{sample}/DoAbsolute.called.ABSOLUTE.plots.pdf"
    # conda:'snakemake'
    threads:8
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    script:
        "../../../../scripts/ABSOLUTE_Battenberg.R"
