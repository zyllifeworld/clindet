#### BRASS workflow
rule SV_brass_bamstat:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bas",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam.bas"
    singularity: config['singularity']['brass']['sif']
    threads:10
    benchmark:
        "{project}/{genome_version}/results/benchmarks/sv/{sample}.bamstat.benchmark.txt"
    shell:
        """
            bam_stats -@ {threads} -i {input.Tum} -o {output.Tum}
            bam_stats -@ {threads} -i {input.NC} -o {output.NC}
        """

rule brass_cnv:
    input:
        cnv=rules.CNA_ASCAT.output.rdata if 'purple' in somatic_cnv_list else rules.paired_purple.output.pp
    output:
        ascat="{project}/{genome_version}/results/cnv/paired/BRASS/{sample}/{sample}.ascat"
    params:
        brass_cnv='purple' if 'purple' in somatic_cnv_list else rules.paired_purple.output.output_dir
    script:
        "../../../../scripts/stats_ascat.R"


rule SV_brass:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        Tumbas="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bas",
        NCbas="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam.bas",
        ascat="{project}/{genome_version}/results/cnv/paired/BRASS/{sample}/{sample}.ascat"
    output:
        log="{project}/{genome_version}/results/sv/paired/BRASS/{sample}/{sample}_brass.log"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        gc=config['singularity']['brass'][genome_version]['gc'],
        b=config['singularity']['brass'][genome_version]['b'],
        d=config['singularity']['brass'][genome_version]['d'],
        cb=config['singularity']['brass'][genome_version]['cb'],
        ct=config['singularity']['brass'][genome_version]['ct'],
        vi=config['singularity']['brass'][genome_version]['vi'],
        mi=config['singularity']['brass'][genome_version]['mi'],
        out_dir="{project}/{genome_version}/results/sv/paired/BRASS/{sample}"
    singularity: config['singularity']['brass']['sif']
    threads:20
    benchmark:
        "{project}/{genome_version}/results/benchmarks/sv/{sample}.brass.benchmark.txt"
    shell:
        """
        mkdir -p {params.out_dir}
        brass.pl -s human -as {wildcards.genome_version} -pr WGS \
        -b {params.b} \
        -c {threads} -o {params.out_dir} \
        -d {params.d} -g {params.ref} \
        -gc {params.gc} -cb {params.cb} \
        -vi {params.vi} -ct {params.ct} -mi {params.mi}\
        -t {input.Tum} \
        -n {input.NC} -ss {input.ascat}
        touch {output.log}
        """
