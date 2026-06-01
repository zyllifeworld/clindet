#### BRASS workflow
rule SV_brass_bamstat:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bas",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam.bas"
    singularity:
        flexible_container_img(config,['singularity','brass','sif'],image_url = config['singularity']['brass']['repo'])
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
        cnv=rules.paired_purple.output.pp if 'purple' in somatic_cnv_list else rules.CNA_ASCAT.output.rdata
    output:
        ascat="{project}/{genome_version}/results/cnv/paired/BRASS/{sample}/{sample}.ascat"
    params:
        brass_cnv='purple' if 'purple' in somatic_cnv_list else 'ascat'
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
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
        gc=config['softwares_params'][genome_version]['brass']['gc'],
        b=config['softwares_params'][genome_version]['brass']['b'],
        d=config['softwares_params'][genome_version]['brass']['d'],
        cb=config['softwares_params'][genome_version]['brass']['cb'],
        ct=config['softwares_params'][genome_version]['brass']['ct'],
        vi=config['softwares_params'][genome_version]['brass']['vi'],
        mi=config['softwares_params'][genome_version]['brass']['mi'],
        out_dir="{project}/{genome_version}/results/sv/paired/BRASS/{sample}"
    singularity:
        flexible_container_img(config,['singularity','brass','sif'],image_url = config['singularity']['brass']['repo'])
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
