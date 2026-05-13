rule facets_pileup:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        pileup="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}.snp_pileup.gz"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}",
        common_vcf=config['resources'][genome_version]['common_vcf'],
    threads: 8
    benchmark:
        "{project}/{genome_version}/results/benchmarks/cnv/{sample}.facets_pileup.benchmark.txt"
    singularity:
        flexible_container_img(config,['singularity','facets','sif'],image_url = config['singularity']['facets']['repo'])
    shell:
        """
        snp-pileup-wrapper.R \
        --vcf-file {params.common_vcf}  \
        --normal-bam {input.NC} \
        --tumor-bam {input.Tum} \
        --output-prefix {params.wd}
        """


rule facets_calling:
    input:
        pileup="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}.snp_pileup.gz"
    output:
        png="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}_purity.png",
        rds="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}_purity.rds",
        gene_level="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}.gene_level.txt",
        arm_level="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}.arm_level.txt"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/facets/{sample}",
        genome=lambda wildcards: get_config_value(config,['singularity', 'facets', wildcards.genome_version, 'genome'],default="hg19"),
    threads: 8
    singularity:
        flexible_container_img(config,['singularity','facets','sif'],image_url = config['singularity']['facets']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/cnv/{sample}.facets_call.benchmark.txt"
    shell:
        """
        run-facets-wrapper.R \
        --counts-file {input.pileup} \
        --sample-id {wildcards.sample} \
        --purity-cval 1000 --cval 500 \
        --everything \
        --directory {params.wd} \
        --genome {params.genome} -fl /usr/local/lib/R/site-library
        """

rule facets_gene_cna:
    input:
        gene_level="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}.gene_level.txt",
    output:
        cna="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}.cna.txt",
    params:
        wd="{project}/{genome_version}/results/cnv/paired/facets/{sample}",
        genome=lambda wildcards: get_config_value(config,['singularity', 'facets', wildcards.genome_version, 'genome'],default="hg19"),
    threads: 1
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    benchmark:
        "{project}/{genome_version}/results/benchmarks/cnv/{sample}.facets_call.benchmark.txt"
    script:
        "../../../../scripts/facet_gene_cna.R"

        
rule facets_annotate_maf:
    input:
        rds="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}_purity.rds",
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf"
    output:
        png="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}_purity.png"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/facets/{sample}",
        genome=lambda wildcards: get_config_value(config,['singularity', 'facets', wildcards.genome_version, 'genome'],default="hg19"),
    threads: 8
    singularity:
        flexible_container_img(config,['singularity','facets','sif'],image_url = config['singularity']['facets']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/cnv/{sample}.facets_annotate_maf.benchmark.txt"
    shell:
        """
        run-facets-wrapper.R \
        --counts-file {input.pileup} \
        --sample-id {wildcards.sample} \
        --purity-cval 1000 --cval 500 \
        --everything \
        --directory {params.wd} \
        --genome {params.genome} -fl /usr/local/lib/R/site-library
        """