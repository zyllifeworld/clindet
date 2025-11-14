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
    singularity:config['singularity']['facets']['sif']
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
        png="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}_purity.cnv.png"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}",
        genome=config['singularity']['facets'][genome_version].get('genome','hg19')
    threads: 8
    singularity:config['singularity']['facets']['sif']
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
        --genome {params.genome}
        """