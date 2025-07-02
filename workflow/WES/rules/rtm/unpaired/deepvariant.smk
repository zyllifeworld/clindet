rule deepvariant_call:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions=get_sample_bed,
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant/{sample}.deepvariant.vcf",
        gvcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant/{sample}.deepvariant.gvcf",
        # dir="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant,
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant",
    threads: 10
    singularity: config['singularity']['deepvariant']['sif']
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WES \
        --ref={input.reference} \
        --reads={input.Tum} \
        --regions {input.regions} \
        --output_vcf={output.vcf} \
        --output_gvcf={output.gvcf} \
        --intermediate_results_dir {params.out_prefix} \
        --num_shards={threads}
        """

rule deepvariant_norm:
    input:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant/{sample}.deepvariant.vcf",
        gvcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant/{sample}.deepvariant.gvcf",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant.vcf"
        # dir="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant,
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant",
    threads: 1
    shell:
        """
        ln -s $(realpath {input.vcf}) {output.vcf}
        """