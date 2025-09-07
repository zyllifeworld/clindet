rule unpaird_deepvariant_call:
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
    singularity: config['singularity']['deepvariant_somatic']['sif']
    shell:
        """
        /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
        --model_type=FFPE_WES_TUMOR_ONLY \
        --ref={input.reference} \
        --reads_tumor={input.Tum} \
        --output_vcf={output.vcf} \
        --sample_name_tumor={wildcards.sample}_T \
        --num_shards={threads} \
        --logging_dir={params.out_prefix}/logs \
        --intermediate_results_dir={params.out_prefix} \
        --use_default_pon_filtering=true 
        """
## 1.9.0 version very low postive-result rate, consider remove
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