rule call_variants_sentieon:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        pon=config['resources'][genome_version]['WGS_PON']
    output:
        temp_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sent_temp.vcf",
        ori_data="{project}/{genome_version}/results/vcf/paired/{sample}/sentieon.ori_data"
        # contam_data="{project}/{genome_version}/results/vcf/paired/{sample}/sentieon.contam_data"
    params:
        temp_directory=config['params']['java']['temp_directory'],
        af_vcf=config['resources'][genome_version]['MUTECT2_VCF'],
        temp_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sent_temp.vcf",
    threads: 10
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.sentieon.benchmark.txt"
    shell:
        """
        {config[softwares][sentieon][call]}  driver -t {threads} -r {input.ref} \
        -i {input.Tum}  \
        -i {input.NC} \
        --algo TNhaplotyper2 --tumor_sample {wildcards.sample}_T \
        --normal_sample {wildcards.sample}_NC \
        --pon {input.pon} \
        {params.temp_vcf} \
        --algo OrientationBias --tumor_sample {wildcards.sample}_T \
        {output.ori_data}
        """


rule filter_sentieon:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        temp_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sent_temp.vcf",
        ori_data="{project}/{genome_version}/results/vcf/paired/{sample}/sentieon.ori_data"
    output:
        flag_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sentieon_flag.vcf",
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sentieon.vcf"
    threads: 10
    conda: flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        {config[softwares][sentieon][call]}  driver -t {threads} \
        -r  {input.ref}  \
        --algo TNfilter --tumor_sample {wildcards.sample}_T \
        --normal_sample {wildcards.sample}_NC \
        -v {input.temp_vcf}   \
        --orientation_priors {input.ori_data}  \
        {output.vcf}
        bcftools filter -i 'FILTER="PASS"'  {output.flag_vcf} > {output.vcf} 
        """