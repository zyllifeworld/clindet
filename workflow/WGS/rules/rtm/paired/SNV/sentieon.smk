localrules: sentieon_build_pon
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
        {output.flag_vcf}
        bcftools filter -i 'FILTER="PASS"'  {output.flag_vcf} > {output.vcf} 
        """

rule sentieon_pon_call_normal:
    input:
        bam="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config["resources"][genome_version]["REFFA"]
    output:
        vcf="{project}/{genome_version}/results/pon/sentieon/vcf/{sample}.vcf.gz",
        tbi="{project}/{genome_version}/results/pon/sentieon/vcf/{sample}.vcf.gz.tbi"
    threads: 20
    shell:
        """
        {config[softwares][sentieon][call]} driver \
            -t {threads} \
            -r {input.ref} \
            -i {input.bam} \
            --algo TNscope \
            --tumor_sample {wildcards.sample}_NC \
            {output.vcf}

        tabix -f -p vcf {output.vcf}
        """

rule sentieon_build_pon:
    input:
        vcfs=expand("{project}/{genome_version}/results/pon/sentieon/vcf/{sample}.vcf.gz",project = project,genome_version = genome_version,sample=paired_samples)
    output:
        vcf="{project}/{genome_version}/results/pon/sentieon/sentieon_pon.vcf.gz",
        tbi="{project}/{genome_version}/results/pon/sentieon/sentieon_pon.vcf.gz.tbi"
    threads: 4
    conda:
        flexible_conda_env(config,["conda", "clindet_main"],env_yaml="envs/clindet.yaml")
    shell:
        """
        bcftools merge \
            --threads {threads} \
            -m all \
            -f PASS,. \
            --force-samples \
            {input.vcfs} \
        | bcftools +fill-AN-AC \
        | bcftools view -i 'SUM(AC)>1' \
        | bcftools sort \
            -Oz \
            -o {output.vcf}

        tabix -f -p vcf {output.vcf}
        """