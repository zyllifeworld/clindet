rule loop_vcf2maf_paired:
    input:
        vcf='{project}/{genome_version}/results/vcf/paired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/paired/{sample}/{caller}.vcf.maf"
    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        name=get_vcf_name,
        ncbi_build=config['softwares']['vcf2maf']['build_version'][genome_version],
        vep_data=config['softwares']['vcf2maf']['vep'][genome_version]['vep_data'],
        vep_path=config['softwares']['vcf2maf']['vep'][genome_version]['vep_path'],
        species=config['softwares']['vcf2maf']['vep'][genome_version]['species'],
        cache_version=config['softwares']['vcf2maf']['vep'][genome_version]['cache_version']
    shell:
        """
        {config[softwares][vcf2maf][call]} --input-vcf {input.vcf} \
        --output-maf {output.maf} --ref-fasta {input.ref} \
        {params.name} \
        --vep-data {params.vep_data} \
        --vep-path {params.vep_path} \
        --vep-fork 40 \
        --vep-overwrite \
        --species {params.species} \
        --ncbi-build {params.ncbi_build} \
        --cache-version {params.cache_version}
        """

rule loop_vcf2maf_unpaired:
    input:
        vcf='{project}/{genome_version}/results/vcf/unpaired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/unpaired/{sample}/{caller}.vcf.maf"
    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        name=get_vcf_name,
        # vep_path=config['softwares']['vcf2maf']['vep'][genome_version]['vep_path'],
        vep_data=config['softwares']['vcf2maf']['vep'][genome_version]['vep_data'],
        ncbi_build=config['softwares']['vcf2maf']['build_version'][genome_version],
        cache_version=config['softwares']['vcf2maf']['vep'][genome_version]['cache_version'],
        species=config['softwares']['vcf2maf']['vep'][genome_version]['species']
    shell:
        """
       {config[softwares][vcf2maf][call]} --input-vcf {input.vcf} \
       --output-maf {output.maf} --ref-fasta {input.ref} \
       {params.name} \
       --vep-data {params.vep_data} \
       --vep-fork 40 \
       --vep-overwrite \
       --species {params.species}
       --ncbi-build {params.ncbi_build} \
       --cache-version {params.cache_version}
        """


rule merge_paired_maf:
    input:
        vcf1=expand("{project}/{genome_version}/results/vcf/paired/{{sample}}/{caller}.vcf",project = project,genome_version = genome_version,caller = caller_list),
        maf1=expand("{project}/{genome_version}/results/maf/paired/{{sample}}/{caller}.vcf.maf",project = project,genome_version = genome_version,caller = caller_list),
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf"
    params:
        dir="{project}/{genome_version}/results/maf/paired/{sample}"
    conda:
        "clindet"
    script:
        "../scripts/merge_maf.R"


rule merge_unpaired_maf:
    input:
        vcf1=expand("{project}/{genome_version}/results/vcf/unpaired/{{sample}}/{caller}.vcf",project = project,genome_version = genome_version,caller = caller_list),
        maf1=expand("{project}/{genome_version}/results/maf/unpaired/{{sample}}/{caller}.vcf.maf",project = project,genome_version = genome_version,caller = caller_list),
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}.maf"
    params:
        dir="{project}/{genome_version}/results/maf/unpaired/{sample}"
    conda:
        "clindet"
    script:
        "../scripts/merge_maf.R"
