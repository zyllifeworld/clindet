rule loop_vcf2maf_germ_paired:
    input:
        vcf='{project}/{genome_version}/results/vcf_germline/paired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf_germline/paired/{sample}/{caller}.vcf.maf"
    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        name=get_vcf_name,
        vep_data=config['softwares_params'][genome_version]['vcf2maf']['vep']['vep_data'],
        ncbi_build=config['softwares_params'][genome_version]['vcf2maf']['build_version'],
        cache_version=config['softwares_params'][genome_version]['vcf2maf']['vep']['cache_version'],
        species=config['softwares_params'][genome_version]['vcf2maf']['vep']['species']
    shell:
        """
       vcf2maf.pl --input-vcf {input.vcf} \
       --output-maf {output.maf} --ref-fasta {input.ref} \
       {params.name} \
       --vep-path $(realpath $(dirname $(which vep))) \
       --vep-data {params.vep_data} \
       --vep-fork 40 \
       --vep-overwrite \
       --species {params.species} \
       --ncbi-build {params.ncbi_build} \
       --cache-version {params.cache_version}
        """

rule loop_vcf2maf_germ_unpaired:
    input:
        vcf='{project}/{genome_version}/results/vcf_germline/unpaired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf_germline/unpaired/{sample}/{caller}.vcf.maf"
    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        name=get_vcf_name,
        vep_data=config['softwares_params'][genome_version]['vcf2maf']['vep']['vep_data'],
        ncbi_build=config['softwares_params'][genome_version]['vcf2maf']['build_version'],
        cache_version=config['softwares_params'][genome_version]['vcf2maf']['vep']['cache_version'],
        species=config['softwares_params'][genome_version]['vcf2maf']['vep']['species']
    shell:
        """
       {config[softwares][vcf2maf][call]} --input-vcf {input.vcf} \
       --output-maf {output.maf} --ref-fasta {input.ref} \
       {params.name} \
       --vep-path $(realpath $(dirname $(which vep))) \
       --vep-data {params.vep_data} \
       --vep-fork 40 \
       --vep-overwrite \
       --species {params.species} \
       --ncbi-build {params.ncbi_build} \
       --cache-version {params.cache_version}
        """


rule merge_paired_germ_maf:
    input:
        vcf1=expand("{project}/{genome_version}/results/vcf_germline/paired/{{sample}}/{caller}.vcf",caller = germ_caller_list,project = project,genome_version = genome_version),
        maf1=expand("{project}/{genome_version}/results/maf_germline/paired/{{sample}}/{caller}.vcf.maf",caller = germ_caller_list,project = project,genome_version = genome_version),
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf_germline/paired/{sample}/merge/{sample}.maf"
    params:
        dir="{project}/{genome_version}/results/maf_germline/paired/{sample}"
    script: merge_maf_script[config["run_type"]]


rule merge_unpaired_germ_maf:
    input:
        vcf1=expand("{project}/{genome_version}/results/vcf_germline/unpaired/{{sample}}/{caller}.vcf",caller = germ_caller_list,project = project,genome_version = genome_version),
        maf1=expand("{project}/{genome_version}/results/maf_germline/unpaired/{{sample}}/{caller}.vcf.maf",caller = germ_caller_list,project = project,genome_version = genome_version),
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf_germline/unpaired/{sample}/merge/{sample}.maf",
        filter_maf="{project}/{genome_version}/results/maf_germline/unpaired/{sample}/merge/{sample}_filter.maf"
    params:
        dir="{project}/{genome_version}/results/maf_germline/unpaired/{sample}"
    script: merge_maf_script[config["run_type"]]
