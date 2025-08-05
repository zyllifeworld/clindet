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
        # vep_path=config['softwares']['vcf2maf']['vep'][genome_version]['vep_path'],
        vep_data=config['softwares']['vcf2maf']['vep'][genome_version]['vep_data'],
        ncbi_build=config['softwares']['vcf2maf']['build_version'][genome_version],
        cache_version=config['softwares']['vcf2maf']['vep'][genome_version]['cache_version'],
        species=config['softwares']['vcf2maf']['vep'][genome_version]['species']
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
    script:
        "../scripts/merge_maf.R"


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
    script:
        "../scripts/merge_maf.R"
