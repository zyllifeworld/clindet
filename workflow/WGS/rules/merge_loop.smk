# rule check_vcf:
#     input:
#         paired=expand("{project}/{genome_version}/results/vcf/unpaired/{{sample}}/{caller}.vcf",caller = caller_list),#sample=unpaired_samples
#         unpaired=expand("{project}/{genome_version}/results/vcf/paired/{{sample}}/{caller}.vcf",caller = caller_list),#,sample = paired_samples
#     output:
#         check="logs/check/{sample}/{caller}.log"
#     shell:
#         """
#         touch {output.check}
#         """

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

# rule loop_annovar_paired:
#     input:
#         vcf='{project}/{genome_version}/results/vcf/paired/{sample}/{caller}.vcf',
#         ref=config['resources']['hg19']['REFFA']
#     output:
#         maf="{project}/{genome_version}/results/maf/paired/{sample}/{caller}.vcf.maf"
#     conda:
#         config['softwares']['vcf2maf']['conda']
#     params:
#         name=get_vcf_name
#     shell:
#         """
#         {config[softwares][vcf2maf][call]} --input-vcf {input.vcf} \
#         --output-maf {output.maf} --ref-fasta {input.ref} \
#         {params.name} \
#         --vep-data /public/home/lijf/env/db/vep \
#         --vep-path /public/home/lijf/env/miniconda3/envs/vep105/bin \
#         --vep-fork 40 \
#         --vep-overwrite \
#         --ncbi-build GRCh37 \
#         --cache-version 105
#         """

rule loop_vcf2maf_unpaired:
    input:
        vcf='{project}/{genome_version}/results/vcf/unpaired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/unpaired/{sample}/{caller}.vcf.maf"
    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        name=get_vcf_name
    shell:
        """
        {config[softwares][vcf2maf][call]} --input-vcf {input.vcf} \
        --output-maf {output.maf} --ref-fasta {input.ref} \
        {params.name} \
        --vep-data /public/home/lijf/env/db/vep \
        --vep-path /public/home/lijf/env/miniconda3/envs/vep105/bin \
        --vep-fork 40 \
        --vep-overwrite \
        --ncbi-build GRCh37 \
        --cache-version 105
        """


rule merge_paired_maf:
    input:
        vcf1=expand("{project}/{genome_version}/results/vcf/paired/{{sample}}/{caller}.vcf",project = project,genome_version = genome_version,caller = caller_list),
        maf1=expand("{project}/{genome_version}/results/maf/paired/{{sample}}/{caller}.vcf.maf",project = project,genome_version = genome_version,caller = caller_list),
        #vcf5="{project}/{genome_version}/results/maf/{sample}/strelkasomatic.vcf.maf",
        #vcf6="{project}/{genome_version}/results/maf/{sample}/strelka.vcf.maf",
        #vcf7="{project}/{genome_version}/results/maf/{sample}/freebayes.vcf.maf",
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf"
        # filter_maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}_filter.maf"
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
        #vcf5="{project}/{genome_version}/results/maf/{sample}/strelkasomatic.vcf.maf",
        #vcf6="{project}/{genome_version}/results/maf/{sample}/strelka.vcf.maf",
        #vcf7="{project}/{genome_version}/results/maf/{sample}/freebayes.vcf.maf",
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}.maf",
        filter_maf="{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}_filter.maf"
    params:
        dir="{project}/{genome_version}/results/maf/unpaired/{sample}"
    conda:
        "clindet"
    script:
        "../scripts/merge_maf.R"

# rule merge_paired_ASCAT_cnv:
#     input:
#         vcf1=expand("{project}/{genome_version}/results/vcf/unpaired/{{sample}}/{caller}.vcf",project = project,genome_version = genome_version,caller = caller_list),
#         maf1=expand("{project}/{genome_version}/results/maf/unpaired/{{sample}}/{caller}.vcf.maf",project = project,genome_version = genome_version,caller = caller_list),
#         #vcf5="{project}/{genome_version}/results/maf/{sample}/strelkasomatic.vcf.maf",
#         #vcf6="{project}/{genome_version}/results/maf/{sample}/strelka.vcf.maf",
#         #vcf7="{project}/{genome_version}/results/maf/{sample}/freebayes.vcf.maf",
#         ref=config['resources']['hg19']['REFFA']
#     output:
#         maf="{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}.maf",
#         filter_maf="{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}_filter.maf"
#     params:
#         dir="{project}/{genome_version}/results/maf/unpaired/{sample}"
#     conda:
#         "clindet"
#     script:
#         "../scripts/merge_maf.R"
