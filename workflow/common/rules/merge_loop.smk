rule vcf_norm:
    input:
        vcf='{project}/{genome_version}/results/vcf/paired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf='{project}/{genome_version}/results/vcf_norm/paired/{sample}/{caller}.vcf',
    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        # gff3="/public/ClinicalExam/lj_sih/projects/project_clindet/reference/b37/Homo_sapiens.GRCh37.87.gff3"
        ## if a gff3 for the genome is provided, 3'rule will be applyed
        gff3=lambda wildcards: get_config_value(
                config, 
                ['resources', wildcards.genome_version, 'GFF'],
                default="",params = '-g'
            ),
    shell:
        """
        bcftools norm -m -both -f {input.ref} {params.gff3} {input.vcf} -O v > {output.vcf}
        """

rule vcf2vcf:
    input:
        vcf='{project}/{genome_version}/results/vcf_norm/paired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf='{project}/{genome_version}/results/vcf2vcf/paired/{sample}/{caller}.vcf',
    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        name=vcf2vcf_name
    shell:
        """
        vcf2vcf.pl --input-vcf {input.vcf} \
        --output-vcf {output.vcf} --ref-fasta {input.ref} \
        {params.name} \
        --add-header '##INFO=<ID=CALLER,Number=1,Type=String,Description="Variant caller">' \
        --add-info 'CALLER={wildcards.caller}' \
        --new-tumor-id   {wildcards.sample}_T \
        --new-normal-id  {wildcards.sample}_NC \
        --retain-info SOMATIC,SS,I16,MQSB,CALLER
        """


merge_maf_script = {
    "wes": "../../WES/scripts/merge_maf.R",
    "wgs": "../../WGS/scripts/merge_maf.R",
    "rna": "../../RNA/scripts/merge_maf.R",
}

## now left align all variants
rule loop_vcf2maf_paired:
    input:
        vcf='{project}/{genome_version}/results/vcf_norm/paired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/paired/{sample}/{caller}.vcf.maf"
    conda:
        config['softwares']['vcf2maf']['conda']
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.{caller}_vcf2maf.benchmark.txt"
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


rule merge_paired_maf:
    input:
       # vcf1=expand("{project}/{genome_version}/results/vcf/paired/{{sample}}/{caller}.vcf",caller = caller_list,project = project,genome_version = genome_version),
        maf1=expand("{project}/{genome_version}/results/maf/paired/{{sample}}/{caller}.vcf.maf",caller = caller_list,project = project,genome_version = genome_version),
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf"
    conda: config['conda']['clindet_main']
    params:
        dir="{project}/{genome_version}/results/maf/paired/{sample}"
    script: merge_maf_script[config["run_type"]]

rule merge_paired_vcf:
    input:
        vcf=expand("{project}/{genome_version}/results/vcf2vcf/paired/{{sample}}/{caller}.vcf",caller = caller_list,project = project,genome_version = genome_version),
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf2vcf/paired/{sample}/merge/{sample}.vcf"
    conda: config['conda']['clindet_main']
    params:
        # dir="{project}/{genome_version}/results/vcf2vcf/paired/{sample}"
    script:
        '../scripts/merge_caller_vcf.py'

rule paired_maf_report:
    input:
       # vcf1=expand("{project}/{genome_version}/results/vcf/paired/{{sample}}/{caller}.vcf",caller = caller_list,project = project,genome_version = genome_version),
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf"
    output:
        report(
            "{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf.report",
            caption="../../../report/somedata.rst",
            category="merge MAF")
    shell:
        "wc -l {input.maf} > {output}"

# add bcftools normal

rule merge_unpaired_maf:
    input:
        vcf1=expand("{project}/{genome_version}/results/vcf/unpaired/{{sample}}/{caller}.vcf",caller = tumor_only_caller,project = project,genome_version = genome_version),
        maf1=expand("{project}/{genome_version}/results/maf/unpaired/{{sample}}/{caller}.vcf.maf",caller = tumor_only_caller,project = project,genome_version = genome_version),
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}.maf",
        # filter_maf="{project}/{genome_version}/results/maf/unpaired/{sample}/merge/{sample}_filter.maf"
    conda: config['conda']['clindet_main']
    params:
        dir="{project}/{genome_version}/results/maf/unpaired/{sample}"
    script: merge_maf_script[config["run_type"]]
