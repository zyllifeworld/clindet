rule vcf_norm:
    input:
        vcf='{project}/{genome_version}/results/vcf/paired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf='{project}/{genome_version}/results/vcf_norm/paired/{sample}/{caller}.vcf',
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        gff3=lambda wildcards: get_config_value(
            config,
            ['resources', wildcards.genome_version, 'GFF'],
            default="",
            params='-g'
        ),
        filter_cmd=lambda wildcards: build_bcftools_filter_cmd(wildcards)
    shell:
        """
        if [ "{wildcards.caller}" = "vardict" ] || [ "{wildcards.caller}" = "VarDict" ]; then

            # Ensure fasta index exists
            if [ ! -s {input.ref}.fai ]; then
                samtools faidx {input.ref}
            fi

            contig_header="{output.vcf}.contigs.header"
            tmp_vcf="{output.vcf}.with_contig.tmp.vcf"

            # Generate contig header from reference fasta index
            awk 'BEGIN{{OFS=""}}{{print "##contig=<ID="$1",length="$2">"}}' \
                {input.ref}.fai > "$contig_header"

            # Only add contig header for VarDict VCF
            bcftools annotate \
                --header-lines "$contig_header" \
                {input.vcf} \
                -Ov \
                -o "$tmp_vcf"

            bcftools norm \
                -f {input.ref} \
                {params.gff3} \
                "$tmp_vcf" \
                -Ou | \
            bcftools view \
                {params.filter_cmd} \
                -Ov \
                -o {output.vcf}

            rm -f "$contig_header" "$tmp_vcf"

        else

            bcftools norm \
                -f {input.ref} \
                {params.gff3} \
                {input.vcf} \
                -Ou | \
            bcftools view \
                {params.filter_cmd} \
                -Ov \
                -o {output.vcf}

        fi
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
        flexible_conda_env(config,['conda','clindet_vep'],env_yaml = 'envs/clindet_vep.yaml')
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
        flexible_conda_env(config,['conda','clindet_vep'],env_yaml = 'envs/clindet_vep.yaml')
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
        maf1=expand("{project}/{genome_version}/results/maf/paired/{{sample}}/{caller}.vcf.maf",caller = caller_list,project = project,genome_version = genome_version),
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf"
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        dir="{project}/{genome_version}/results/maf/paired/{sample}"
    script: merge_maf_script[config["run_type"]]


rule paired_maf_report:
    input:
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
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        dir="{project}/{genome_version}/results/maf/unpaired/{sample}"
    script: merge_maf_script[config["run_type"]]


rule merge_paired_vcf:
    input:
        vcfs=expand("{project}/{genome_version}/results/vcf_norm/paired/{{sample}}/{caller}.vcf",caller = caller_list,project = project,genome_version = genome_version),
        ref=config['resources'][genome_version]['REFFA']
    output:
        merged_vcf="{project}/{genome_version}/results/vcf_norm/paired/{sample}/merge/{sample}.vcf"
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        dir="{project}/{genome_version}/results/maf/paired/{sample}"
    script: "../scripts/merge_caller_vcfs.py"


rule merge_unpaired_vcf:
    input:
        vcfs=expand("{project}/{genome_version}/results/vcf/unpaired/{{sample}}/{caller}.vcf",caller = tumor_only_caller,project = project,genome_version = genome_version),
        ref=config['resources'][genome_version]['REFFA']
    output:
        merged_vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/merge/{sample}.vcf"
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        dir="{project}/{genome_version}/results/maf/paired/{sample}"
    script: "../scripts/merge_caller_vcfs.py"