rule vcf_norm_germline:
    input:
        vcf='{project}/{genome_version}/results/vcf_germline/{sample_type}/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf='{project}/{genome_version}/results/vcf_germline_norm/{sample_type}/{sample}/{caller}.vcf',
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


rule loop_vcf2maf_germ_paired:
    input:
        vcf='{project}/{genome_version}/results/vcf_germline_norm/paired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf_germline/paired/{sample}/{caller}.vcf.maf"
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

rule loop_vcf2maf_germ_unpaired:
    input:
        vcf='{project}/{genome_version}/results/vcf_germline/unpaired/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/maf_germline/unpaired/{sample}/{caller}.vcf.maf"
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
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        dir="{project}/{genome_version}/results/maf_germline/unpaired/{sample}"
    script: merge_maf_script[config["run_type"]]
