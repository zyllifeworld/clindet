rule loop_vcf2maf_rna:
    input:
        vcf='{project}/{genome_version}/results/mut/vcf/{sample}/{caller}.vcf',
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/mut/maf/{sample}/{caller}.vcf.maf"
    conda:
        flexible_conda_env(config,['conda','clindet_vep'],env_yaml = 'envs/clindet_vep.yaml')
    params:
        name=get_vcf_name,
        # vep_path=config['softwares_params'][genome_version]['vcf2maf']['vep'][genome_version]['vep_path'],
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

rule merge_rna_maf:
    input:
        vcf1=expand("{project}/{genome_version}/results/mut/vcf/{{sample}}/{caller}.vcf",caller = rna_caller_list,project = project,genome_version = genome_version),
        maf1=expand("{project}/{genome_version}/results/mut/maf/{{sample}}/{caller}.vcf.maf",caller = rna_caller_list,project = project,genome_version = genome_version),
        ref=config['resources'][genome_version]['REFFA']
    output:
        maf="{project}/{genome_version}/results/mut/maf/{sample}/merge/{sample}.maf",
    params:
        dir="{project}/{genome_version}/results/mut/maf/{sample}"
    script:
        "../scripts/merge_maf.R"
