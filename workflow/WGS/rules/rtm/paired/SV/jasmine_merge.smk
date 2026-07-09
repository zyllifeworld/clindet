## I consolidated the preprocessing steps for all structural variant callers into this smk file,
## then used Jasmine to generate a consensus set of SVs.

### for BRASS
### for fast combine only include translocation, if you want all types SVs, change the command to
### bcftools view -i 'SVTYPE="BND"' 
rule brass_pre_merge:
    input:
        log="{project}/{genome_version}/results/sv/paired/BRASS/{sample}/{sample}_brass.log"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/brass.vcf"
    conda: flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/sv/paired/BRASS/{sample}/{sample}_T_vs_{sample}_NC.annot.vcf.gz"
    shell:
        """
        bcftools view -i 'INFO/SVCLASS="translocation" || INFO/SVCLASS="inversion"' {params.vcf} -o {output.vcf} -O v
        """

## for delly
rule delly_pre_merge:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/delly.vcf"
    conda: flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        ref=config['resources'][genome_version]['REFFA'],
    shell:
        """
        bcftools annotate -x INFO/END,INFO/SVMETHOD,INFO/CONSENSUS -i 'SVTYPE="BND" || SVTYPE="INV"'  {input.vcf} -o {output.vcf} -O v
        """

## for svaba
rule svaba_pre_merge:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/svaba.vcf"
    conda: flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        ref=config['resources'][genome_version]['REFFA']
    shell:
        """
        bcftools view -i 'SVTYPE="BND" || SVTYPE="INV"' {input.vcf} -o {output.vcf} -O v
        """

## for gridss
rule gridss_pre_merge:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf.bgz"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/gridss.vcf"
    conda: flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        ref=config['resources'][genome_version]['REFFA']
    shell:
        """
        bcftools view -i 'SVTYPE="BND" || SVTYPE="INV"'  {input.vcf} -o {output.vcf} -O v
        """

### Manta
rule manta_pre_merge:
    input:
        tamp="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}-Manta.log",
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/Manta.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/somaticSV.vcf.gz"
    conda: flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        bcftools view -i 'SVTYPE="BND" & FILTER="PASS"' {input.vcf} -o {output.vcf} -O v
        """

### jasmine merge
rule jasmine_merge:
    input:
        vcf=expand("{project}/{genome_version}/results/sv/paired/merge/{{sample}}/{caller}.vcf",project = project,genome_version = genome_version,caller = somatic_sv_list),
    output:
        file_list="{project}/{genome_version}/results/sv/paired/merge/{sample}/merge_filelist.txt",
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/merge.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        wd="{project}/{genome_version}/results/sv/paired/merge/{sample}"
    singularity:
        flexible_container_img(config,['singularity','jasmine','sif'],image_url = config['singularity']['jasmine']['repo'])
    shell:
        """
        for f in $(echo "{input}")
        do
            echo "$f" >> "{output[0]}"
        done
        jasmine --normalize_type --allow_intrasample  file_list={output[0]} out_file={output.vcf}
        """
