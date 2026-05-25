try:
    pindel_normal_panel = config['softwares_params'][genome_version]['cgppindel'].get("normal_panel", False)
except KeyError:
    # 如果任何一级键不存在，返回默认值
    pindel_normal_panel = False

if not pindel_normal_panel:
    ### pindel_normal_panel
    rule build_fake_bam:
        input:
            ref=config['resources'][genome_version]['REFFA'],
        output:
            bam="{project}/{genome_version}/results/recal/pindel_fake.bam"
        conda: 
            flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
        shell:
            """
            (samtools dict {input.ref} | cut -f 1-4 && echo -e '@RG\tID:1\tSM:FAKE') | samtools view -S -bo {output.bam} -
            samtools index {output.bam}
            """

    rule PI_NP: 
        input:
            Tum="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
            NC="{project}/{genome_version}/results/recal/pindel_fake.bam",
        output:
            log="{project}/{genome_version}/analysis/pindel_normal/log/{sample}_pindel_NC.log",
            vcf="{project}/{genome_version}/analysis/pindel_normal/{sample}/FAKE_vs_{sample}_NC.vcf.gz"
        threads: 20
        params:
            ref=config['resources'][genome_version]['REFFA'],
            out_dir='{project}/{genome_version}/analysis/pindel_normal/{sample}',
            simrep= lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'simrep'],
                        default=""
                    ),
            genes=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'genes'],
                        default=""
                    ),
            seqtype = 'WXS' if config["run_type"] == 'wes' else "WGS"
        singularity:
            config['singularity']['cgppindel']['sif']
        shell:
            """
            pindel.pl -noflag \
            -reference {params.ref} \
            -simrep {params.simrep} \
            -genes {params.genes} \
            -assembly {genome_version} \
            -species Human \
            -seqtype {params.seqtype} \
            -tumour {input.NC} \
            -normal {input.Tum} \
            -outdir {params.out_dir} \
            -cpus {threads} > {output.log}
            """


    ### pindel_normal_panel_build
    rule PI_UM: 
        input:
            vcfs=expand("{project}/{genome_version}/analysis/pindel_normal/log/{sample}_pindel_NC.log",sample = paired_samples,project = project,genome_version = genome_version)
        output:
            '{project}/{genome_version}/analysis/normalPanel/pindel_{sample}.gff3.gz'# add a pseudo wildcards.sample for slurm log
        threads: 20
        params:
            ref=config['resources'][genome_version]['REFFA'],
            gff3='{project}/{genome_version}/analysis/normalPanel/pindel_{sample}'
        singularity:
            flexible_container_img(config,['singularity','cgppindel','sif'],image_url = config['singularity']['cgppindel']['repo'])
        shell:
            """
            pindel_np_from_vcf.pl -o {params.gff3} -samp_id NORMAL {wildcards.project}/{wildcards.genome_version}/analysis/pindel_normal/*/*.vcf.gz        
            """

    rule PI_call: 
        input:
            Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
            NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
            NP_gff3='{project}/{genome_version}/analysis/normalPanel/pindel_{sample}.gff3.gz'
        output:
            out_dir=directory('{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel'),
            log='{project}/{genome_version}/results/logs/paired/cgppindel_{sample}.log',
            vcf='{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
        threads: 20
        params:
            ref=config['resources'][genome_version]['REFFA'],
            simrep=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'simrep'],
                        default=""
                    ),
            genes=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'genes'],
                        default=""
                    ),
            filter=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'WES', 'filter'],
                        default=""
                    ),
            softfil=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'softfil'],
                        default=""
                    ),
            species=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'species'],
                        default=""
                    ),
            seqtype = 'WXS' if config["run_type"] == 'wes' else "WGS"
        singularity:
            flexible_container_img(config,['singularity','cgppindel','sif'],image_url = config['singularity']['cgppindel']['repo'])
        benchmark:
            "{project}/{genome_version}/results/benchmarks/snv/{sample}.cgppindel.benchmark.txt"
        shell:
            """
            pindel.pl \
            -reference {params.ref} \
            -simrep {params.simrep} \
            -genes {params.genes} \
            -exclude chrUn%,NC_007605,hs37d5,GL% \
            -unmatched {input.NP_gff3} \
            -filter {params.filter} \
            -softfil {params.softfil} \
            -assembly {wildcards.genome_version} \
            -species {params.species} \
            -seqtype {params.seqtype} \
            -tumour {input.Tum} \
            -normal {input.NC} \
            -outdir {output.out_dir} \
            -cpus {threads} > {output.log}
            """
else:
### pindel_normal_panel_pre-exist,just call
    rule PI_call: 
        input:
            Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
            NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
            NP_gff3=config['softwares_params'][genome_version]['cgppindel']['normal_panel']
        output:
            out_dir=directory('{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel'),
            log='{project}/{genome_version}/results/logs/paired/cgppindel_{sample}.log',
            vcf='{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
        threads: 20
        params:
            ref=config['resources'][genome_version]['REFFA'],
            simrep=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'simrep'],
                        default=""
                    ),
            genes=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'genes'],
                        default=""
                    ),
            filter=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'WES', 'filter'],
                        default=""
                    ),
            softfil=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'softfil'],
                        default=""
                    ),
            species=lambda wildcards: get_config_value(
                        config,
                        ['softwares_params',  wildcards.genome_version, 'cgppindel', 'species'],
                        default=""
                    ),
            seqtype = 'WXS' if config["run_type"] == 'wes' else "WGS"
        singularity:
            flexible_container_img(config,['singularity','cgppindel','sif'],image_url = config['singularity']['cgppindel']['repo'])
        benchmark:
            "{project}/{genome_version}/results/benchmarks/snv/{sample}.cgppindel.benchmark.txt"
        shell:
            """
            pindel.pl \
            -reference {params.ref} \
            -simrep {params.simrep} \
            -genes {params.genes} \
            -exclude chrUn% \
            -unmatched {input.NP_gff3} \
            -filter {params.filter} \
            -softfil {params.softfil} \
            -assembly {wildcards.genome_version} \
            -species {params.species} \
            -seqtype {params.seqtype} \
            -tumour {input.Tum} \
            -normal {input.NC} \
            -outdir {output.out_dir} \
            -cpus {threads} > {output.log}
            """

rule PI_ggz: 
    input:
        log='{project}/{genome_version}/results/logs/paired/cgppindel_{sample}.log',
        germ_bed='{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.germline.bed'
    output:
        log='{project}/{genome_version}/logs/paired/germline_bed_{sample}.log'
    threads: 20
    singularity:
        flexible_container_img(config,['singularity','cgppindel','sif'],image_url = config['singularity']['cgppindel']['repo'])
    shell:
        """
        bgzip -c {input.germ_bed} > {input.germ_bed}.gz
        tabix -p {input.germ_bed}.gz
        touch {output.log}
        """

rule cgppindel_filter_somatic:
    input:
        vcf='{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel.vcf"
    threads: 1
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    params:
        caller='cgppindel'
    shell:
        """
        bcftools view -e 'FILTER~"FF010"' {input.vcf} -Ov -o {output.vcf}
        """
        # loose filter rules, only filter FF010
        # 'FF010' 'Variant must not exist within the Unmatched Normal Panel'
        # 'FF001' 'Pass if Mt > Wt Reads: Likely GERMLINE',
        # 'FF003'  'Tum low call count strand bias check'
        # 'FF006' 'Small call excessive repeat check: Fail if Length <= 4 and Repeats > 9'
        # 'FF008' Wildtype contamination: Fail when wt reads > 5% mt reads.
        # 'FF017' 'Variant must not overlap with a simple repeat'