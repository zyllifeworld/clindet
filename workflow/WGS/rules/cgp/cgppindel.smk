pindel_normal_panel = config['singularity']['cgppindel'][genome_version]['WGS'].get("normal_panel", False)
if not pindel_normal_panel:
    ### pindel_normal_panel
    rule build_fake_bam:
        input:
            ref=config['resources'][genome_version]['REFFA'],
        output:
            bam="{project}/{genome_version}/results/recal/pindel_fake.bam"
        conda: config['conda']['clindet_main']# use samtools from clindet env
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
            simrep=config['singularity']['cgppindel'][genome_version]['simrep'],
            genes=config['singularity']['cgppindel'][genome_version]['genes']
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
            -seqtype WGS \
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
            gff3='analysis/normalPanel/pindel_{sample}'
        singularity:
            config['singularity']['cgppindel']['sif']
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
            log='{project}/{genome_version}/results/logs/paired/cgppindel_{sample}.log'
        threads: 20
        params:
            ref=config['resources'][genome_version]['REFFA'],
            simrep=config['singularity']['cgppindel'][genome_version]['simrep'],
            genes=config['singularity']['cgppindel'][genome_version]['genes'],
            filter=config['singularity']['cgppindel'][genome_version]['WGS']['filter'],
            softfil=config['singularity']['cgppindel'][genome_version]['softfil'],
            species=config['singularity']['cgppindel'][genome_version]['species']
        singularity:
            config['singularity']['cgppindel']['sif']
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
            -seqtype WGS \
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
            NP_gff3=config['singularity']['cgppindel'][genome_version]['WGS']['normal_panel']
        output:
            out_dir=directory('{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel'),
            log='{project}/{genome_version}/results/logs/paired/cgppindel_{sample}.log',
            vcf='{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
        threads: 20
        params:
            ref=config['resources'][genome_version]['REFFA'],
            simrep=config['singularity']['cgppindel'][genome_version]['simrep'],
            genes=config['singularity']['cgppindel'][genome_version]['genes'],
            filter=config['singularity']['cgppindel'][genome_version]['WGS']['filter'],
            softfil=config['singularity']['cgppindel'][genome_version]['softfil'],
            species=config['singularity']['cgppindel'][genome_version]['species']
        singularity:
            config['singularity']['cgppindel']['sif']
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
            -seqtype WGS \
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
        config['singularity']['cgppindel']['sif']
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
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel_filter.vcf"
    threads: 1
    conda: config['conda']['clindet_main']
    params:
        caller='cgppindel',
        vcf='{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
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