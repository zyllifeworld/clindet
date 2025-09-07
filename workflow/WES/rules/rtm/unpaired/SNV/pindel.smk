# run with fake bam
rule build_unpaired_fake_bam:
    input:
        ref=config['resources'][genome_version]['REFFA'],
    output:
        bam="{project}/{genome_version}/results/vcf/unpaired/{sample}/pindel/{sample}-fake-NC.bam",
    conda:'clindet'# use samtools from clindet env
    shell:
        """
        (samtools dict {input.ref} | cut -f 1-4 && echo -e '@RG\tID:1\tSM:{wildcards.sample}_NC') | samtools view -S -bo {output.bam} -
        samtools index {output.bam}
        """

rule unpaired_PI_call: 
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/vcf/unpaired/{sample}/pindel/{sample}-fake-NC.bam",
        NP_gff3=config['singularity']['cgppindel'][genome_version]['WES']['normal_panel']
    output:
        out_dir=directory('{project}/{genome_version}/results/vcf/unpaired/{sample}/cgppindel'),
        log='{project}/{genome_version}/results/vcf/unpaired/{sample}/cgppindel_{sample}.log',
        vcf='{project}/{genome_version}/results/vcf/unpaired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
    threads: 20
    params:
        ref=config['resources'][genome_version]['REFFA'],
        simrep=config['singularity']['cgppindel'][genome_version]['simrep'],
        genes=config['singularity']['cgppindel'][genome_version]['genes'],
        filter=config['singularity']['cgppindel'][genome_version]['WES']['filter'],
        softfil=config['singularity']['cgppindel'][genome_version]['softfil'],
        species=config['singularity']['cgppindel'][genome_version]['species']
    singularity:
        config['singularity']['cgppindel']['sif']
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
        -seqtype WXS \
        -tumour {input.Tum} \
        -normal {input.NC} \
        -outdir {output.out_dir} \
        -cpus {threads} > {output.log}
        """

rule PI_ggz_unpaired: 
    input:
        log='{project}/{genome_version}/results/vcf/unpaired/{sample}/cgppindel_{sample}.log',
        germ_bed='{project}/{genome_version}/results/vcf/unpaired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.germline.bed'
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
## filter an format DP and AD tag
rule cgppindel_filter_somatic_unpaired:
    input:
        vcf='{project}/{genome_version}/results/vcf/unpaired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/cgppindel.vcf"
    threads: 1
    params:
        caller='cgppindel',
        vcf='{project}/{genome_version}/results/vcf/unpaired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
    shell:
        """
        bcftools view -e 'FILTER~"FF010"' {input.vcf} -Ov -o {output.vcf}
        """
        # bcftools view -e 'FILTER~"FF010" || FILTER~"FF009"' {input.vcf} -Ov -o {output.vcf}