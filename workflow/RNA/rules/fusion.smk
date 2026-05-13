rule arriba_fusion:
    input:
        bam="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam",
        stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log"
    output:
        tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion.tsv"
        # dis_tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion_discarded.tsv"
    # conda:
    #     config['softwares']['samtools']['conda']
    params:
        temp_directory=config['params']['java']['temp_directory'],
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
        blacklist=config['softwares_params'][genome_version]['arriba']['database']['blacklist'],
        known_fus=config['softwares_params'][genome_version]['arriba']['database']['known_fusions'],
        pro_dom=config['softwares_params'][genome_version]['arriba']['database']['protein_domains']
    singularity:
        config['singularity']['arriba']['sif']
    shell:
        """
            /arriba_v2.4.0/arriba \
            -x {input.bam}  -o {output.tsv} \
            -a {params.ref} -g {params.gtf} \
            -b {params.blacklist} -k {params.known_fus} \
            -t {params.known_fus} \
            -p {params.pro_dom}    
        """

rule arriba_draw:
    input:
        # bam="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam",
        bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.sort.bam",
        stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log",
        tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion.tsv"
    output:
        pdf="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion.pdf"
    params:
        temp_directory=config['params']['java']['temp_directory'],
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
        c_band=config['softwares_params'][genome_version]['arriba']['database']['cytobands'],
        blacklist=config['softwares_params'][genome_version]['arriba']['database']['blacklist'],
        known_fus=config['softwares_params'][genome_version]['arriba']['database']['known_fusions'],
        pro_dom=config['softwares_params'][genome_version]['arriba']['database']['protein_domains']
    singularity:
        config['singularity']['arriba']['sif']
    shell:
        """   
            /arriba_v2.4.0/draw_fusions.R \
            --fusions={input.tsv}\
            --alignments={input.bam} \
            --output={output.pdf} \
            --annotation={params.gtf} \
            --cytobands={params.c_band} \
            --proteinDomains={params.pro_dom}
        """

rule TRUST4_TBCR:
    input:
        bam="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam",
        stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log"
    output:
        "{project}/{genome_version}/results/IG/TRUST4/{sample}_report.tsv"
    conda:
        flexible_conda_env(config,['conda','rna'],env_yaml = 'envs/rsem.yaml')
    threads:10
    params:
        temp_directory=config['params']['java']['temp_directory'],
        ref=config['softwares_params'][genome_version]['trust4']['ref'],
        f=config['softwares_params'][genome_version]['trust4']['f'],
        oref="{project}/{genome_version}/results/IG/TRUST4/{sample}"
    shell:
        """
        {config[softwares][trust4][call]} \
            -b {input.bam}  -f {params.f} \
            -o {params.oref} --ref {params.ref}  -t {threads}
        """
