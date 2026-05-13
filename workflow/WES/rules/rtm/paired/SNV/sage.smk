rule paired_sage:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref_genome=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.vcf.gz",
    params:
        high_confidence_bed=config['softwares_params'][genome_version]['hmftools']['sage']['high_confidence_bed'],
        ensembl_data_dir=config['softwares_params'][genome_version]['hmftools']['sage']['ensembl_data_dir'],
        coverage_bed=config['softwares_params'][genome_version]['hmftools']['sage']['coverage_bed'],
        hotspots=config['softwares_params'][genome_version]['hmftools']['sage']['hotspots'],
        panel_bed=config['softwares_params'][genome_version]['hmftools']['sage']['panel_bed'],
        ref_genome_version=config['softwares_params'][genome_version]['hmftools']['sage']['ref_genome_version'],
        xms=2000,
        xmx=2000
    threads: 30
    conda:
        flexible_conda_env(config,['conda','hmftools'],env_yaml = 'envs/hmftools.yaml')
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.sage.benchmark.txt"
    shell:
        """
        sage \
        -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -reference {wildcards.sample}_NC -reference_bam {input.NC} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -threads {threads} \
        -skip_msi_jitter \
        -coverage_bed {params.coverage_bed} \
        -hotspots  {params.hotspots} \
        -panel_bed {params.panel_bed} \
        -high_confidence_bed {params.high_confidence_bed} \
        -output_vcf {output.vcf} -write_bqr_plot
        """

rule pave_anno_sage:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.vcf.gz",
        ref_genome=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.pave.vcf.gz",
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    conda:
        flexible_conda_env(config,['conda','hmftools'],env_yaml = 'envs/hmftools.yaml')
    # singularity:"https://depot.galaxyproject.org/singularity/hmftools-pave:1.7.1--hdfd78af_0" # snakemake not work well for conda & singularity
    params:
        wd="{project}/{genome_version}/results/vcf/paired/{sample}/sage",
        driver_gene_panel=config['softwares_params'][genome_version]['hmftools']['purple']['driver_gene_panel'],
        high_confidence_bed=config['softwares_params'][genome_version]['hmftools']['sage']['high_confidence_bed'],
        ensembl_data_dir=config['softwares_params'][genome_version]['hmftools']['sage']['ensembl_data_dir'],
        coverage_bed=config['softwares_params'][genome_version]['hmftools']['sage']['coverage_bed'],
        hotspots=config['softwares_params'][genome_version]['hmftools']['sage']['hotspots'],
        panel_bed=config['softwares_params'][genome_version]['hmftools']['sage']['panel_bed'],
        ref_genome_version=config['softwares_params'][genome_version]['hmftools']['sage']['ref_genome_version'],
    threads: 8
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.pave.benchmark.txt"
    resources:
        mem_mb=lambda wildcards, input: max(100 * input.size_files_mb[0], 1000) # 100 times vcf file size
    shell:
        """
        pave -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        -sample {wildcards.sample} \
        -vcf_file {input.vcf} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -driver_gene_panel {params.driver_gene_panel} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -output_dir {params.wd} \
        -threads {threads}
        """


rule sage_filter_pass:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.vcf.gz",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    threads: 1
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf} 
        """

rule paired_pave_filter_pass:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage/{sample}.sage.pave.vcf.gz",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/pave.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    threads: 1
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf} 
        """

#### Sage germline calling
rule paired_sage_germline:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref_genome=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf_germline/paired/{sample}/sage/{sample}.sage.vcf.gz",
    params:
        high_confidence_bed=config['softwares_params'][genome_version]['hmftools']['sage']['high_confidence_bed'],
        ensembl_data_dir=config['softwares_params'][genome_version]['hmftools']['sage']['ensembl_data_dir'],
        coverage_bed=config['softwares_params'][genome_version]['hmftools']['sage']['coverage_bed'],
        hotspots=config['softwares_params'][genome_version]['hmftools']['sage']['hotspots'],
        panel_bed=config['softwares_params'][genome_version]['hmftools']['sage']['panel_bed'],
        ref_genome_version=config['softwares_params'][genome_version]['hmftools']['sage']['ref_genome_version'],
        xms=2000,
        xmx=2000
    threads: 30
    conda:
        flexible_conda_env(config,['conda','hmftools'],env_yaml = 'envs/hmftools.yaml')
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.sage.benchmark.txt"
    shell:
        """
        sage \
        -tumor {wildcards.sample} -tumor_bam {input.NC} \
        -reference {wildcards.sample}_T -reference_bam {input.Tum} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -threads {threads} \
        -skip_msi_jitter \
        -coverage_bed {params.coverage_bed} \
        -hotspots  {params.hotspots} \
        -panel_bed {params.panel_bed} \
        -high_confidence_bed {params.high_confidence_bed} \
        -output_vcf {output.vcf} -write_bqr_plot
        """

rule sage_germline_filter_pass:
    input:
        vcf="{project}/{genome_version}/results/vcf_germline/paired/{sample}/sage/{sample}.sage.vcf.gz",
    output:
        vcf="{project}/{genome_version}/results/vcf_germline/paired/{sample}/sage.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    threads: 1
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf} 
        """