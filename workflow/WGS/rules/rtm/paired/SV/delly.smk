#### delly workflow
rule SV_delly:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        bcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.bcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
    resources:
        shell_exec="sh"
    benchmark:
        "{project}/{genome_version}/results/benchmarks/sv/{sample}.delly.benchmark.txt"
    singularity:
        flexible_container_img(config,['singularity','delly','sif'],image_url = config['singularity']['delly']['repo'])
    shell:
        """
        delly sr \
            -g {params.ref} \
            -q 20 \
            -s 15 \
            -o {output.bcf} \
            {input.Tum} {input.NC}
        """

rule SV_delly_sample_tsv:
    input:
        bcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.bcf"
    output:
        tsv="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.samples.tsv"
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        bcftools query -l {input.bcf} > {output.tsv}.samples
        awk 'NR==1{{print $1"\ttumor"}}
             NR==2{{print $1"\tcontrol"}}' \
             {output.tsv}.samples > {output.tsv}

        rm -f {output.tsv}.samples
        """

rule SV_delly_filter_somatic:
    input:
        bcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.bcf",
        samples="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.samples.tsv",
    output:
        bcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_somatic.bcf"
    params:
        min_af=0.05
    benchmark:
        "{project}/{genome_version}/results/benchmarks/sv/{sample}.delly.somatic_filter.benchmark.txt"
    singularity:
        flexible_container_img(config,['singularity','delly','sif'],image_url = config['singularity']['delly']['repo'])
    shell:
        """
        delly filter \
            -p \
            -f somatic \
            -a {params.min_af} \
            -s {input.samples} \
            -o {output.bcf} \
            {input.bcf}
        """

rule SV_delly_to_vcf:
    input:
        bcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_somatic.bcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_somatic.vcf",
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        bcftools index -f {input.bcf}
        bcftools view \
            -f PASS \
            -O v \
            -o {output.vcf} \
            {input.bcf}
        """

rule SV_delly_germ:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        sv="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_germ.vcf",
    params:
        ref=config['resources'][genome_version]['REFFA']
    resources:
        shell_exec="sh"
    benchmark:
        "{project}/{genome_version}/results/benchmarks/sv/{sample}.dellygerm.benchmark.txt"
    singularity:
        flexible_container_img(config,['singularity','delly','sif'],image_url = config['singularity']['delly']['repo'])
    shell:
        """
        delly call -g {params.ref} {input.NC} -q 10 -s 15 -n > {output}
        """

rule delly_filter:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_somatic.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf}
        """

rule delly2bnd:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter_bnd.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        script=workflow.source_path("../../../../scripts/delly2bnd.py")
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        python {params.script} -v {input.vcf} -r {params.ref} -o {output.vcf}
        """

sansa_config = config['softwares'].get('sansa',{}).get(genome_version, False)
if sansa_config:
    rule SV_sansa_annodelly:
        input:
            vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter.vcf"
        output:
            anno="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_anno_{sample}.bcf",
            query="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/query_{sample}.tsv.gz"
        params:
            ref=config['resources'][genome_version]['REFFA'],
            db=config['softwares']['sansa'][genome_version]['db'],
            g=config['softwares']['sansa'][genome_version]['g'],
            t=10000
        shell:
            """
            {config[softwares][sansa][call]} annotate -i Name  -g {params.g} -t {params.t} \
            -a {output.anno} -o {output.query} {input.vcf} 
            """
