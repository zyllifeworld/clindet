#### delly workflow
rule SV_delly:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        sv="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf",
    params:
        ref=config['resources'][genome_version]['REFFA'],
    benchmark:
        "{project}/{genome_version}/results/benchmarks/sv/{sample}.delly.benchmark.txt"
    singularity: config['singularity']['delly']['sif']
    shell:
        """
        delly call -g {params.ref}  {input.Tum} {input.NC} > {output}
        """

rule SV_delly_germ:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        sv="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_germ.vcf",
    params:
        ref=config['resources'][genome_version]['REFFA']
    benchmark:
        "{project}/{genome_version}/results/benchmarks/sv/{sample}.dellygerm.benchmark.txt"
    singularity: config['singularity']['delly']['sif']
    shell:
        """
        {config[softwares][delly][call]} call -g {params.ref} {input.NC} -q 10 -s 15 -n > {output}
        """

rule delly_filter:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    conda: config['conda']['clindet_main']
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf}
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
