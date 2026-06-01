#### gridss workflow
# param check section
def get_gridss_blacklist(wildcards):
    blacklist = config.get('softwares_params', {}).get(genome_version, {}).get('gridss', {}).get('blacklist', '')
    if not blacklist:
        r_b = ""
    else:
        r_b = " -b " + config['softwares_params'][genome_version]['gridss']['blacklist']
    return r_b

def get_gridss_pondir(wildcards):
    pondir = config.get('softwares_params', {}).get(genome_version, {}).get('gridss', {}).get('pondir', '')
    if not pondir:
        r_pon = ""
    else:
        r_pon = " --pondir " + config['softwares_params'][genome_version]['gridss']['pondir']
    return r_pon

rule SV_gridss:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/SV_gridss_{sample}.vcf",
    params:
        ref=config['resources'][genome_version]['REFFA'],
        wd="{project}/{genome_version}/results/sv/paired/gridss/{sample}",
        blacklist=get_gridss_blacklist
    singularity:
        flexible_container_img(config,['singularity','gridss','sif'],image_url = config['singularity']['gridss']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/sv/{sample}.gridss.benchmark.txt"
    shell:
        """
        gridss \
        -r {input.ref} \
        -j /usr/local/share/gridss-2.13.2-2/gridss.jar \
        -o {output.vcf} \
        -w {params.wd} \
        {params.blacklist} \
        {input.NC} \
        {input.Tum}
        """

# You can config pon based on https://github.com/PapenfussLab/gridss/issues/605
rule SV_gridss_filter:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/SV_gridss_{sample}.vcf"
    output:
        hvcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf.bgz",
        fullvcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_and_low_confidence_somatic.vcf.bgz",
    params:
        ref=config['resources'][genome_version]['REFFA'],
        wd="{project}/{genome_version}/results/sv/paired/gridss/{sample}",
        pondir=get_gridss_pondir,
        script=Path(str(workflow.current_basedir) + '/../../../../scripts/gridss/scripts'),
        # gridss will auto bgzip file,so will add bgz suffix,use params
        hvcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf",
        fullvcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_and_low_confidence_somatic.vcf"
    singularity:
        flexible_container_img(config,['singularity','gridss','sif'],image_url = config['singularity']['gridss']['repo'])
    shell:
        """
        /usr/local/share/gridss-2.13.2-2/gridss_somatic_filter \
        {params.pondir} \
        --input {input.vcf} \
        --output {params.hvcf} \
        --fulloutput {params.fullvcf} \
        -s {params.script} \
        -c {params.script} \
        -n 1 \
        -t 2
        """

#change vcf sample name for purple use
rule gridss_rename_tumor:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf.bgz",
    output:
        vcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic_rename.vcf.bgz",
    singularity:
        flexible_container_img(config,['singularity','gridss','sif'],image_url = config['singularity']['gridss']['repo'])
    params:
        vcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic_rename.vcf",
    script:
        "../../../../scripts/gridss/scripts/change_vcf_sample_name.R"