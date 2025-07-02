#### gridss workflow
# param check section
def get_gridss_blacklist(wildcards):
    if config['singularity']['gridss'][genome_version]['blacklist'] == "":
        r_b = ""
    else:
        r_b = " -b " + config['singularity']['gridss'][genome_version]['blacklist']
    return r_b

def get_gridss_pondir(wildcards):
    if config['singularity']['gridss'][genome_version]['pondir'] == "":
        r_pon = ""
    else:
        r_pon = " --pondir " + config['singularity']['gridss'][genome_version]['pondir']
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
        # if else like  (exprs ? c1 : c2) in C++/javascript is OK too,eg on next line
        # blacklist=("" if config['singularity']['gridss'][genome_version]['blacklist'] == "" else " -b " + config['singularity']['gridss'][genome_version]['blacklist'])
    singularity: config['singularity']['gridss']['sif']
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
        hvcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf",
        # if else like  (exprs ? c1 : c2) in C++/javascript is OK too,eg on next line
        # blacklist=("" if config['singularity']['gridss'][genome_version]['blacklist'] == "" else " -b " + config['singularity']['gridss'][genome_version]['blacklist'])
    singularity: config['singularity']['gridss']['sif']
    shell:
        """
        /usr/local/share/gridss-2.13.2-2/gridss_somatic_filter \
        {params.pondir} \
        --input {input.vcf} \
        --output {params.hvcf} \
        --fulloutput {output.fullvcf} \
        -s /public/ClinicalExam/lj_sih/softwares/gridss/scripts/ \
        -c /public/ClinicalExam/lj_sih/softwares/gridss/scripts/ \
        -n 1 \
        -t 2
        """

#change vcf sample name for purple use
rule gridss_rename_tumor:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf.bgz",
    output:
        vcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic_rename.vcf.bgz",
    singularity: config['singularity']['gridss']['sif']
    params:
        vcf="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic_rename.vcf",
    script:
        "../../../../scripts/gridss/scripts/change_vcf_sample_name.R"