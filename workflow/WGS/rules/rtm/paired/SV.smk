#### delly workflow
include: "SV/delly.smk"
#### BRASS work flow
ascat_config = config['softwares']['ascat_wgs'].get(genome_version, False)
if ascat_config and genome_version in ['b37','hg38']:
    include: "SV/BRASS.smk"

rule manta_pre_merge:
    input:
        tamp="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}-Manta.log",
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_manta.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/somaticSV.vcf.gz"
    shell:
        """
        zcat {params.vcf}  | bcftools annotate -i 'SVTYPE="BND" && FILTER="PASS" ' > {output.vcf}
        """


#### gridss workflow
include:"SV/gridss.smk"
#### svaba workflow
include:"SV/svaba.smk"

rule svanno_manta:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/somaticSV.vcf.gz"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/manta/{sample}/{sample}_SVAnnotate.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
    shell:
        """
        {config[softwares][gatk4][call]} SVAnnotate \
        -V {input.vcf} \
        --protein-coding-gtf {params.gtf} \
        -O {output.vcf}
        """


rule svanno_delly:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/{sample}_SVAnnotate.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
    shell:
        """
        {config[softwares][gatk4][call]} SVAnnotate \
        -V {input.vcf} \
        --protein-coding-gtf {params.gtf} \
        -O {output.vcf}
        """
#### igcaller workflow igcaller for B-cell
igcaller_config = config['singularity'].get('igcaller',{}).get('sif', False)
if igcaller_config:
    include:"SV/igcaller.smk"
### jasmine merge
rule jasmine_merge:
    input:
    brass="{project}/{genome_version}/results/sv/paired/BRASS/{sample}/{sample}_brass.log"  if 'BRASS' in somatic_sv_list else '',
    delly="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf"   if 'delly'  in somatic_sv_list else '',
    gridss="{project}/{genome_version}/results/sv/paired/gridss/{sample}/high_confidence_somatic.vcf.bgz" if 'gridss' in somatic_sv_list else '',
    svaba="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf"  if 'svaba' in somatic_sv_list else '', 
    manta="{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/somaticSV.vcf.gz"  if 'Manta' in somatic_sv_list else '',
    output:
        file_list="{project}/{genome_version}/results/sv/paired/merge/{sample}/merge_filelist.txt",
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_merge.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    shell:
        """
        set -xe
        # create a file list of VCF files
        for f in $(echo "{input}")
        do
            echo "$f" >> "{output[0]}"
        done
        jasmine --output_genotypes --allow_intrasample  file_list={output[0]} out_file={output.vcf}
        """
