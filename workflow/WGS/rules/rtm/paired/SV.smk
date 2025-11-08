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
        delly="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_delly.vcf",
        brass="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_brass.vcf",
        manta="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_manta.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_merge.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    shell:
        """
        /public/ClinicalExam/lj_sih/softwares/Jasmine/jasmine --output_genotypes --allow_intrasample  --comma_filelist  file_list={input.brass},{input.manta},{input.delly} out_file={output.vcf}
        """

# rule msi: 
#     input:
#         bed=config['bed_top100'],
#         list=config['list'],
#         tumor='bam/tumor.sort.mkdup.bam',
#     output:
#         'msi/sample.msi.result'
#     params:
#         '-l 1 -p 1 -q 1 -s 1 -b'
#     threads:
#         16
#     log:
#         'logs/msi.log'
#     singularity:
#         '/data_sas/ypu/git_repository/HEMECDx/test_data/tumorOnly/msisensor_v0.6.sif'
#     shell:
#         "msisensor msi {params} {threads} -e {input.bed} -d {input.list} "
#         "-t {input.tumor} -o {output} 1>{log} 2>&1 "



#rule SV_lumpy:
#    input:
#        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#    output:
        # sv="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf",
#    params:
#        ref=config['resources'][genome_version]['REFFA'],
#    shell:
#        """
#        {config[softwares][delly][call]} call -g {params.ref} {input.Tum} {input.NC} > {output}
#        
#        samtools fixmate -O bam - - | samtools sort -@ 30 -O bam -o {output.NC_spl}
#        samtools fixmate -O bam - - | samtools sort -@ 30 -O bam -o {output.NC_dis}
#        samtools fixmate -O bam - - | samtools sort -@ 30 -O bam -o {output.T_spl}
#        samtools fixmate -O bam - - | samtools sort -@ 30 -O bam -o {output.T_dis}
#       """


#rule SV_lumpy:
#    input:
#        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#    output:
#        sv="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf",
#    params:
#        ref=config['resources'][genome_version]['REFFA'],
#    shell:
#        """
#        {config[softwares][delly][call]} call -g {params.ref} {input.Tum} {input.NC} > {output}
#        """

