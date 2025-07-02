#### svaba workflow
rule SV_svaba:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        dbsnp_indel=config['resources'][genome_version]['DBSNP_INDEL'],
    output:
        sv="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf",
        sv_uf="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.unfiltered.somatic.sv.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        wd="{project}/{genome_version}/results/sv/paired/svaba/{sample}"
    threads:20
    singularity: config['singularity']['svaba']['sif']
    shell:
        """
        svaba run -t {input.Tum} -n {input.NC} -p {threads} \
        -D {input.dbsnp_indel} -a {params.wd}/{wildcards.sample} -G {params.ref}
        """

rule anno_svaba:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf"
    output:
        anno="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}_SV_anno.bcf",
        query="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}_query.tsv.gz"
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


rule svanno_svaba:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}_SVAnnotate.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        db=config['softwares']['sansa'][genome_version]['db'],
        # gtf=config['resources'][genome_version]['GTF'],
        gtf='/public/ClinicalExam/lj_sih/projects/project_clindet/reference/b37/gencode.v19.annotation.gtf'
    shell:
        """
        {config[softwares][gatk4][call]} SVAnnotate \
        -V {input.vcf} \
        --protein-coding-gtf {params.gtf} \
        -O {output.vcf}
        """
