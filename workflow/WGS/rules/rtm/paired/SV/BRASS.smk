#### BRASS work flow
rule SV_brass_bamstat:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bas",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam.bas"
    singularity: config['singularity']['brass']['sif']
    shell:
        """
            bam_stats -i {input.Tum} -o {output.Tum}
            bam_stats -i {input.NC} -o {output.NC}
        """

rule SV_brass:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        Tumbas="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bas",
        NCbas="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam.bas",
        ascat="{project}/{genome_version}/results/cnv/paired/BRASS/{sample}/{sample}.ascat"
    output:
        log="{project}/{genome_version}/results/sv/paired/BRASS/{sample}/{sample}_brass.log"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        gc=config['singularity']['brass'][genome_version]['gc'],
        b=config['singularity']['brass'][genome_version]['b'],
        d=config['singularity']['brass'][genome_version]['d'],
        cb=config['singularity']['brass'][genome_version]['cb'],
        ct=config['singularity']['brass'][genome_version]['ct'],
        vi=config['singularity']['brass'][genome_version]['vi'],
        mi=config['singularity']['brass'][genome_version]['mi'],
        out_dir="{project}/{genome_version}/results/sv/paired/BRASS/{sample}"
    singularity: config['singularity']['brass']['sif']
    threads:20
    shell:
        """
        mkdir -p {params.out_dir}
        brass.pl -s human -as {wildcards.genome_version} -pr WGS \
        -b {params.b} \
        -c {threads} -o {params.out_dir} \
        -d {params.d} -g {params.ref} \
        -gc {params.gc} -cb {params.cb} \
        -vi {params.vi} -ct {params.ct} -mi {params.mi}\
        -t {input.Tum} \
        -n {input.NC} -ss {input.ascat}
        touch {output.log}
        """

rule brass_pre_merge:
    input:
        log="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/{sample}_brass.log"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_brass.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/sv/paired/BRASS/{sample}/{sample}_T_vs_{sample}_NC.annot.vcf.gz"
    shell:
        """
        zcat {params.vcf} | bcftools annotate  -i 'SVTYPE="BND"' > {output.vcf}
        """