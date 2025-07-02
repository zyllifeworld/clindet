rule paired_sage:
    input:
        # Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        # indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        ref_genome=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}.sage.vcf.gz",
    params:
        output_dir="{project}/{genome_version}/results/cnv/paired/purple/{sample}/amber",
        loci=config['singularity']['hmftools'][genome_version]['amber']['loci'],
        high_confidence_bed=config['singularity']['hmftools'][genome_version]['sage']['high_confidence_bed'],
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['sage']['ensembl_data_dir'],
        coverage_bed=config['singularity']['hmftools'][genome_version]['sage']['coverage_bed'],
        hotspots=config['singularity']['hmftools'][genome_version]['sage']['hotspots'],
        panel_bed=config['singularity']['hmftools'][genome_version]['sage']['panel_bed'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['sage']['ref_genome_version'],
    threads: 30
    singularity:config['singularity']['hmftools']['sif']
    shell:
        """
        sage \
        -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -reference {wildcards.sample}_NC -reference_bam {input.NC} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -threads {threads} \
        -coverage_bed {params.coverage_bed} \
        -hotspots  {params.hotspots} \
        -panel_bed {params.panel_bed} \
        -high_confidence_bed {params.high_confidence_bed} \
        -output_vcf {output.vcf}
        """

rule sage_filter_pass:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}.sage.vcf.gz",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sage.vcf"
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    params:
        ref=config['resources'][genome_version]['REFFA']
    threads: 1
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf} 
        """


