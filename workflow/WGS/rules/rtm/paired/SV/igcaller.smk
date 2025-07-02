#### igcaller workflow igcaller for B-cell
rule SV_igcaller:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        dbsnp_indel=config['resources'][genome_version]['DBSNP_INDEL']
    output:
        sv="{project}/{genome_version}/results/sv/paired/igcaller/{sample}/{sample}-T_IgCaller/{sample}-T_output_filtered.tsv"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        wd="{project}/{genome_version}/results/sv/paired/igcaller/{sample}"
    threads:20
    singularity: config['singularity']['igcaller']['sif']
    shell:
        """
            python3 $IGCALLER_DIR/IgCaller.py \
            --inputsFolder $IGCALLER_DIR/IgCaller_reference_files/ \
            --genomeVersion {wildcards.genome_version} \
            --chromosomeAnnotation ucsc \
            --bamN {input.NC} \
            --bamT {input.Tum} \
            --refGenome {params.ref} \
            --outputPath {params.wd} \
            -@ {threads}
        """