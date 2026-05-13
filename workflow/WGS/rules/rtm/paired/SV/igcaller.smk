def igcaller_extract_pp(wildcards, **kwargs):
    project = wildcards.project
    genome_version = wildcards.genome_version
    sample = wildcards.sample
    infile = f"{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_purity.ploidy.tsv"
    try:
        import pandas as pd
        pp_df = pd.read_csv(infile,sep="\t",engine="python")
        pp = list(pp_df.loc[0,['purity','ploidy']])
        return {
            "purity":pp[0],
            "ploidy":pp[1],
        }
    except ImportError:
        raise WorkflowError("Pandas is required to extract checksum from file.")

#### igcaller workflow igcaller for B-cell
if genome_version in ['hg19','b37','hg38','hg38_EBV']:
    rule SV_igcaller:
        input:
            Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
            NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
            dbsnp_indel=config['resources'][genome_version]['DBSNP_INDEL'],
            pp="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_purity.ploidy.tsv"
        output:
            sv="{project}/{genome_version}/results/sv/paired/igcaller/{sample}/{sample}-T_IgCaller/{sample}-T_output_filtered.tsv"
        params:
            ref=config['resources'][genome_version]['REFFA'],
            wd="{project}/{genome_version}/results/sv/paired/igcaller/{sample}",
            pp=igcaller_extract_pp,
            genome_version= config['softwares_params'][genome_version]['igcaller']['genome'],
            chromosomeAnnotation = config['softwares_params'][genome_version]['igcaller']['chromosomeAnnotation']
        threads:10
        singularity: config['singularity']['igcaller']['sif']
        benchmark:
            "{project}/{genome_version}/results/benchmarks/sv/{sample}.igcaller.benchmark.txt"
        shell:
            """
                python3 $IGCALLER_DIR/IgCaller.py \
                --inputsFolder $IGCALLER_DIR/IgCaller_reference_files/ \
                --genomeVersion {params.genome_version} \
                --chromosomeAnnotation {params.chromosomeAnnotation} \
                -p {params[pp][purity]} -kmb yes \
                --bamN {input.NC} \
                --bamT {input.Tum} \
                --refGenome {params.ref} \
                --outputPath {params.wd} \
                -@ {threads}
            """