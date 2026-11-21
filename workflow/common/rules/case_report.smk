##
def moalmanac_extract_pp(wildcards, **kwargs):
    """
    Resolve tumor purity/ploidy for MoAlmanac.

    Preference order:
      1. ASCAT  -> results/cnv/paired/ascat/{sample}/{sample}_purity.ploidy.tsv
      2. FACETS -> results/cnv/paired/facets/{sample}/{sample}.txt
         (purity model, i.e. the row with the largest cval)
    When neither yields a usable value, source is "NONE" and the shell block
    aborts with an explicit error.
    """
    import pandas as pd

    project = wildcards.project
    genome_version = wildcards.genome_version
    sample = wildcards.sample

    ascat_infile = f"{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_purity.ploidy.tsv"
    facets_infile = f"{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}.txt"

    def valid(purity, ploidy):
        try:
            return float(purity) > 0 and float(ploidy) > 0
        except (TypeError, ValueError):
            return False

    # 1) Prefer ASCAT when it found a solution
    if os.path.exists(ascat_infile):
        try:
            pp_df = pd.read_csv(ascat_infile, sep="\t", engine="python")
            purity, ploidy = pp_df.loc[0, "purity"], pp_df.loc[0, "ploidy"]
            if valid(purity, ploidy):
                return {"purity": float(purity), "ploidy": float(ploidy), "source": "ASCAT"}
        except (pd.errors.EmptyDataError, pd.errors.ParserError, KeyError, IndexError, OSError):
            pass  # empty/malformed ASCAT output -> fall back to FACETS

    # 2) Fall back to the FACETS purity model (largest cval, normally 1000)
    if os.path.exists(facets_infile):
        try:
            pp_df = pd.read_csv(facets_infile, sep="\t", engine="python")
            pp = pp_df.sort_values("cval", ascending=False).iloc[0]
            purity, ploidy = pp["purity"], pp["ploidy"]
            if valid(purity, ploidy):
                return {"purity": float(purity), "ploidy": float(ploidy), "source": "FACETS"}
        except (pd.errors.EmptyDataError, pd.errors.ParserError, KeyError, IndexError, OSError):
            pass

    return {"purity": "NA", "ploidy": "NA", "source": "NONE"}


# rule pcgr_anno:
#     input:
#         vcf="{project}/{genome_version}/results/vcf_norm/paired/{sample}/merge/{sample}.vcf",
#         ref=config['resources'][genome_version]['REFFA']
#     output:
#         "{project}/{genome_version}/results/vcf_norm/paired/{sample}/merge/{sample}.maf",
#     conda:
#         config['softwares']['vcf2maf']['conda']
#     params:
#         out_dir="{project}/{genome_version}/results/report/paired/{sample}/pcgr",
#         refdata_dir=config['softwares_params'][genome_version]['pcgr']['refdata_dir']
#     shell:
#         """
#         pcgr \
#             --refdata_dir /Users/you/dir2/data \
#             --vep_dir /Users/you/dir3/vep/.vep \
#             --output_dir /Users/you/dir4/pcgr_outputs \
#             --sample_id {wildcards.sample} \
#             --genome_assembly  \
#             --input_vcf {input.vcf} \
#             --tumor_dp_tag TDP \
#             --tumor_af_tag TVAF \
#             --assay WES \
#             --vcf2maf \
#             --estimate_signatures \
#             --estimate_msi \
#             --estimate_tmb \
#             --force_overwrite
#         """

rule split_maf_snp_indel:
    input:
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf"
    output:
        snp="{project}/{genome_version}/results/maf/paired/{sample}/split/{sample}.snp.maf",
        indel="{project}/{genome_version}/results/maf/paired/{sample}/split/{sample}.indel.maf"
    conda:
        flexible_conda_env(config,['conda', 'clindet_main'],env_yaml="envs/clindet.yaml")
    script:
        "../scripts/split_maf.R"


rule moalmanac_annotation:
    input:
        snp="{project}/{genome_version}/results/maf/paired/{sample}/split/{sample}.snp.maf",
        indel="{project}/{genome_version}/results/maf/paired/{sample}/split/{sample}.indel.maf",
        cna="{project}/{genome_version}/results/cnv/paired/facets/{sample}/{sample}.cna.txt",
        pp="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_purity.ploidy.tsv"
    output:
        html="{project}/{genome_version}/results/report/{sample}/moalmanac/{sample}.report.html",
        moalmanac_dir=directory("{project}/{genome_version}/results/report/{sample}/moalmanac/")
    singularity:
        flexible_container_img(config,['singularity','moalmanac','sif'],image_url = config['singularity']['moalmanac']['repo'])
    params:
        pp=moalmanac_extract_pp,
        outdir_abs=lambda wc, output: os.path.abspath(output.moalmanac_dir)
    shell:
        """
        if [ "{params[pp][source]}" = "NONE" ]; then
            echo "ERROR: {wildcards.sample} has no usable purity/ploidy: ASCAT result is missing/empty and FACETS result is missing/invalid" >&2
            exit 1
        fi
        echo "[moalmanac] purity/ploidy source: {params[pp][source]}" >&2
        python /moalmanac/moalmanac.py \
            --patient_id "{wildcards.sample}" \
            --tumor_type "SKCM" \
            --stage "Metastatic" \
            --description "Example profile for interpretation with the Molecular Oncology Almanac" \
            --snv_handle {input.snp} \
            --indel_handle {input.indel} \
            --called_cn_handle {input.cna} \
            --purity {params[pp][purity]} \
            --ploidy {params[pp][ploidy]} \
            --config /moalmanac/config.ini \
            --dbs /moalmanac/annotation-databases.ini \
            --output_directory {params.outdir_abs} \
            --preclinical-dbs /moalmanac/preclinical-databases.ini
        """