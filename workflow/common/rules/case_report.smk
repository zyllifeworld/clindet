## 
def moalmanac_extract_pp(wildcards, **kwargs):
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