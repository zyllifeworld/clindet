## 
rule maf2vcf:
    input:
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf",
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.vcf"
    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        out_dir="{project}/{genome_version}/results/maf/paired/{sample}/merge"
    shell:
        """
        maf2vcf.pl --input-maf {input.maf} \
        --output-dir {params.out_dir} --ref-fasta {input.ref} \
        --output-vcf {output.vcf}  
        sed -i 's/ID=AD,Number=R/ID=AD,Number=2/g' {output.vcf} 
        """

rule pcgr_anno:
    input:
        vcf="{project}/{genome_version}/results/vcf2vcf/paired/{sample}/merge/{sample}.vcf",
        ref=config['resources'][genome_version]['REFFA']
    output:

    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        out_dir="{project}/{genome_version}/results/report/paired/{sample}/pcgr",
        refdata_dir=config['softwares_params'][genome_version]['pcgr']['refdata_dir'],
        vep_dir=
    shell:
        """
        pcgr \
            --refdata_dir /Users/you/dir2/data \
            --vep_dir /Users/you/dir3/vep/.vep \
            --output_dir /Users/you/dir4/pcgr_outputs \
            --sample_id {wildcards.sample} \
            --genome_assembly {} \
            --input_vcf {input.vcf} \
            --tumor_dp_tag TDP \
            --tumor_af_tag TVAF \
            --assay WES \
            --vcf2maf \
            --estimate_signatures \
            --estimate_msi \
            --estimate_tmb \
            --force_overwrite
        """

rule split_maf_snp_indel:
    input:
        maf="{project}/{genome_version}/results/maf/{sample}.maf"
    output:
        snp="{project}/{genome_version}/results/maf_split/{sample}.snp.maf",
        indel="{project}/{genome_version}/results/maf_split/{sample}.indel.maf"
conda: flexible_conda_env(config,['conda', 'clindet_main'],default="envs/clindet.yaml"),
    script:
        "../scripts/split_maf.R"
    
rule moalmanac_annotation:
    input:
        snp="{project}/{genome_version}/results/maf_split/{sample}.snp.maf",
        indel="{project}/{genome_version}/results/maf_split/{sample}.indel.maf"
    output:
    conda: flexible_conda_env(config,['conda', 'clindet_main'],default="envs/clindet.yaml"),
    script:
        """
        python /moalmanac/moalmanac.py \
        --patient_id "example" \
        --tumor_type "SKCM" \
        --stage "Metastatic" \
        --description "Example profile for interpretation with the Molecular Oncology Almanac" \
        --snv_handle "../example_data/example_patient.capture.somatic.snvs.maf" \
        --indel_handle "../example_data/example_patient.capture.somatic.indels.maf" \
        --bases_covered_handle "../example_data/example_patient.capture.somatic.coverage.txt" \
        --called_cn_handle "../example_data/example_patient.capture.somatic.called.cna.txt" \
        --fusion_handle "../example_data/example_patient.capture.somatic.seg.annotated" \
        --germline_handle "../example_data/example_patient.rna.star.fusions.txt" \
        --validation_handle "../example_data/example_patient.rna.somatic.snvs.maf" \
        --purity 0.85 \
        --ploidy 4.02 \
        --ms_status "msih" \
        --wgd \
        --config config.ini \
        --dbs annotation-databases.ini \
        --preclinical-dbs preclinical-databases.ini
    """"""