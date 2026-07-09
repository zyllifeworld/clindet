"""Setup rules for Clindet bootstrap assets.

This file is meant to be included or run as a standalone setup workflow.
It converts the original bootstrap shell script into restartable Snakemake
rules while keeping the same resource layout under `resources/` and `build_log/`.
"""

shell.executable("/bin/bash")
SOFT_DIR = "resources/softwares"
BUILD_LOG = "build_log/b37"
REF_B37 = "resources/ref_genome/b37"

rule build_b37_ref:
    input:
        # expand(f"{BUILD_LOG}/{{env_name}}.env.done", env_name=sorted(ENV_SPECS)),
        f"{BUILD_LOG}/download_b37_hmftools.log",
        f"{BUILD_LOG}/download_b37_gatk.log",
        f"{BUILD_LOG}/download_gatk.log",
        f"{BUILD_LOG}/download_ascat.log",
        f"{BUILD_LOG}/download_sanger.log",
        f"{BUILD_LOG}/download_vep.log",
        f"{BUILD_LOG}/mass_config.log",
        f"{BUILD_LOG}/rna_edit_vcf.log",
        f"{BUILD_LOG}/rsem_star_index.log",
        f"{BUILD_LOG}/star_index.log",
        f"{BUILD_LOG}/kallisto_salmon_index.log",
        f"{BUILD_LOG}/download_bed_files.log",
        f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta",
        f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict",
        f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai",
        f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta.amb",
        f"{REF_B37}/Homo_sapiens_assembly19.dbsnp138.vcf",
        f"{REF_B37}/Homo_sapiens.GRCh37.87.gtf",
        f"{REF_B37}/Homo_sapiens.GRCh37.87.gff3",
        f"{REF_B37}/Homo_sapiens.GRCh37.cdna.all.fa",
        f"{REF_B37}/Homo_sapiens.GRCh37.cdna.all.fa.gz"


rule download_b37_reference:
    output:
        f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta",
        f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict",
        f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai",
    log:
        f"{BUILD_LOG}/download_b37_reference.log"
    conda:
        flexible_conda_env(config,['conda','gsutils'],env_yaml = 'envs/gsutils.yaml')
    shell:
        r"""
        mkdir -p {REF_B37}
        gsutil -m cp -r -n \
          "gs://hmf-public/HMFtools-Resources/ref_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta" \
          "gs://hmf-public/HMFtools-Resources/ref_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict" \
          "gs://hmf-public/HMFtools-Resources/ref_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai" \
          {REF_B37} &> {log}
        """


rule download_b37_hmftools:
    output:
        touch(f"{BUILD_LOG}/download_b37_hmftools.log"),
    log:
        f"{BUILD_LOG}/download_b37_hmftools.run.log"
    conda:
        flexible_conda_env(config,['conda','gsutils'],env_yaml = 'envs/gsutils.yaml')
    shell:
        r"""
        mkdir -p {BUILD_LOG}
        mkdir -p {REF_B37}/hmf_pipeline_resources
        {{
            gsutil -m cp \
              "gs://hmf-public/HMFtools-Resources/pipeline/oncoanalyser/2.0/37/hmf_panel_resources.tso500.37_v2.0.0--3.tar.gz" \
              "gs://hmf-public/HMFtools-Resources/pipeline/oncoanalyser/2.0/37/hmf_pipeline_resources.37_v2.0.0--3.tar.gz" \
              {REF_B37}
            wget -P {REF_B37}/hmf_pipeline_resources -c https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/GermlineHetPon.37.vcf.gz
            wget -P {REF_B37}/hmf_pipeline_resources -c https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/Amber.snpcheck.37.vcf
            tar -xzvf {REF_B37}/hmf_panel_resources.tso500.37_v2.0.0--3.tar.gz --strip-components 1 -C {REF_B37}/hmf_pipeline_resources/
            tar -xzvf {REF_B37}/hmf_pipeline_resources.37_v2.0.0--3.tar.gz --strip-components 1 -C {REF_B37}/hmf_pipeline_resources/
        }} &> {log}
        touch {output}
        """


rule download_b37_gatk:
    output:
        dbsnp_vcf=f"{REF_B37}/Homo_sapiens_assembly19.dbsnp138.vcf",
        dbsnp_idx=f"{REF_B37}/Homo_sapiens_assembly19.dbsnp138.vcf.idx",
        indels_gz=f"{REF_B37}/1000G_phase1.indels.b37.vcf.gz",
        indels_tbi=f"{REF_B37}/1000G_phase1.indels.b37.vcf.gz.tbi",
        marker=touch(f"{BUILD_LOG}/download_b37_gatk.log"),
    log:
        f"{BUILD_LOG}/download_b37_gatk.run.log"
    conda:
        flexible_conda_env(config,['conda','gsutils'],env_yaml = 'envs/gsutils.yaml')
    shell:
        r"""
        mkdir -p {BUILD_LOG}
        {{
            gsutil -m cp -r -n \
              "gs://gcp-public-data--broad-references/hg19/v0/1000G_omni2.5.b37.vcf.gz" \
              "gs://gcp-public-data--broad-references/hg19/v0/1000G_omni2.5.b37.vcf.gz.tbi" \
              "gs://gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz" \
              "gs://gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz.tbi" \
              "gs://gcp-public-data--broad-references/hg19/v0/1000G_reference_panel" \
              "gs://gcp-public-data--broad-references/hg19/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz" \
              "gs://gcp-public-data--broad-references/hg19/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz.tbi" \
              "gs://gcp-public-data--broad-references/hg19/v0/ExomeContam.vcf" \
              "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dbsnp138.vcf" \
              "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dbsnp138.vcf.idx" \
              {REF_B37}

            gsutil -m cp -n -r \
                "gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf" \
                "gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf.idx" \
                "gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf" \
                "gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx" \
                "gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf" \
                "gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx" \
              {REF_B37}

            gsutil -m cp -n \
                "gs://gatk-legacy-bundles/b37/1000G_omni2.5.b37.vcf" \
                "gs://gatk-legacy-bundles/b37/1000G_omni2.5.b37.vcf.gz" \
                "gs://gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf.gz" \
                "gs://gatk-legacy-bundles/b37/1000G_phase3_v4_20130502.sites.vcf.gz" \
                "gs://gatk-legacy-bundles/b37/1000G_phase3_v4_20130502.sites.vcf.gz.tbi" \
                "gs://gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz" \
                "gs://gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz.md5" \
                "gs://gatk-legacy-bundles/b37/dbsnp_138.b37.vcf" \
                "gs://gatk-legacy-bundles/b37/hapmap_3.3.b37.vcf" \
              {REF_B37}

            gunzip -c {output.indels_gz} > {REF_B37}/1000G_phase1.indels.b37.vcf
            bgzip -k -f -o {output.indels_gz} {REF_B37}/1000G_phase1.indels.b37.vcf
            tabix -f {output.indels_gz}
        }} &> {log}
        """


rule download_gatk_software:
    output:
        directory(f"{SOFT_DIR}/gatk"),
        marker=touch(f"{BUILD_LOG}/download_gatk.log"),
    log:
        f"{BUILD_LOG}/download_gatk.run.log"
    shell:
        r"""
        mkdir -p {SOFT_DIR}
        if [ ! -x {SOFT_DIR}/gatk/gatk ]; then
            GATK_version="4.6.2.0"
            wget -P {SOFT_DIR} -c https://github.com/broadinstitute/gatk/releases/download/${{GATK_version}}/gatk-${{GATK_version}}.zip
            unzip -o {SOFT_DIR}/gatk-${{GATK_version}}.zip -d {SOFT_DIR}
            mv -f {SOFT_DIR}/gatk-${{GATK_version}} {SOFT_DIR}/gatk
        fi
        echo "GATK ready" > {log}
        """


rule download_ascat_refdata:
    output:
        touch(f"{BUILD_LOG}/download_ascat.log")
    log:
        f"{BUILD_LOG}/download_ascat.run.log"
    shell:
        r"""
        mkdir -p {REF_B37}/ASCAT/WES
        mkdir -p {REF_B37}/ASCAT/WGS
        {{
            wget -P {REF_B37}/ASCAT/WES -c https://zenodo.org/records/14008443/files/G1000_alleles_WES_hg19.zip
            wget -P {REF_B37}/ASCAT/WES -c https://zenodo.org/records/14008443/files/G1000_loci_WES_hg19.zip
            wget -P {REF_B37}/ASCAT/WES -c https://zenodo.org/records/14008443/files/GC_G1000_WES_hg19.zip
            wget -P {REF_B37}/ASCAT/WES -c https://zenodo.org/records/14008443/files/RT_G1000_WES_hg19.zip
            unzip -o -d {REF_B37}/ASCAT/WES {REF_B37}/ASCAT/WES/G1000_alleles_WES_hg19.zip
            unzip -o -d {REF_B37}/ASCAT/WES {REF_B37}/ASCAT/WES/G1000_loci_WES_hg19.zip
            unzip -o -d {REF_B37}/ASCAT/WES {REF_B37}/ASCAT/WES/GC_G1000_WES_hg19.zip
            unzip -o -d {REF_B37}/ASCAT/WES {REF_B37}/ASCAT/WES/RT_G1000_WES_hg19.zip

            wget -P {REF_B37}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/G1000_alleles_WGS_hg19.zip
            wget -P {REF_B37}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/G1000_loci_WGS_hg19.zip
            wget -P {REF_B37}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/GC_G1000_WGS_hg19.zip
            wget -P {REF_B37}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/RT_G1000_WGS_hg19.zip
            unzip -o -d {REF_B37}/ASCAT/WGS {REF_B37}/ASCAT/WGS/G1000_alleles_WGS_hg19.zip
            unzip -o -d {REF_B37}/ASCAT/WGS {REF_B37}/ASCAT/WGS/G1000_loci_WGS_hg19.zip
            unzip -o -d {REF_B37}/ASCAT/WGS {REF_B37}/ASCAT/WGS/GC_G1000_WGS_hg19.zip
            unzip -o -d {REF_B37}/ASCAT/WGS {REF_B37}/ASCAT/WGS/RT_G1000_WGS_hg19.zip
        }} &> {log}
        touch {output}
        """


rule download_sanger_refdata:
    output:
        touch(f"{BUILD_LOG}/download_sanger.log")
    log:
        f"{BUILD_LOG}/download_sanger.run.log"
    shell:
        r"""
        mkdir -p {REF_B37}/Sanger
        {{
            wget -P {REF_B37}/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/SNV_INDEL_ref_GRCh37d5-fragment.tar.gz
            wget -P {REF_B37}/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/CNV_SV_ref_GRCh37d5_brass6+.tar.gz
            wget -P {REF_B37}/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/core_ref_GRCh37d5.tar.gz
            wget -P {REF_B37}/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz
            wget -P {REF_B37}/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/support-files/cgpPindel/cgpPindel_CPBI_refarea.tar.gz
            wget -P {REF_B37}/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/support-files/CPIB/caveman/cgpCaVEManWrapper_CPBI_refarea.tar.gz

            tar -zxvf {REF_B37}/Sanger/SNV_INDEL_ref_GRCh37d5-fragment.tar.gz -C {REF_B37}/Sanger
            tar -zxvf {REF_B37}/Sanger/CNV_SV_ref_GRCh37d5_brass6+.tar.gz -C {REF_B37}/Sanger
            tar -zxvf {REF_B37}/Sanger/core_ref_GRCh37d5.tar.gz -C {REF_B37}/Sanger
            tar -zxvf {REF_B37}/Sanger/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz -C {REF_B37}/Sanger
            tar -zxvf {REF_B37}/Sanger/cgpPindel_CPBI_refarea.tar.gz -C {REF_B37}/Sanger
            tar -zxvf {REF_B37}/Sanger/cgpCaVEManWrapper_CPBI_refarea.tar.gz -C {REF_B37}/Sanger

            mkdir -p {REF_B37}/Sanger/SNV_INDEL_ref/caveman/flagging
            cp {REF_B37}/Sanger/VAGrENT_ref_GRCh37d5_ensembl_75/vagrent/gene_regions.bed* {REF_B37}/Sanger/SNV_INDEL_ref/caveman/flagging/
            cp {REF_B37}/Sanger/VAGrENT_ref_GRCh37d5_ensembl_75/vagrent/codingexon_regions.sub.bed* {REF_B37}/Sanger/SNV_INDEL_ref/caveman/flagging/
        }} &> {log}
        touch {output}
        """


rule download_vep_cache:
    output:
        touch(f"{BUILD_LOG}/download_vep.log")
    log:
        f"{BUILD_LOG}/download_vep.run.log"
    params:
        cache_num=str(config["softwares_params"]['b37']["vcf2maf"]["vep"]["cache_version"]).replace("v", ""),
        cache_dir=f"{REF_B37}/vep/v{str(config['softwares_params']['b37']['vcf2maf']['vep']['cache_version']).replace('v', '')}",
        vep_dir=f"{REF_B37}/vep"
    shell:
        """
        (
            set -euo pipefail

            mkdir -p {params.cache_dir}

            wget \
                -P {params.cache_dir} \
                -c https://ftp.ensembl.org/pub/release-{params.cache_num}/variation/indexed_vep_cache/homo_sapiens_vep_{params.cache_num}_GRCh37.tar.gz

            tar -xzvf \
                {params.cache_dir}/homo_sapiens_vep_{params.cache_num}_GRCh37.tar.gz \
                -C {params.vep_dir}
        ) &> {log}

        touch {output}
        """


rule download_grch37_gtf:
    output:
        f"{REF_B37}/Homo_sapiens.GRCh37.87.gtf",
    log:
        f"{BUILD_LOG}/download_grch37_gtf.run.log"
    shell:
        r"""
        mkdir -p {REF_B37}
        wget -P {REF_B37} -c https://ftp.ensembl.org/pub/grch37/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
        gunzip -c {REF_B37}/Homo_sapiens.GRCh37.87.gtf.gz > {output}
        echo "GTF ready" > {log}
        """

rule download_grch37_gff3:
    output:
        f"{REF_B37}/Homo_sapiens.GRCh37.87.gff3",
    log:
        f"{BUILD_LOG}/download_grch37_gff.run.log"
    shell:
        r"""
        mkdir -p {REF_B37}
        wget -P {REF_B37} -c https://ftp.ensembl.org/pub/grch37/release-114/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
        gunzip -c {REF_B37}/Homo_sapiens.GRCh37.87.gff3.gz > {output}
        echo "GFF3 ready" > {log}
        """

rule download_cdna_fasta:
    output:
        fa=f"{REF_B37}/Homo_sapiens.GRCh37.cdna.all.fa",
        gz=f"{REF_B37}/Homo_sapiens.GRCh37.cdna.all.fa.gz",
    log:
        f"{BUILD_LOG}/download_cdna_fasta.run.log"
    shell:
        r"""
        mkdir -p {REF_B37}
        wget -P {REF_B37} -c https://ftp.ensembl.org/pub/grch37/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz
        gunzip -c {REF_B37}/Homo_sapiens.GRCh37.cdna.all.fa.gz > {output.fa}
        echo "cDNA ready" > {log}
        """


rule download_rna_edit_vcf:
    output:
        vcf_gz=f"{REF_B37}/b37.RNAediting.vcf.gz",
        tbi=f"{REF_B37}/b37.RNAediting.vcf.gz.tbi",
        rename_tsv=f"{REF_B37}/hg19_to_g1kv37.tsv",
        marker=touch(f"{BUILD_LOG}/rna_edit_vcf.log"),
    log:
        f"{BUILD_LOG}/rna_edit_vcf.run.log"
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        r"""
        mkdir -p {REF_B37}
        wget -P {REF_B37} -c https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/rna_editing/GRCh37.RNAediting.vcf.gz
        wget -P {REF_B37} -c https://raw.githubusercontent.com/lindenb/jvarkit/refs/heads/master/src/main/resources/chromnames/hg19_to_g1kv37.tsv
        gunzip -c {REF_B37}/GRCh37.RNAediting.vcf.gz > {REF_B37}/GRCh37.RNAediting.vcf
        bgzip -f {REF_B37}/GRCh37.RNAediting.vcf
        tabix -f {REF_B37}/GRCh37.RNAediting.vcf.gz
        bcftools annotate --rename-chrs {REF_B37}/hg19_to_g1kv37.tsv {REF_B37}/GRCh37.RNAediting.vcf.gz -Ov -o {REF_B37}/b37.RNAediting.vcf
        bgzip -f {REF_B37}/b37.RNAediting.vcf
        tabix -f {REF_B37}/b37.RNAediting.vcf.gz
        echo "RNA editing VCF ready" > {log}
        """


rule build_bwa_index:
    input:
        fasta=f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta",
    output:
        f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta.amb",
    log:
        f"{BUILD_LOG}/bwa_index.run.log"
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        r"""
        bwa index {input.fasta}
        echo "BWA index ready" > {log}
        """


rule build_star_index:
    input:
        fasta=f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta",
        gtf=f"{REF_B37}/Homo_sapiens.GRCh37.87.gtf",
    output:
        directory(f"{REF_B37}/STAR/b37"),
        marker=touch(f"{BUILD_LOG}/star_index.log"),
    log:
        f"{BUILD_LOG}/star_index.run.log"
    conda:
        flexible_conda_env(config,['conda','rna'],env_yaml = 'envs/rsem.yaml')
    shell:
        r"""
        mkdir -p {REF_B37}/STAR/b37
        STAR \
          --runThreadN 20 \
          --runMode genomeGenerate \
          --genomeFastaFiles {input.fasta} \
          --sjdbOverhang 100 --genomeSAindexNbases 2 \
          --sjdbGTFfile {input.gtf} \
          --genomeDir {REF_B37}/STAR/b37
        echo "STAR index ready" > {log}
        """


rule build_rsem_reference:
    input:
        fasta=f"{REF_B37}/Homo_sapiens.GRCh37.GATK.illumina.fasta",
        gtf=f"{REF_B37}/Homo_sapiens.GRCh37.87.gtf",
    output:
        marker=touch(f"{BUILD_LOG}/rsem_star_index.log"),
    log:
        f"{BUILD_LOG}/rsem_star_index.run.log"
    conda:
        flexible_conda_env(config,['conda','rna'],env_yaml = 'envs/rsem.yaml')
    shell:
        r"""
        mkdir -p {REF_B37}/RSEM
        rsem-prepare-reference \
          --gtf {input.gtf} \
          --star -p 20 \
          {input.fasta} \
          {REF_B37}/RSEM/b37
        echo "RSEM reference ready" > {log}
        """


rule build_kallisto_salmon_index:
    conda:
        flexible_conda_env(config,['conda','rna'],env_yaml = 'envs/rsem.yaml')
    input:
        fasta=f"{REF_B37}/Homo_sapiens.GRCh37.cdna.all.fa",
        gz=f"{REF_B37}/Homo_sapiens.GRCh37.cdna.all.fa.gz",
    output:
        kallisto=directory(f"{REF_B37}/kallisto/b37"),
        salmon=directory(f"{REF_B37}/salmon/b37"),
        marker=touch(f"{BUILD_LOG}/kallisto_salmon_index.log"),
    log:
        f"{BUILD_LOG}/kallisto_salmon_index.run.log"
    conda:
        flexible_conda_env(config,['conda','rna'],env_yaml = 'envs/rsem.yaml')
    shell:
        r"""
        mkdir -p {REF_B37}/kallisto/b37
        mkdir -p {REF_B37}/salmon/b37
        kallisto index -i {REF_B37}/kallisto/b37/b37 {input.fasta} -t 20
        salmon index -t {input.gz} -i {REF_B37}/salmon/b37/b37
        echo "kallisto/salmon index ready" > {log}
        """


rule install_clindet_extras:
    input:
        dbsnp_vcf=f"{REF_B37}/Homo_sapiens_assembly19.dbsnp138.vcf",
    output:
        touch(f"{BUILD_LOG}/mass_config.log"),
    log:
        f"{BUILD_LOG}/mass_config.run.log"
    shell:
        r"""
        mkdir -p {BUILD_LOG}
        {{
            echo "Add extra tools to clindet env"
            echo "Create dbsnp bgzip and indel resources"
            bgzip -k -f -o {REF_B37}/Homo_sapiens_assembly19.dbsnp138.vcf.gz {input.dbsnp_vcf}
            tabix -f {REF_B37}/Homo_sapiens_assembly19.dbsnp138.vcf.gz
            bcftools view -v indels --write-index -Oz -o {REF_B37}/Homo_sapiens_assembly19.dbsnp138.indel.vcf.gz {REF_B37}/Homo_sapiens_assembly19.dbsnp138.vcf.gz
            tabix -f {REF_B37}/Homo_sapiens_assembly19.dbsnp138.indel.vcf.gz
            echo "Clone TRUST4"
            if [ ! -d "{SOFT_DIR}/TRUST4" ]; then
                git clone https://github.com/liulab-dfci/TRUST4.git {SOFT_DIR}/TRUST4
            fi
            echo "Download common VCF"
            wget -P {REF_B37} -c https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz
            wget -P {REF_B37} -c https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz.tbi
        }} &> {log}
        touch {output}
        """

rule download_mutation_anno_bed_b37:
    output:
        output=f"{BUILD_LOG}/download_bed_files.log"
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        r"""
        mkdir -p {REF_B37}/bed/giab
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh37@all/Union/GRCh37_alldifficultregions.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh37@all/Union/GRCh37_alllowmapandsegdupregions.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh37@all/FunctionalTechnicallyDifficult/GRCh37_BadPromoters.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh37@all/XY/GRCh37_chrX_PAR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh37@all/XY/GRCh37_chrX_XTR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh37@all/XY/GRCh37_chrY_PAR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh37@all/XY/GRCh37_chrY_XTR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh37@all/OtherDifficult/GRCh37_KIR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh37@all/OtherDifficult/GRCh37_MHC.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh37@all/OtherDifficult/GRCh37_KIR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz
        touch {output}
        """