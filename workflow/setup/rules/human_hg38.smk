import os

REF_Hg38 = "resources/ref_genome/hg38"
BUILD_LOG_hg38 = "build_log/hg38"

rule build_hg38_ref:
    input:
        f"{BUILD_LOG_hg38}/download_hg38.log",
        f"{BUILD_LOG_hg38}/download_hg38_hmftools.log",
        f"{BUILD_LOG_hg38}/download_hg38_gatk.log",
        f"{BUILD_LOG_hg38}/download_ascat.log",
        f"{BUILD_LOG_hg38}/download_sanger.log",
        f"{BUILD_LOG_hg38}/download_vep.log",
        f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta.amb",
        f"{BUILD_LOG_hg38}/mass_config.log",
        f"{BUILD_LOG_hg38}/rsem_star_index.log",
        f"{BUILD_LOG_hg38}/star_index.log",
        f"{BUILD_LOG_hg38}/kallisto_salmon_index.log",
        f"{REF_Hg38}/Homo_sapiens.GRCh38.115.chr_prefix.gff3"

rule download_hg38_genome:
    output:
        fasta=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta",
        dictf=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta.dict",
        fai=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta.fai",
        done=f"{BUILD_LOG_hg38}/download_hg38.log"
    shell:
        r"""
        mkdir -p {REF_Hg38} {BUILD_LOG_hg38}
        gsutil -m cp -r -n \
          gs://hmf-public/HMFtools-Resources/ref_genome/38/GRCh38_masked_exclusions_alts_hlas.fasta \
          gs://hmf-public/HMFtools-Resources/ref_genome/38/GRCh38_masked_exclusions_alts_hlas.fasta.dict \
          gs://hmf-public/HMFtools-Resources/ref_genome/38/GRCh38_masked_exclusions_alts_hlas.fasta.fai \
          {REF_Hg38}
        touch {output.done}
        """

rule download_hmftools:
    output:
        panel_tar=f"{REF_Hg38}/hmf_panel_resources.tso500.38_v2.2.0--3.tar.gz",
        pipeline_tar=f"{REF_Hg38}/hmf_pipeline_resources.38_v2.2.0--3.tar.gz",
        amber1=f"{REF_Hg38}/hmf_pipeline_resources/GermlineHetPon.hg38.vcf.gz",
        amber2=f"{REF_Hg38}/hmf_pipeline_resources/GermlineHetPon.hg38.snpcheck.vcf.gz",
        amber3=f"{REF_Hg38}/hmf_pipeline_resources/Amber.snpcheck.38.vcf",
        done=f"{BUILD_LOG_hg38}/download_hg38_hmftools.log"
    shell:
        r"""
        mkdir -p {REF_Hg38}/hmf_pipeline_resources {BUILD_LOG_hg38}

        gsutil -m cp -r -n \
          gs://hmf-public/HMFtools-Resources/pipeline/oncoanalyser/2.2/38/hmf_panel_resources.tso500.38_v2.2.0--3.tar.gz \
          gs://hmf-public/HMFtools-Resources/pipeline/oncoanalyser/2.2/38/hmf_pipeline_resources.38_v2.2.0--3.tar.gz \
          {REF_Hg38}

        wget -P {REF_Hg38}/hmf_pipeline_resources -c \
          https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/GermlineHetPon.hg38.vcf.gz

        wget -P {REF_Hg38}/hmf_pipeline_resources -c \
          https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/GermlineHetPon.hg38.snpcheck.vcf.gz

        wget -P {REF_Hg38}/hmf_pipeline_resources -c \
          https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/Amber.snpcheck.38.vcf

        tar -xzvf {REF_Hg38}/hmf_panel_resources.tso500.38_v2.2.0--3.tar.gz \
          --strip-components 1 \
          -C {REF_Hg38}/hmf_pipeline_resources/

        tar -xzvf {REF_Hg38}/hmf_pipeline_resources.38_v2.2.0--3.tar.gz \
          --strip-components 1 \
          -C {REF_Hg38}/hmf_pipeline_resources/

        touch {output.done}
        """

rule download_gatk_resources:
    output:
        omni=f"{REF_Hg38}/1000G_omni2.5.hg38.vcf.gz",
        omni_tbi=f"{REF_Hg38}/1000G_omni2.5.hg38.vcf.gz.tbi",
        phase1=f"{REF_Hg38}/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        phase1_tbi=f"{REF_Hg38}/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
        axiom=f"{REF_Hg38}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
        axiom_tbi=f"{REF_Hg38}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi",
        dbsnp_vcf=f"{REF_Hg38}/Homo_sapiens_assembly38.dbsnp138.vcf",
        dbsnp_idx=f"{REF_Hg38}/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
        pon=f"{REF_Hg38}/1000g_pon.hg38.vcf.gz",
        pon_tbi=f"{REF_Hg38}/1000g_pon.hg38.vcf.gz.tbi",
        gnomad=f"{REF_Hg38}/af-only-gnomad.hg38.vcf.gz",
        gnomad_tbi=f"{REF_Hg38}/af-only-gnomad.hg38.vcf.gz.tbi",
        mills=f"{REF_Hg38}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        mills_tbi=f"{REF_Hg38}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
        dbsnp138=f"{REF_Hg38}/dbsnp_138.hg38.vcf.gz",
        dbsnp138_tbi=f"{REF_Hg38}/dbsnp_138.hg38.vcf.gz.tbi",
        common=f"{REF_Hg38}/00-common_all.vcf.gz",
        common_tbi=f"{REF_Hg38}/00-common_all.vcf.gz.tbi",
        done=f"{BUILD_LOG_hg38}/download_hg38_gatk.log"
    shell:
        r"""
        mkdir -p {REF_Hg38} {BUILD_LOG_hg38}

        gsutil -m cp -r -n \
            gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz \
            gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi \
            gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
            gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi \
            gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
            gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi \
            gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf \
            gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx \
            {REF_Hg38}

        gsutil -m cp -n -r \
            gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz \
            gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi \
            gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
            gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi \
            {REF_Hg38}

        wget -P {REF_Hg38} -c \
          http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        wget -P {REF_Hg38} -c \
          http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
        wget -P {REF_Hg38} -c \
          http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/dbsnp_138.hg38.vcf.gz
        wget -P {REF_Hg38} -c \
          http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/dbsnp_138.hg38.vcf.gz.tbi
        wget -P {REF_Hg38} -c \
          https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz
        wget -P {REF_Hg38} -c \
          https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz.tbi

        touch {output.done}
        """

rule download_ascat:
    output:
        wes_done=f"{REF_Hg38}/ASCAT/WES/.done",
        wgs_done=f"{REF_Hg38}/ASCAT/WGS/.done",
        done=f"{BUILD_LOG_hg38}/download_ascat.log"
    shell:
        r"""
        mkdir -p {REF_Hg38}/ASCAT/WES {REF_Hg38}/ASCAT/WGS {BUILD_LOG_hg38}

        wget -P {REF_Hg38}/ASCAT/WES -c https://zenodo.org/records/14008443/files/G1000_alleles_WES_hg38.zip
        wget -P {REF_Hg38}/ASCAT/WES -c https://zenodo.org/records/14008443/files/G1000_loci_WES_hg38.zip
        wget -P {REF_Hg38}/ASCAT/WES -c https://zenodo.org/records/14008443/files/GC_G1000_WES_hg38.zip
        wget -P {REF_Hg38}/ASCAT/WES -c https://zenodo.org/records/14008443/files/RT_G1000_WES_hg38.zip

        unzip -o -d {REF_Hg38}/ASCAT/WES {REF_Hg38}/ASCAT/WES/G1000_alleles_WES_hg38.zip
        unzip -o -d {REF_Hg38}/ASCAT/WES {REF_Hg38}/ASCAT/WES/G1000_loci_WES_hg38.zip
        unzip -o -d {REF_Hg38}/ASCAT/WES {REF_Hg38}/ASCAT/WES/GC_G1000_WES_hg38.zip
        unzip -o -d {REF_Hg38}/ASCAT/WES {REF_Hg38}/ASCAT/WES/RT_G1000_WES_hg38.zip

        for i in $(seq 1 22) X; do
            sed -i 's/^/chr/' {REF_Hg38}/ASCAT/WES/G1000_lociAll_hg38/G1000_loci_hg38_chr${{i}}.txt
        done

        wget -P {REF_Hg38}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/G1000_alleles_WGS_hg38.zip
        wget -P {REF_Hg38}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/G1000_loci_WGS_hg38.zip
        wget -P {REF_Hg38}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/GC_G1000_WGS_hg38.zip
        wget -P {REF_Hg38}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/RT_G1000_WGS_hg38.zip

        unzip -o -d {REF_Hg38}/ASCAT/WGS {REF_Hg38}/ASCAT/WGS/G1000_alleles_WGS_hg38.zip
        unzip -o -d {REF_Hg38}/ASCAT/WGS {REF_Hg38}/ASCAT/WGS/G1000_loci_WGS_hg38.zip
        unzip -o -d {REF_Hg38}/ASCAT/WGS {REF_Hg38}/ASCAT/WGS/GC_G1000_WGS_hg38.zip
        unzip -o -d {REF_Hg38}/ASCAT/WGS {REF_Hg38}/ASCAT/WGS/RT_G1000_WGS_hg38.zip

        for i in $(seq 1 22) X; do
            sed -i 's/^/chr/' {REF_Hg38}/ASCAT/WGS/G1000_loci_hg38_chr${{i}}.txt
        done

        touch {output.wes_done} {output.wgs_done} {output.done}
        """

rule download_sanger:
    output:
        flagging_bed=f"{REF_Hg38}/Sanger/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment/caveman/flagging/gene_regions.bed",
        done=f"{BUILD_LOG_hg38}/download_sanger.log"
    shell:
        r"""
        mkdir -p {REF_Hg38}/Sanger {BUILD_LOG_hg38}

        wget -P {REF_Hg38}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz
        wget -P {REF_Hg38}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz
        wget -P {REF_Hg38}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz
        wget -P {REF_Hg38}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz
        wget -P {REF_Hg38}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/support-files/cgpPindel/cgpPindel_CPBI_refarea.tar.gz
        wget -P {REF_Hg38}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/support-files/CPIB/caveman/cgpCaVEManWrapper_CPBI_refarea.tar.gz

        tar -zxvf {REF_Hg38}/Sanger/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz -C {REF_Hg38}/Sanger
        tar -zxvf {REF_Hg38}/Sanger/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz -C {REF_Hg38}/Sanger
        tar -zxvf {REF_Hg38}/Sanger/core_ref_GRCh38_hla_decoy_ebv.tar.gz -C {REF_Hg38}/Sanger
        tar -zxvf {REF_Hg38}/Sanger/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz -C {REF_Hg38}/Sanger

        cp {REF_Hg38}/Sanger/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91/vagrent/gene_regions.bed* \
           {REF_Hg38}/Sanger/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment/caveman/flagging/

        cp {REF_Hg38}/Sanger/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91/vagrent/codingexon_regions.sub.bed* \
           {REF_Hg38}/Sanger/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment/caveman/flagging/

        touch {output.done}
        """

rule download_vep_cache_hg38:
    output:
        touch(f"{BUILD_LOG_hg38}/download_vep.log")
    log:
        f"{BUILD_LOG_hg38}/download_vep.run.log"
    params:
        cache_num=str(config["softwares_params"]['hg38']["vcf2maf"]["vep"]["cache_version"]).replace("v", ""),
        cache_dir=f"{REF_Hg38}/vep/v{str(config['softwares_params']['hg38']['vcf2maf']['vep']['cache_version']).replace('v', '')}",
        vep_dir=f"{REF_Hg38}/vep"
    shell:
        """
        (
            set -euo pipefail

            mkdir -p {params.cache_dir}

            wget \
                -P {params.cache_dir} \
                -c https://ftp.ensembl.org/pub/release-{params.cache_num}/variation/indexed_vep_cache/homo_sapiens_vep_{params.cache_num}_GRCh38.tar.gz

            tar -xzvf \
                {params.cache_dir}/homo_sapiens_vep_{params.cache_num}_GRCh37.tar.gz \
                -C {params.vep_dir}
        ) &> {log}

        touch {output}
        """

rule bwa_index:
    input:
        f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta"
    output:
        amb=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta.amb",
        ann=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta.ann",
        bwt=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta.bwt",
        pac=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta.pac",
        sa=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta.sa"
    threads: 8
    shell:
        r"""
        bwa index {input}
        """

rule mass_config:
    input:
        fasta=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta",
        dbsnp_gz=f"{REF_Hg38}/dbsnp_138.hg38.vcf.gz",
        bwa_done=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta.amb"
    output:
        dbsnp_vcf=f"{REF_Hg38}/dbsnp_138.hg38.vcf",
        dictf=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.dict",
        done=f"{BUILD_LOG_hg38}/mass_config.log"
    shell:
        r"""
        mkdir -p {BUILD_LOG_hg38}
        bgzip -k -d -o {output.dbsnp_vcf} {input.dbsnp_gz}
        bcftools view -v indels --write-index -Oz -o {REF_Hg38}/dbsnp_138.hg38.indel.vcf.gz {input.dbsnp_gz}
        tabix -f {REF_Hg38}/dbsnp_138.hg38.indel.vcf.gz
        echo "Clone TRUST4"
        if [ ! -d "{SOFT_DIR}/TRUST4" ]; then
            git clone https://github.com/liulab-dfci/TRUST4.git {SOFT_DIR}/TRUST4
        fi
        touch {output.done}
        """

rule download_gtf:
    output:
        gtf_gz=f"{REF_Hg38}/gencode.v49.annotation.gtf.gz",
        gtf=f"{REF_Hg38}/gencode.v49.annotation.gtf"
    shell:
        r"""
        mkdir -p {REF_Hg38}
        wget -P {REF_Hg38} -c \
          https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
        gzip -cd {output.gtf_gz} > {output.gtf}
        """

rule download_ensembel_gff3:
    output:
        gff_gz=f"{REF_Hg38}/Homo_sapiens.GRCh38.115.chr.gff3.gz",
        gff=f"{REF_Hg38}/Homo_sapiens.GRCh38.115.chr_prefix.gff3"
    shell:
        r"""
        mkdir -p {REF_Hg38}
        wget -P {REF_Hg38} -c https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gff3.gz
        zcat {output.gff_gz} | \
        awk 'BEGIN{FS=OFS="\t"}
            /^#/ {print; next}
            {
              if ($1 ~ /^(1[0-9]|2[0-2]|[1-9]|X|Y)$/) $1="chr"$1;
              else if ($1=="MT") $1="chrM";
              print
            }' > {output.gff}
        """

rule rsem_star_index:
    input:
        fasta=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta",
        gtf=f"{REF_Hg38}/gencode.v49.annotation.gtf",
        gff=f"{REF_Hg38}/Homo_sapiens.GRCh38.115.chr_prefix.gff3"
    output:
        done=f"{BUILD_LOG_hg38}/rsem_star_index.log"
    threads: 20
    shell:
        r"""
        mkdir -p {REF_Hg38}/RSEM {BUILD_LOG_hg38}
        rsem-prepare-reference \
            --gtf {input.gtf} \
            --star -p {threads} \
            {input.fasta} \
            {REF_Hg38}/RSEM/hg38
        touch {output.done}
        """

rule star_index:
    input:
        fasta=f"{REF_Hg38}/GRCh38_masked_exclusions_alts_hlas.fasta",
        gtf=f"{REF_Hg38}/gencode.v49.annotation.gtf"
    output:
        done=f"{BUILD_LOG_hg38}/star_index.log"
    threads: 20
    shell:
        r"""
        mkdir -p {REF_Hg38}/STAR/hg38 {BUILD_LOG_hg38}
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeFastaFiles {input.fasta} \
            --sjdbOverhang 100 \
            --genomeSAindexNbases 2 \
            --sjdbGTFfile {input.gtf} \
            --genomeDir {REF_Hg38}/STAR/hg38
        touch {output.done}
        """

rule download_cdna:
    output:
        fa_gz=f"{REF_Hg38}/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        fa=f"{REF_Hg38}/Homo_sapiens.GRCh38.cdna.all.fa"
    shell:
        r"""
        mkdir -p {REF_Hg38}
        wget -P {REF_Hg38} -c \
          https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        gzip -cd {output.fa_gz} > {output.fa}
        """

rule kallisto_salmon_index:
    input:
        fa=f"{REF_Hg38}/Homo_sapiens.GRCh38.cdna.all.fa",
        fa_gz=f"{REF_Hg38}/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    output:
        done=f"{BUILD_LOG_hg38}/kallisto_salmon_index.log"
    threads: 20
    shell:
        r"""
        mkdir -p {REF_Hg38}/kallisto {REF_Hg38}/salmon {BUILD_LOG_hg38}

        kallisto index \
            -i {REF_Hg38}/kallisto/hg38 \
            {input.fa} \
            -t {threads}

        salmon index \
            -t {input.fa_gz} \
            -i {REF_Hg38}/salmon/hg38

        touch {output.done}
        """

rule download_mutation_anno_bed_hg38:
    output:
        output=f"{BUILD_LOG_hg38}/download_bed_files.log"
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        r"""
        mkdir -p {REF_B37}/bed/giab
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/Union/GRCh38_alllowmapandsegdupregions.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/FunctionalTechnicallyDifficult/GRCh38_BadPromoters.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/XY/GRCh38_chrX_PAR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/XY/GRCh38_chrX_XTR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/XY/GRCh38_chrY_PAR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/XY/GRCh38_chrY_XTR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/OtherDifficult/GRCh38_KIR.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/OtherDifficult/GRCh38_MHC.bed.gz
        wget -P {REF_B37}/bed/giab -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/OtherDifficult/GRCh38_KIR.bed.gz

        wget -P {REF_B37}/bed/giab -c https://zenodo.org/records/16755940/files/hg38.pm151b-v3.easy.bed.gz
        wget -P {REF_B37}/bed/giab -c https://zenodo.org/records/16755940/files/hg38.longdust.bed.gz
        wget https://github.com/parklab/SMaHT_Regional_Categorization/raw/main/SMaHT_easy_hg38.bed.gz
        wget https://github.com/parklab/SMaHT_Regional_Categorization/raw/main/SMaHT_difficult_hg38.bed.gz
        wget https://github.com/parklab/SMaHT_Regional_Categorization/raw/main/SMaHT_extreme_hg38.bed.gz

        touch {output}
        """