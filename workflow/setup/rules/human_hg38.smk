import os

REFDIR = "resources/ref_genome/hg38"
LOGDIR = "build_log/hg38"

rule build_hg38_ref:
    input:
        f"{LOGDIR}/download_hg38.log",
        f"{LOGDIR}/download_hg38_hmftools.log",
        f"{LOGDIR}/download_hg38_gatk.log",
        f"{LOGDIR}/download_ascat.log",
        f"{LOGDIR}/download_sanger.log",
        f"{LOGDIR}/download_vep.log",
        f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta.amb",
        f"{LOGDIR}/mass_config.log",
        f"{LOGDIR}/rsem_star_index.log",
        f"{LOGDIR}/star_index.log",
        f"{LOGDIR}/kallisto_salmon_index.log"

rule download_hg38_genome:
    output:
        fasta=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta",
        dictf=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta.dict",
        fai=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta.fai",
        done=f"{LOGDIR}/download_hg38.log"
    shell:
        r"""
        mkdir -p {REFDIR} {LOGDIR}
        gsutil -m cp -r -n \
          gs://hmf-public/HMFtools-Resources/ref_genome/38/GRCh38_masked_exclusions_alts_hlas.fasta \
          gs://hmf-public/HMFtools-Resources/ref_genome/38/GRCh38_masked_exclusions_alts_hlas.fasta.dict \
          gs://hmf-public/HMFtools-Resources/ref_genome/38/GRCh38_masked_exclusions_alts_hlas.fasta.fai \
          {REFDIR}
        touch {output.done}
        """

rule download_hmftools:
    output:
        panel_tar=f"{REFDIR}/hmf_panel_resources.tso500.38_v2.2.0--3.tar.gz",
        pipeline_tar=f"{REFDIR}/hmf_pipeline_resources.38_v2.2.0--3.tar.gz",
        amber1=f"{REFDIR}/hmf_pipeline_resources/GermlineHetPon.hg38.vcf.gz",
        amber2=f"{REFDIR}/hmf_pipeline_resources/GermlineHetPon.hg38.snpcheck.vcf.gz",
        amber3=f"{REFDIR}/hmf_pipeline_resources/Amber.snpcheck.38.vcf",
        done=f"{LOGDIR}/download_hg38_hmftools.log"
    shell:
        r"""
        mkdir -p {REFDIR}/hmf_pipeline_resources {LOGDIR}

        gsutil -m cp -r -n \
          gs://hmf-public/HMFtools-Resources/pipeline/oncoanalyser/2.2/38/hmf_panel_resources.tso500.38_v2.2.0--3.tar.gz \
          gs://hmf-public/HMFtools-Resources/pipeline/oncoanalyser/2.2/38/hmf_pipeline_resources.38_v2.2.0--3.tar.gz \
          {REFDIR}

        wget -P {REFDIR}/hmf_pipeline_resources -c \
          https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/GermlineHetPon.hg38.vcf.gz

        wget -P {REFDIR}/hmf_pipeline_resources -c \
          https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/GermlineHetPon.hg38.snpcheck.vcf.gz

        wget -P {REFDIR}/hmf_pipeline_resources -c \
          https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/Amber.snpcheck.38.vcf

        tar -xzvf {REFDIR}/hmf_panel_resources.tso500.38_v2.2.0--3.tar.gz \
          --strip-components 1 \
          -C {REFDIR}/hmf_pipeline_resources/

        tar -xzvf {REFDIR}/hmf_pipeline_resources.38_v2.2.0--3.tar.gz \
          --strip-components 1 \
          -C {REFDIR}/hmf_pipeline_resources/

        touch {output.done}
        """

rule download_gatk_resources:
    output:
        omni=f"{REFDIR}/1000G_omni2.5.hg38.vcf.gz",
        omni_tbi=f"{REFDIR}/1000G_omni2.5.hg38.vcf.gz.tbi",
        phase1=f"{REFDIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        phase1_tbi=f"{REFDIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
        axiom=f"{REFDIR}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
        axiom_tbi=f"{REFDIR}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi",
        dbsnp_vcf=f"{REFDIR}/Homo_sapiens_assembly38.dbsnp138.vcf",
        dbsnp_idx=f"{REFDIR}/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
        pon=f"{REFDIR}/1000g_pon.hg38.vcf.gz",
        pon_tbi=f"{REFDIR}/1000g_pon.hg38.vcf.gz.tbi",
        gnomad=f"{REFDIR}/af-only-gnomad.hg38.vcf.gz",
        gnomad_tbi=f"{REFDIR}/af-only-gnomad.hg38.vcf.gz.tbi",
        mills=f"{REFDIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        mills_tbi=f"{REFDIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
        dbsnp138=f"{REFDIR}/dbsnp_138.hg38.vcf.gz",
        dbsnp138_tbi=f"{REFDIR}/dbsnp_138.hg38.vcf.gz.tbi",
        common=f"{REFDIR}/00-common_all.vcf.gz",
        common_tbi=f"{REFDIR}/00-common_all.vcf.gz.tbi",
        done=f"{LOGDIR}/download_hg38_gatk.log"
    shell:
        r"""
        mkdir -p {REFDIR} {LOGDIR}

        gsutil -m cp -r -n \
            gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz \
            gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi \
            gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
            gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi \
            gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
            gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi \
            gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf \
            gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx \
            {REFDIR}

        gsutil -m cp -n -r \
            gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz \
            gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi \
            gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
            gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi \
            {REFDIR}

        wget -P {REFDIR} -c \
          http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        wget -P {REFDIR} -c \
          http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
        wget -P {REFDIR} -c \
          http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/dbsnp_138.hg38.vcf.gz
        wget -P {REFDIR} -c \
          http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/dbsnp_138.hg38.vcf.gz.tbi
        wget -P {REFDIR} -c \
          https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz
        wget -P {REFDIR} -c \
          https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz.tbi

        touch {output.done}
        """

rule download_ascat:
    output:
        wes_done=f"{REFDIR}/ASCAT/WES/.done",
        wgs_done=f"{REFDIR}/ASCAT/WGS/.done",
        done=f"{LOGDIR}/download_ascat.log"
    shell:
        r"""
        mkdir -p {REFDIR}/ASCAT/WES {REFDIR}/ASCAT/WGS {LOGDIR}

        wget -P {REFDIR}/ASCAT/WES -c https://zenodo.org/records/14008443/files/G1000_alleles_WES_hg38.zip
        wget -P {REFDIR}/ASCAT/WES -c https://zenodo.org/records/14008443/files/G1000_loci_WES_hg38.zip
        wget -P {REFDIR}/ASCAT/WES -c https://zenodo.org/records/14008443/files/GC_G1000_WES_hg38.zip
        wget -P {REFDIR}/ASCAT/WES -c https://zenodo.org/records/14008443/files/RT_G1000_WES_hg38.zip

        unzip -o -d {REFDIR}/ASCAT/WES {REFDIR}/ASCAT/WES/G1000_alleles_WES_hg38.zip
        unzip -o -d {REFDIR}/ASCAT/WES {REFDIR}/ASCAT/WES/G1000_loci_WES_hg38.zip
        unzip -o -d {REFDIR}/ASCAT/WES {REFDIR}/ASCAT/WES/GC_G1000_WES_hg38.zip
        unzip -o -d {REFDIR}/ASCAT/WES {REFDIR}/ASCAT/WES/RT_G1000_WES_hg38.zip

        for i in $(seq 1 22) X; do
            sed -i 's/^/chr/' {REFDIR}/ASCAT/WES/G1000_lociAll_hg38/G1000_loci_hg38_chr${{i}}.txt
        done

        wget -P {REFDIR}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/G1000_alleles_WGS_hg38.zip
        wget -P {REFDIR}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/G1000_loci_WGS_hg38.zip
        wget -P {REFDIR}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/GC_G1000_WGS_hg38.zip
        wget -P {REFDIR}/ASCAT/WGS -c https://zenodo.org/records/14008443/files/RT_G1000_WGS_hg38.zip

        unzip -o -d {REFDIR}/ASCAT/WGS {REFDIR}/ASCAT/WGS/G1000_alleles_WGS_hg38.zip
        unzip -o -d {REFDIR}/ASCAT/WGS {REFDIR}/ASCAT/WGS/G1000_loci_WGS_hg38.zip
        unzip -o -d {REFDIR}/ASCAT/WGS {REFDIR}/ASCAT/WGS/GC_G1000_WGS_hg38.zip
        unzip -o -d {REFDIR}/ASCAT/WGS {REFDIR}/ASCAT/WGS/RT_G1000_WGS_hg38.zip

        for i in $(seq 1 22) X; do
            sed -i 's/^/chr/' {REFDIR}/ASCAT/WGS/G1000_loci_hg38_chr${{i}}.txt
        done

        touch {output.wes_done} {output.wgs_done} {output.done}
        """

rule download_sanger:
    output:
        flagging_bed=f"{REFDIR}/Sanger/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment/caveman/flagging/gene_regions.bed",
        done=f"{LOGDIR}/download_sanger.log"
    shell:
        r"""
        mkdir -p {REFDIR}/Sanger {LOGDIR}

        wget -P {REFDIR}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz
        wget -P {REFDIR}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz
        wget -P {REFDIR}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz
        wget -P {REFDIR}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz
        wget -P {REFDIR}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/support-files/cgpPindel/cgpPindel_CPBI_refarea.tar.gz
        wget -P {REFDIR}/Sanger -c \
          https://ftp.sanger.ac.uk/pub/cancer/support-files/CPIB/caveman/cgpCaVEManWrapper_CPBI_refarea.tar.gz

        tar -zxvf {REFDIR}/Sanger/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz -C {REFDIR}/Sanger
        tar -zxvf {REFDIR}/Sanger/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz -C {REFDIR}/Sanger
        tar -zxvf {REFDIR}/Sanger/core_ref_GRCh38_hla_decoy_ebv.tar.gz -C {REFDIR}/Sanger
        tar -zxvf {REFDIR}/Sanger/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz -C {REFDIR}/Sanger

        cp {REFDIR}/Sanger/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91/vagrent/gene_regions.bed* \
           {REFDIR}/Sanger/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment/caveman/flagging/

        cp {REFDIR}/Sanger/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91/vagrent/codingexon_regions.sub.bed* \
           {REFDIR}/Sanger/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment/caveman/flagging/

        touch {output.done}
        """

rule download_vep:
    output:
        cache_tar=f"{REFDIR}/vep/v110/homo_sapiens_vep_110_GRCh38.tar.gz",
        done=f"{LOGDIR}/download_vep.log"
    shell:
        r"""
        mkdir -p {REFDIR}/vep/v110 {LOGDIR}
        wget -P {REFDIR}/vep/v110 -c \
          https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
        tar -xzvf {REFDIR}/vep/v110/homo_sapiens_vep_110_GRCh38.tar.gz -C {REFDIR}/vep/
        touch {output.done}
        """

rule bwa_index:
    input:
        f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta"
    output:
        amb=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta.amb",
        ann=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta.ann",
        bwt=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta.bwt",
        pac=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta.pac",
        sa=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta.sa"
    threads: 8
    shell:
        r"""
        bwa index {input}
        """

rule mass_config:
    input:
        fasta=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta",
        dbsnp_gz=f"{REFDIR}/dbsnp_138.hg38.vcf.gz",
        bwa_done=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta.amb"
    output:
        dbsnp_vcf=f"{REFDIR}/dbsnp_138.hg38.vcf",
        dictf=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.dict",
        done=f"{LOGDIR}/mass_config.log"
    shell:
        r"""
        mkdir -p {LOGDIR}
        bgzip -k -d -o {output.dbsnp_vcf} {input.dbsnp_gz}
        resources/softwares/gatk/gatk CreateSequenceDictionary -R {input.fasta}
        touch {output.done}
        """

rule download_gtf:
    output:
        gtf_gz=f"{REFDIR}/gencode.v49.annotation.gtf.gz",
        gtf=f"{REFDIR}/gencode.v49.annotation.gtf"
    shell:
        r"""
        mkdir -p {REFDIR}
        wget -P {REFDIR} -c \
          https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
        gzip -cd {output.gtf_gz} > {output.gtf}
        """

rule download_ensembel_gff3:
    output:
        gff_gz=f"{REFDIR}/Homo_sapiens.GRCh38.115.chr.gff3.gz",
        gff=f"{REFDIR}/Homo_sapiens.GRCh38.115.chr_prefix.gff3"
    shell:
        r"""
        mkdir -p {REFDIR}
        wget -P {REFDIR} -c https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gff3.gz
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
        fasta=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta",
        gtf=f"{REFDIR}/gencode.v49.annotation.gtf",
        gff=f"{REFDIR}/Homo_sapiens.GRCh38.115.chr_prefix.gff3"
    output:
        done=f"{LOGDIR}/rsem_star_index.log"
    threads: 20
    shell:
        r"""
        mkdir -p {REFDIR}/RSEM {LOGDIR}
        rsem-prepare-reference \
            --gtf {input.gtf} \
            --star -p {threads} \
            {input.fasta} \
            {REFDIR}/RSEM/hg38
        touch {output.done}
        """

rule star_index:
    input:
        fasta=f"{REFDIR}/GRCh38_masked_exclusions_alts_hlas.fasta",
        gtf=f"{REFDIR}/gencode.v49.annotation.gtf"
    output:
        done=f"{LOGDIR}/star_index.log"
    threads: 20
    shell:
        r"""
        mkdir -p {REFDIR}/STAR/hg38 {LOGDIR}
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeFastaFiles {input.fasta} \
            --sjdbOverhang 100 \
            --genomeSAindexNbases 2 \
            --sjdbGTFfile {input.gtf} \
            --genomeDir {REFDIR}/STAR/hg38
        touch {output.done}
        """

rule download_cdna:
    output:
        fa_gz=f"{REFDIR}/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        fa=f"{REFDIR}/Homo_sapiens.GRCh38.cdna.all.fa"
    shell:
        r"""
        mkdir -p {REFDIR}
        wget -P {REFDIR} -c \
          https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        gzip -cd {output.fa_gz} > {output.fa}
        """

rule kallisto_salmon_index:
    input:
        fa=f"{REFDIR}/Homo_sapiens.GRCh38.cdna.all.fa",
        fa_gz=f"{REFDIR}/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    output:
        done=f"{LOGDIR}/kallisto_salmon_index.log"
    threads: 20
    shell:
        r"""
        mkdir -p {REFDIR}/kallisto {REFDIR}/salmon {LOGDIR}

        kallisto index \
            -i {REFDIR}/kallisto/hg38 \
            {input.fa} \
            -t {threads}

        salmon index \
            -t {input.fa_gz} \
            -i {REFDIR}/salmon/hg38

        touch {output.done}
        """