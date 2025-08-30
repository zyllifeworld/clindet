set -euo pipefail
RED='\e[31m'
RED_B='\e[41m'
NC='\e[0m'
GREEN='\e[32m'
GREEN_B='\e[42m'
bold='\e[1;43m'
mkdir -p build_log

### parallel conda install
cleanup() {
        pkill -P $$
        kill 0
}

for sig in INT QUIT HUP TERM; do
        trap "
            cleanup
            trap - $sig EXIT
            kill -s $sig "'"$$"' "$sig"
done
declare -a pids;

# check conda
command='conda'
if type $command >/dev/null 2>&1; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${command} ${NC} existed, Start download and config."
else
    echo -e "${RED_B}ERROE:  ${NC} ${bold} ${command} ${NC} Not exited !!!!, please install it."
    exit 1 
fi

# check sdk
command='sdk'
if type $command >/dev/null 2>&1; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${command} ${NC} existed, Start download and config."
else
    echo -e "${RED_B}ERROE:  ${NC} ${bold} ${command} ${NC} Not exited !!!!, please install it."
    exit 1 
fi

# check singularity
command='singularity'
if type $command >/dev/null 2>&1; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${command} ${NC} existed, Start download and config."
else
    echo -e "${RED_B}ERROE:  ${NC} ${bold} ${command} ${NC} Not exited !!!!, please install it."
    exit 1 
fi

eval "$(conda shell.bash hook)"

idx=0;
## check conda env and build conda env
ENV_NAME="clindet"
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${ENV_NAME} ${NC} existed."  & pids[$idx]=$!;
	idx=$(($idx+1));
else
    echo -e "${RED_B} ${ENV_NAME} ${NC} not exist，Start build..."
    conda env create -f envs/clindet.yaml -y  &>build_log/${ENV_NAME}.log &  pids[$idx]=$!;
	idx=$(($idx+1));
fi

ENV_NAME="hmftools"
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${ENV_NAME} ${NC} existed." & pids[$idx]=$!;
	idx=$(($idx+1));
else
    echo -e "${RED_B} ${ENV_NAME} ${NC} not exist，Start build..."
    conda env create -f envs/hmftools.yaml -y &>build_log/${ENV_NAME}.log &  pids[$idx]=$!;
	idx=$(($idx+1));
fi

ENV_NAME="clindet_vep"
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${ENV_NAME} ${NC} existed." & pids[$idx]=$!;
	idx=$(($idx+1));
else
    echo -e "${RED_B} ${ENV_NAME} ${NC} not exist，Start build..."
    conda env create -f envs/clindet_vep.yaml -y  &>build_log/${ENV_NAME}.log &  pids[$idx]=$!;
	idx=$(($idx+1));
fi

ENV_NAME="gsutil"
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${ENV_NAME} ${NC} existed." & pids[$idx]=$!;
	idx=$(($idx+1));
else
    echo -e "${RED_B} ${ENV_NAME} ${NC} not exist，Start build..."
    conda env create -f envs/gsutils.yaml -y &>build_log/${ENV_NAME}.log &  pids[$idx]=$!;
	idx=$(($idx+1));
fi

ENV_NAME="strelka"
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${ENV_NAME} ${NC} existed." & pids[$idx]=$!;
	idx=$(($idx+1));
else
    echo -e "${RED_B} ${ENV_NAME} ${NC} not exist，Start build..."
    conda env create -f envs/strelka.yaml -y &>build_log/${ENV_NAME}.log &  pids[$idx]=$!;
	idx=$(($idx+1));
fi

ENV_NAME="clindet_rsem"
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${ENV_NAME} ${NC} existed." & pids[$idx]=$!;
	idx=$(($idx+1));
else
    echo -e "${RED_B} ${ENV_NAME} ${NC} not exist，Start build..."
    conda env create -f envs/rsem.yaml -y &>build_log/${ENV_NAME}.log &  pids[$idx]=$!;
	idx=$(($idx+1));
fi

ENV_NAME="cancer_report"
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${ENV_NAME} ${NC} existed." & pids[$idx]=$!;
	idx=$(($idx+1));
else
    echo -e "${RED_B} ${ENV_NAME} ${NC} not exist，Start build..."
    conda env create -f envs/cancer_report.yaml -y &>build_log/${ENV_NAME}.log &  pids[$idx]=$!;
	idx=$(($idx+1));
fi

ENV_NAME="clindet_multiqc"
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${ENV_NAME} ${NC} existed." & pids[$idx]=$!;
	idx=$(($idx+1));
else
    echo -e "${RED_B} ${ENV_NAME} ${NC} not exist，Start build..."
    conda env create -f envs/multiqc.yaml -y &>build_log/${ENV_NAME}.log &  pids[$idx]=$!;
	idx=$(($idx+1));
fi

### wait install pids
fail=0
for pid in ${pids[*]}; do
        wait $pid || fail=$(( fail + 1 ));
done

if [ "$fail" -gt 0 ]; then
    # 打印红色错误信息并退出
    echo -e "${RED_B} ERROR: Conda build fail ${NC}, See log!" >&2
    exit 1
fi

# echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${ENV_NAME} ${NC} existed." & pids[$idx]=$!;
touch build_log/conda_env_setup.log & echo -e "${GREEN_B} All Conda env built OK ${NC}"

## custome install R package tidyverse/ascat/ for ascat (R >4.3)
# conda activate clindet
# conda install libxml2 pandoc fontconfig freetype -y &>>build_log/clindet.log
# r_mirror="https://cloud.r-project.org" ## can change to nearest location mirror
# R -q -e "install.packages(c('tidyverse','BiocManager'),repos = c(CRAN = '${r_mirror}'))" &>>build_log/clindet.log
# R -q -e 'BiocManager::install(c("devtools", "splines", "readr", "doParallel", "ggplot2", "RColorBrewer", "gridExtra", "gtools", "parallel", "igordot/copynumber", "VariantAnnotation", "GenomicRanges","IRanges"))' &>>build_log/clindet.log
# R -q -e 'devtools::install_github("VanLoo-lab/ascat/ASCAT")' &>>build_log/clindet.log

## consider remove this report section, packages hard to install under "the great wall" !
conda activate cancer_report
echo "Install some R packages of env cancer_report ..."
if [ ! -f "build_log/cancer_report_install_r.log" ]; then
    r_mirror="https://cloud.r-project.org" ## can change to nearest location mirror
    R -q -e "install.packages(c('BiocManager'),repos = c(CRAN = '${r_mirror}'))" &>>build_log/cancer_report.log
    R -q -e 'devtools::install_github("umccr/gpgr")' &>>build_log/cancer_report.log
    R -q -e "install.packages(c('details','DT','kableExtra','patchwork'),repos = c(CRAN = '${r_mirror}'))" &>>build_log/cancer_report.log
    R -q -e 'BiocManager::install("GenomicFeatures")' &>>build_log/cancer_report.log
    # R -q -e 'devtools::install_github("umccr/sigrap")' &>>build_log/cancer_report.log
    echo -e "${GREEN_B} Install custome R packages of env cancer_report finished, continue ....${NC}"
else
    echo -e "${GREEN_B} pull zenodo singularity, continue ${NC}"
fi

# pull image from zenodo/or build by singularity
mkdir -p resources/containers
echo "Beginning singularity image pull ..."
if [ ! -f "build_log/pull_zenodo.log" ]; then
    wget -P resources/containers -c https://zenodo.org/records/15787887/files/pindel.sif
    wget -P resources/containers -c https://zenodo.org/records/15787887/files/brass634.sif
    wget -P resources/containers -c https://zenodo.org/records/15787887/files/caveman153.sif
    wget -P resources/containers -c https://zenodo.org/records/15787887/files/muse230.sif
    wget -P resources/containers -c https://zenodo.org/records/15787887/files/conpair_latest.sif
    wget -P resources/containers -c https://zenodo.org/records/15787887/files/svaba.sif
    touch build_log/pull_zenodo.log
else
    echo -e "${GREEN_B} pull zenodo singularity, continue ${NC}"
fi


conda activate gsutil
mkdir -p resources/ref_genome/b37

# check gsutil
command='gsutil'
if type $command >/dev/null 2>&1; then
    echo -e "${GREEN_B}Check  OK ${NC}, ${bold} ${command} ${NC} existed, Start download and config."
else
    echo -e "${RED_B}ERROE:  ${NC} ${bold} ${command} ${NC} Not exited !!!!"
    exit 1 
fi


## Download b37 genome
echo "Beginning b37 genome fasta Download ..."
if [ ! -f "build_log/download_b37.log" ]; then
    echo "Start download from hmf ..."
    gsutil -m cp -r -n \
      "gs://hmf-public/HMFtools-Resources/ref_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta" \
      "gs://hmf-public/HMFtools-Resources/ref_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict" \
      "gs://hmf-public/HMFtools-Resources/ref_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai" \
      resources/ref_genome/b37
    echo "Done"
    touch build_log/download_b37.log
else
    echo -e "${GREEN_B} already download b37 reference, continue ${NC}"
fi


## Download b37 hmftools softwares data

echo "Beginning b37 hmftools data Download ..."
if [ ! -f "build_log/download_b37_hmftools.log" ]; then
    echo "Start download from hmf ..."
    gsutil -m cp \
      "gs://hmf-public/HMFtools-Resources/pipeline/oncoanalyser/2.0/37/hmf_panel_resources.tso500.37_v2.0.0--3.tar.gz" \
      "gs://hmf-public/HMFtools-Resources/pipeline/oncoanalyser/2.0/37/hmf_pipeline_resources.37_v2.0.0--3.tar.gz" \
      resources/ref_genome/b37
    echo -e "${GREEN_B} Download Done ${NC}"

    wget -P resources/ref_genome/b37/hmf_pipeline_resources -c https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/GermlineHetPon.37.vcf.gz
    wget -P resources/ref_genome/b37/hmf_pipeline_resources -c https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/Amber.snpcheck.37.vcf
    
    mkdir -p resources/ref_genome/b37/hmf_pipeline_resources
    tar -xzvf resources/ref_genome/b37/hmf_panel_resources.tso500.37_v2.0.0--3.tar.gz --strip-components 1 -C resources/ref_genome/b37/hmf_pipeline_resources/
    tar -xzvf resources/ref_genome/b37/hmf_pipeline_resources.37_v2.0.0--3.tar.gz --strip-components 1 -C resources/ref_genome/b37/hmf_pipeline_resources/

    echo -e "${GREEN_B} Decompression Done ${NC}"
    touch build_log/download_b37_hmftools.log
else
    echo -e "${GREEN_B} already download hmftools b37 resources, continue ${NC}"
fi



## Download b37 GATK tools data
echo "Beginning b37 GATK files Download ..."
if [ ! -f "build_log/download_b37_gatk.log" ]; then
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
      resources/ref_genome/b37

    gsutil -m cp -n -r \
        "gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf" \
        "gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf.idx" \
        "gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf" \
        "gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx" \
        "gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf" \
        "gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx" \
      resources/ref_genome/b37

    ## GATK gatk-legacy-bundles always broke, user can try another FTP OR USE gatk-legacy-bundles !!!
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz.md5
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz.md5
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.gz
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.gz.md5
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.gz
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.gz.md5
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase3_v4_20130502.sites.vcf.gz
    # wget -P resources/ref_genome/b37 -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase3_v4_20130502.sites.vcf.gz.md5

    # wget -P resources/ref_genome/b37 -c https://figshare.com/ndownloader/files/45168565 -O 
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
        resources/ref_genome/b37

    echo -e "${GREEN_B} download GATK tools b37 resources Done ${NC}"
    touch build_log/download_b37_gatk.log
else
    echo -e "${GREEN_B} already download GATK tools b37 resources, continue ${NC}"
fi

## Convert gzip to bgzip
echo "Beginning convert b37 GATK gzip files to bgzip ..."
conda activate clindet
if [ ! -f "build_log/b37_convert_to_bgzip.log" ]; then
    gzip -d resources/ref_genome/b37/1000G_phase1.indels.b37.vcf.gz
    bgzip -k -o resources/ref_genome/b37/1000G_phase1.indels.b37.vcf.gz resources/ref_genome/b37/1000G_phase1.indels.b37.vcf
    tabix resources/ref_genome/b37/1000G_phase1.indels.b37.vcf.gz

    echo -e "${GREEN_B} Convert to b37 Done ${NC}"
    touch build_log/b37_convert_to_bgzip.log
else
    echo -e "${GREEN_B} already download GATK tools b37 resources, continue ${NC}"
fi
# Download b37 reference data and softwares config data 


# Download GATK4 softwares
echo "Beginning GATK Download ..."
if [ ! -f "build_log/download_gatk.log" ]; then
    mkdir -p resources/softwares
    ## instal java first this version gatk work on at least java 17.0.15 
    ## sdk install java 17.0.15-tem
    ## sdk default java 17.0.15-tem
    GATK_version="4.6.2.0"
    wget -P resources/softwares -c https://github.com/broadinstitute/gatk/releases/download/${GATK_version}/gatk-${GATK_version}.zip
    unzip resources/softwares/gatk-${GATK_version}.zip -d resources/softwares
    mv resources/softwares/gatk-${GATK_version} resources/softwares/gatk

    echo -e "${GREEN_B} GATK tools Downloaded ${NC}"

    touch build_log/download_gatk.log
else
    echo -e "${GREEN_B} already download GATK tools softwares, continue ${NC}"
fi

## download ASCAT refdata
echo "Beginning ASCAT config files  Download ..."
if [ ! -f "build_log/download_ascat.log" ]; then
    mkdir -p resources/ref_genome/b37/ASCAT/WES
    wget -P resources/ref_genome/b37/ASCAT/WES -c https://zenodo.org/records/14008443/files/G1000_alleles_WES_hg19.zip
    wget -P resources/ref_genome/b37/ASCAT/WES -c https://zenodo.org/records/14008443/files/G1000_loci_WES_hg19.zip
    wget -P resources/ref_genome/b37/ASCAT/WES -c https://zenodo.org/records/14008443/files/GC_G1000_WES_hg19.zip
    wget -P resources/ref_genome/b37/ASCAT/WES -c https://zenodo.org/records/14008443/files/RT_G1000_WES_hg19.zip

    unzip -d resources/ref_genome/b37/ASCAT/WES resources/ref_genome/b37/ASCAT/WES/G1000_alleles_WES_hg19.zip
    unzip -d resources/ref_genome/b37/ASCAT/WES resources/ref_genome/b37/ASCAT/WES/G1000_loci_WES_hg19.zip
    unzip -d resources/ref_genome/b37/ASCAT/WES resources/ref_genome/b37/ASCAT/WES/GC_G1000_WES_hg19.zip
    unzip -d resources/ref_genome/b37/ASCAT/WES resources/ref_genome/b37/ASCAT/WES/RT_G1000_WES_hg19.zip


    wget -P resources/ref_genome/b37/ASCAT/WGS -c https://zenodo.org/records/14008443/files/G1000_alleles_WGS_hg19.zip
    wget -P resources/ref_genome/b37/ASCAT/WGS -c https://zenodo.org/records/14008443/files/G1000_loci_WGS_hg19.zip
    wget -P resources/ref_genome/b37/ASCAT/WGS -c https://zenodo.org/records/14008443/files/GC_G1000_WGS_hg19.zip
    wget -P resources/ref_genome/b37/ASCAT/WGS -c https://zenodo.org/records/14008443/files/RT_G1000_WGS_hg19.zip


    unzip -d resources/ref_genome/b37/ASCAT/WGS resources/ref_genome/b37/ASCAT/WGS/G1000_alleles_WGS_hg19.zip
    unzip -d resources/ref_genome/b37/ASCAT/WGS resources/ref_genome/b37/ASCAT/WGS/G1000_loci_WGS_hg19.zip
    unzip -d resources/ref_genome/b37/ASCAT/WGS resources/ref_genome/b37/ASCAT/WGS/GC_G1000_WGS_hg19.zip
    unzip -d resources/ref_genome/b37/ASCAT/WGS resources/ref_genome/b37/ASCAT/WGS/RT_G1000_WGS_hg19.zip
    echo -e "${GREEN_B} ASCAT configs Downloaded ${NC}"

    touch build_log/download_ascat.log
else
    echo -e "${GREEN_B} already download GATK tools softwares, continue ${NC}"
fi
# Debugging settings
## Download CaVEMan,BRASS config data
echo "Beginning CaVEMan,BRASS config files  Download ..."
if [ ! -f "build_log/download_sanger.log" ]; then
    mkdir -p resources/ref_genome/b37/Sanger
    wget -P resources/ref_genome/b37/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/SNV_INDEL_ref_GRCh37d5-fragment.tar.gz
    wget -P resources/ref_genome/b37/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/CNV_SV_ref_GRCh37d5_brass6+.tar.gz
    wget -P resources/ref_genome/b37/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/core_ref_GRCh37d5.tar.gz
    wget -P resources/ref_genome/b37/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz
    wget -P resources/ref_genome/b37/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/support-files/cgpPindel/cgpPindel_CPBI_refarea.tar.gz
    wget -P resources/ref_genome/b37/Sanger -c https://ftp.sanger.ac.uk/pub/cancer/support-files/CPIB/caveman/cgpCaVEManWrapper_CPBI_refarea.tar.gz


    tar -zxvf resources/ref_genome/b37/Sanger/SNV_INDEL_ref_GRCh37d5-fragment.tar.gz -C resources/ref_genome/b37/Sanger
    tar -zxvf resources/ref_genome/b37/Sanger/CNV_SV_ref_GRCh37d5_brass6+.tar.gz -C resources/ref_genome/b37/Sanger
    tar -zxvf resources/ref_genome/b37/Sanger/core_ref_GRCh37d5.tar.gz -C resources/ref_genome/b37/Sanger
    tar -zxvf resources/ref_genome/b37/Sanger/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz -C resources/ref_genome/b37/Sanger
    tar -zxvf resources/ref_genome/b37/Sanger/cgpPindel_CPBI_refarea.tar.gz -C resources/ref_genome/b37/Sanger
    tar -zxvf resources/ref_genome/b37/Sanger/cgpCaVEManWrapper_CPBI_refarea.tar.gz -C resources/ref_genome/b37/Sanger

    ### copy coding snp and indel to 
    cp resources/ref_genome/b37/Sanger/VAGrENT_ref_GRCh37d5_ensembl_75/vagrent/gene_regions.bed* \
      resources/ref_genome/b37/Sanger/SNV_INDEL_ref/caveman/flagging/

    cp resources/ref_genome/b37/Sanger/VAGrENT_ref_GRCh37d5_ensembl_75/vagrent/codingexon_regions.sub.bed* \
      resources/ref_genome/b37/Sanger/SNV_INDEL_ref/caveman/flagging/

    touch build_log/download_sanger.log
else
    echo -e "${GREEN_B} already download CaVEMan,BRASS, continue ${NC}"
fi



### Download VEP data, version 110
echo "Beginning VEP config files  Download ..."
if [ ! -f "build_log/download_vep.log" ]; then
    mkdir -p resources/ref_genome/b37/vep/v110
    wget -P resources/ref_genome/b37/vep/v110 -c https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh37.tar.gz

    tar -xzvf resources/ref_genome/b37/vep/v110/homo_sapiens_vep_110_GRCh37.tar.gz -C resources/ref_genome/b37/vep/
    touch build_log/download_vep.log
else
    echo -e "${GREEN_B} already download VEP cache ${NC}"
fi

### BWA index
echo "Beginning BWA genome index build ..."
if [ ! -f "resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta.amb" ]; then
    conda activate clindet 
    bwa index resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta
else
    echo -e "${GREEN_B} already built Genome BWA index ${NC}"
fi

### STAR index
echo "Beginning STAR genome index build ..."
if [ ! -f "resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta.amb" ]; then
    conda activate clindet_rsem
    bwa index resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta
else
    echo -e "${GREEN_B} already built Genome BWA index ${NC}"
fi

### do some mass config
conda activate clindet 
echo "Do some mass config ..."
if [ ! -f "build_log/mass_config.log" ]; then
    ### install some package for snakemake slurm and freec
    conda install ucsc-fasplit &>>build_log/clindet.log
    pip install snakemake-executor-plugin-cluster-generic &>>build_log/clindet.log
    pip install configparser &>>build_log/clindet.log

    ### dbsnp bgzip
    bgzip -k -o resources/ref_genome/b37/Homo_sapiens_assembly19.dbsnp138.vcf.gz resources/ref_genome/b37/Homo_sapiens_assembly19.dbsnp138.vcf
    tabix resources/ref_genome/b37/Homo_sapiens_assembly19.dbsnp138.vcf.gz

    ### gatk CreateSequenceDictionary
    resources/softwares/gatk/gatk CreateSequenceDictionary -R resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta
    ### download TRUST4 ref
    git clone https://github.com/liulab-dfci/TRUST4.git resources/softwares/TRUST4
    ### download delly

    ### download varscan2
    ### pull strelka

    ### download common_vcf
    wget -P resources/ref_genome/b37 -c https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz
    wget -P resources/ref_genome/b37 -c https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz.tbi
else

fi

# bwa index resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta

### Config freeC split fasta
# mkdir -p resources/softwares/UCSC_tools
# wget -P resources/softwares/UCSC_tools -c https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSplit
# wget -P resources/ref_genome/b37 -c http://xfer.curie.fr/get/QKFgcU5caZd/hg19_snp137.SingleDiNucl.1based.txt.gz

# chmod 755 resources/softwares/UCSC_tools/faSplit
# mkdir -p resources/ref_genome/b37/fasta
# resources/softwares/UCSC_tools/faSplit byname resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta resources/ref_genome/b37/fasta/
# head -n 22 resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai > resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta.auto.fai

### RSME STAR index
conda activate clindet_rsem
echo "Beginning RSEM star indexing ..."
if [ ! -f "build_log/rsem_star_index.log" ]; then
    wget -P resources/ref_genome/b37/ -c https://ftp.ensembl.org/pub/grch37/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
    gzip -c -d resources/ref_genome/b37/Homo_sapiens.GRCh37.87.gtf.gz > resources/ref_genome/b37/Homo_sapiens.GRCh37.87.gtf
    rsem-prepare-reference \
    --gtf resources/ref_genome/b37/Homo_sapiens.GRCh37.87.gtf \
    --star -p 20 \
    resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta \
    resources/ref_genome/b37/RSEM/b37
    touch build_log/rsem_star_index.log
else
    echo -e "${GREEN_B} already built RSEM STAR index ${NC}"
fi

conda activate clindet_rsem
echo "Beginning STAR indexing ..."
if [ ! -f "build_log/star_index.log" ]; then
    STAR \
    --runThreadN 20 \
    --runMode genomeGenerate \
    --genomeFastaFiles resources/ref_genome/b37/Homo_sapiens.GRCh37.GATK.illumina.fasta \
    --sjdbOverhang 100 --genomeSAindexNbases 2 \
    --sjdbGTFfile resources/ref_genome/b37/Homo_sapiens.GRCh37.87.gtf \
    --genomeDir  resources/ref_genome/b37/STAR/b37 

    touch build_log/star_index.log
else
    echo -e "${GREEN_B} already built STAR index ${NC}"
fi

### build kallisto index
echo "Beginning kallisto & salmon indexing ..."
if [ ! -f "build_log/kallisto_salmon_index.log" ]; then
    wget -P resources/ref_genome/b37/ -c https://ftp.ensembl.org/pub/grch37/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz
    gzip -c -d resources/ref_genome/b37/Homo_sapiens.GRCh37.cdna.all.fa.gz > resources/ref_genome/b37/Homo_sapiens.GRCh37.cdna.all.fa
    mkdir -p resources/ref_genome/b37/kallisto
    kallisto index -i resources/ref_genome/b37/kallisto/b37 resources/ref_genome/b37/Homo_sapiens.GRCh37.cdna.all.fa -t 20
    ### build salmon index
    mkdir -p resources/ref_genome/b37/salmon
    salmon index -t resources/ref_genome/b37/Homo_sapiens.GRCh37.cdna.all.fa.gz -i resources/ref_genome/b37/salmon/b37

    touch build_log/kallisto_salmon_index.log
else
    echo -e "${GREEN_B} already built STAR index ${NC}"
fi

PWD=$(pwd)
echo " \n \n"
echo -e "Now that Clindet has been set up, replace ${RED_B} '/AbsoPath/of/clindet/folder' ${NC}in the workflow/config/config_local_test.yaml file with the absolute path to your Clindet folder."
echo -e "${GREEN_B} ${PWD} ${NC}"

