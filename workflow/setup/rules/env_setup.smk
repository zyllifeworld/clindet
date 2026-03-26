SOFT_DIR = "resources/softwares"
CONTAINER_DIR = "resources/containers"

ENV_SPECS = {
    "clindet": "envs/clindet.yaml",
    "hmftools": "envs/hmftools.yaml",
    "clindet_vep": "envs/clindet_vep.yaml",
    "gsutil": "envs/gsutils.yaml",
    "strelka": "envs/strelka.yaml",
    "clindet_rsem": "envs/rsem.yaml",
    "cancer_report": "envs/cancer_report.yaml",
    "clindet_multiqc": "envs/multiqc.yaml",
}


def env_done(env_name):
    return f"{BUILD_LOG}/{env_name}.env.done"


def env_yaml(wildcards):
    return ENV_SPECS[wildcards.env_name]
rule append_analysis_envs:
    input:
        f"{BUILD_LOG}/prereqs.ok",
        expand(f"{BUILD_LOG}/{{env_name}}.env.done", env_name=sorted(ENV_SPECS)),
        f"{BUILD_LOG}/cancer_report_install_r.log",
        f"{BUILD_LOG}/pull_zenodo.log"

rule prerequisites:
    output:
        touch(f"{BUILD_LOG}/prereqs.ok")
    log:
        f"{BUILD_LOG}/prereqs.run.log"
    shell:
        r"""
        mkdir -p {BUILD_LOG}
        for cmd in conda sdk singularity; do
            if ! command -v "$cmd" >/dev/null 2>&1; then
                echo "ERROR: $cmd is not installed or not on PATH." >&2
                exit 1
            fi
        done
        echo "Prerequisites OK" > {log}
        """

rule create_conda_env:
    input:
        env_yaml
    output:
        touch(f"{BUILD_LOG}" + "/{env_name}.env.done")
    log:
        f"{BUILD_LOG}" + "/{env_name}.env.install.log"
    params:
        env_name=lambda wc: wc.env_name
    shell:
        r"""
        mkdir -p {BUILD_LOG}
        if conda env list | awk '$1 !~ /^#/ {{print $1}}' | grep -qx "{params.env_name}"; then
            echo "{params.env_name} already exists" > {log}
        else
            conda env create -f {input} -y &> {log}
        fi
        """


rule install_cancer_report_r:
    input:
        env_done("cancer_report")
    output:
        touch(f"{BUILD_LOG}/cancer_report_install_r.log")
    log:
        f"{BUILD_LOG}/cancer_report_install_r.run.log"
    shell:
        r"""
        mkdir -p {BUILD_LOG}
        {{
            echo "Install custom R packages for cancer_report"
            conda run -n cancer_report R -q -e "install.packages(c('BiocManager'), repos = c(CRAN = 'https://cloud.r-project.org'))"
            conda run -n cancer_report R -q -e 'devtools::install_github("umccr/gpgr")'
            conda run -n cancer_report R -q -e "install.packages(c('details','DT','kableExtra','patchwork'), repos = c(CRAN = 'https://cloud.r-project.org'))"
            conda run -n cancer_report R -q -e 'BiocManager::install("GenomicFeatures")'
        }} &> {log}
        """


rule download_zenodo_containers:
    output:
        touch(f"{BUILD_LOG}/pull_zenodo.log")
    shell:
        r"""
        mkdir -p {CONTAINER_DIR}
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/15787887/files/pindel.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/15787887/files/brass634.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/15787887/files/caveman153.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/15787887/files/muse230.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/15787887/files/conpair_latest.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/15787887/files/svaba.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/16892396/files/deepsomatic_160.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/17963718/files/delly_v1.7.2.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/17963718/files/facets-suite-dev.img
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/17963718/files/jasminesv.sif
        touch {output}
        """

