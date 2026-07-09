CONTAINER_DIR = "resources/containers"
build_software_log = "build_log/envs"
configfile: "workflow/config/conf/softwares_params.yaml"
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
    return f"{build_software_log}/{env_name}.env.done"

def env_yaml(wildcards):
    return ENV_SPECS[wildcards.env_name]


rule prerequisites:
    output:
        touch(f"{build_software_log}/prereqs.ok")
    log:
        f"{build_software_log}/prereqs.run.log"
    shell:
        r"""
        mkdir -p {build_software_log}
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
        touch(f"{build_software_log}" + "/{env_name}.env.done")
    log:
        f"{build_software_log}" + "/{env_name}.env.install.log"
    params:
        env_name=lambda wc: wc.env_name
    shell:
        r"""
        mkdir -p {build_software_log}
        if conda env list | awk '$1 !~ /^#/ {{print $1}}' | grep -qx "{params.env_name}"; then
            echo "{params.env_name} already exists" > {log}
        else
            conda env create -f {input} -y &> {log}
        fi
        """

rule download_zenodo_containers:
    output:
        touch(f"{build_software_log}/pull_zenodo.log")
    shell:
        r"""
        mkdir -p {CONTAINER_DIR}
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/pindel.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/cgpwgs.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/caveman153.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/brass634.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/control-freec-11.6.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/facets-suite-208.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/igcaller-1.2.1.simg
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/delly_v1.7.2.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/jasminesv.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/svaba.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/arriba240.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/lofreq215.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/muse230.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/deepvariant_deepsomatic-1.9.0.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/gridss2_13_2.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/conpair_latest.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/lancet2.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/ngs-bit-2025.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/moalmanac_latest.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/octopus.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/gatk4_462.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/freebayes.sif
        wget -P {CONTAINER_DIR} -c https://zenodo.org/records/20783116/files/vardict_java.sif
        touch {output}
        """

