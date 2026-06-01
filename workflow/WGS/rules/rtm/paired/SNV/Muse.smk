rule muse_call:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}.MuSE.txt"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}",
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','muse','sif'],image_url = config['singularity']['muse']['repo'])
    shell:
        """
        if [ -x /MuSE/bin/MuSE ]; then
            MUSE_BIN=/MuSE/bin/MuSE
        elif command -v MuSE >/dev/null 2>&1; then
            MUSE_BIN=MuSE
        else
            echo "ERROR: Cannot find MuSE executable. Tried /MuSE/bin/MuSE and MuSE in PATH." >&2
            exit 1
        fi

        $MUSE_BIN call -f {params.ref} -O {params.out_prefix} -n {threads} {input.Tum} {input.NC}
        """

rule muse_sump:
    input:
        txt="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}.MuSE.txt"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/muse.vcf"
    params:
        dbsnp=config['resources'][genome_version]['DBSNP_GZ']
    threads: 10
    singularity: 
        flexible_container_img(config,['singularity','muse','sif'],image_url = config['singularity']['muse']['repo'])
    shell:
        """
        if [ -x /MuSE/bin/MuSE ]; then
            MUSE_BIN=/MuSE/bin/MuSE
        elif command -v MuSE >/dev/null 2>&1; then
            MUSE_BIN=MuSE
        else
            echo "ERROR: Cannot find MuSE executable. Tried /MuSE/bin/MuSE and MuSE in PATH." >&2
            exit 1
        fi

        $MUSE_BIN sump -I {input.txt} -O {output.vcf} -G -n {threads} -D {params.dbsnp}
        """
