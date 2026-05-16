rule fastp_trim:
    input:
        unpack(get_rna_fastq)
    output:
        R1=temp("{project}/{genome_version}/results/trimmed/{sample}_R1.fastq.gz"),
        R2=temp("{project}/{genome_version}/results/trimmed/{sample}_R2.fastq.gz"),
        html="{project}/{genome_version}/results/trimmed/fastp/{sample}-fastp.html",
        json="{project}/{genome_version}/results/trimmed/fastp/{sample}-fastp.json"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 2
    conda:
        flexible_conda_env(config,['conda','clindet_mian'],env_yaml = 'envs/clindet.yaml')
    shell:
        """fastp --thread {threads} \
            -i {input.R1} -I {input.R2} \
            -w 16 -Q -c -L \
            -h {output.html} -j {output.json}\
            -o {output.R1} -O {output.R2}
        """


