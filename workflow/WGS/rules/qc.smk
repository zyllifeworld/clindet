rule fastp_normal_sample:
    input:
        unpack(get_normal_fastq)
    output:
        R1=temp("{project}/{genome_version}/results/trimmed/{sample}-NC_R1.fastq.gz"),
        R2=temp("{project}/{genome_version}/results/trimmed/{sample}-NC_R2.fastq.gz"),
        html="{project}/{genome_version}/results/trimmed/fastp/{sample}-NC-fastp.html",
        json="{project}/{genome_version}/results/trimmed/fastp/{sample}-NC-fastp.json"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 10
    conda:
        config['softwares']['fastp']['conda']
    shell:
        """fastp --thread {threads} \
            -i {input.R1} -I {input.R2} \
            -w 8 -Q -c -L \
            -h {output.html} -j {output.json} \
            -o {output.R1} -O {output.R2}
        """

rule fastp_tumor_sample:
    input:
        unpack(get_tumor_fastq)
    output:
        R1=temp("{project}/{genome_version}/results/trimmed/{sample}-T_R1.fastq.gz"),
        R2=temp("{project}/{genome_version}/results/trimmed/{sample}-T_R2.fastq.gz"),
        html="{project}/{genome_version}/results/trimmed/fastp/{sample}-T-fastp.html",
        json="{project}/{genome_version}/results/trimmed/fastp/{sample}-T-fastp.json"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 2
    conda:
        config['softwares']['fastp']['conda']
    shell:
        """fastp --thread {threads} \
            -i {input.R1} -I {input.R2} \
            -w 8 -Q -c -L \
            -h {output.html} -j {output.json} \
            -o {output.R1} -O {output.R2}
        """

### check paired Sample Swap with conpair
## GATK3 likely not work with genome sequence file with fasta suffix
conpair_marker_defalut=config['singularity']['conpair'][genome_version]['marker']
rule conpair_pileup:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        Tum_pileup="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T.pileup",
        NC_pileup="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-NC.pileup",
    params:
        ref=config['resources'][genome_version]['REFFA'],
        marker = '' if isinstance(conpair_marker_defalut, bool) and conpair_marker_defalut is True else '-M ' + config['singularity']['conpair'][genome_version]['marker']
    singularity: config['singularity']['conpair']['sif']
    shell:
        """
        /Conpair-0.2/scripts/run_gatk_pileup_for_sample.py -R {params.ref} -B {input.Tum} -O {output.Tum_pileup} {params.marker}
        /Conpair-0.2/scripts/run_gatk_pileup_for_sample.py -R {params.ref} -B {input.NC} -O {output.NC_pileup} {params.marker}
        """

rule conpair_concordance:
    input:
        Tum_pileup="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T.pileup",
        NC_pileup="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-NC.pileup",
    output:
        txt="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T_concordance.txt",
    singularity: config['singularity']['conpair']['sif']
    params:
        marker = '' if isinstance(conpair_marker_defalut, bool) and conpair_marker_defalut is True else '-M ' + config['singularity']['conpair'][genome_version]['marker']
    shell:
        """
        /Conpair-0.2/scripts/verify_concordance.py -T {input.Tum_pileup} -N {input.NC_pileup} --outfile {output.txt}
        """

rule conpair_contamination:
    input:
        Tum_pileup="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T.pileup",
        NC_pileup="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-NC.pileup",
    output:
        txt="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T_contamination.txt",
    params:
        marker = '' if isinstance(conpair_marker_defalut, bool) and conpair_marker_defalut is True else '-M ' + config['singularity']['conpair'][genome_version]['marker']
    singularity: config['singularity']['conpair']['sif']
    shell:
        """
        /Conpair-0.2/scripts/estimate_tumor_normal_contamination.py -T {input.Tum_pileup} -N {input.NC_pileup} --outfile {output.txt} {params.marker}
        """

rule conpair:
    input:
        contam="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T_contamination.txt",
        concord="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T_concordance.txt",
    output:
        touch('{project}/{genome_version}/logs/paired/conpair/{sample}.done')