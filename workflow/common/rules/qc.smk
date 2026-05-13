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
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mapping/{sample}_NC.fastp.benchmark.txt"
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
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
    threads: 10
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mapping/{sample}_T.fastp.benchmark.txt"
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """fastp --thread {threads} \
            -i {input.R1} -I {input.R2} \
            -w 8 -Q -c -L \
            -h {output.html} -j {output.json} \
            -o {output.R1} -O {output.R2}
        """

### check paired Sample Swap with conpair
## GATK3 likely not work with genome sequence file with fasta suffix
## unable while not human data
conpair_config = config['singularity']['conpair'].get(genome_version, False)
if conpair_config:
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
        singularity: 
            flexible_container_img(config,['singularity','conpair','sif'],image_url = config['singularity']['conpair']['repo'])
        benchmark:
            "{project}/{genome_version}/results/benchmarks/conpair/{sample}.pileup.benchmark.txt"
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
        singularity:
            flexible_container_img(config,['singularity','conpair','sif'],image_url = config['singularity']['conpair']['repo'])
        params:
            marker = '' if isinstance(conpair_marker_defalut, bool) and conpair_marker_defalut is True else '-M ' + config['singularity']['conpair'][genome_version]['marker']
        benchmark:
            "{project}/{genome_version}/results/benchmarks/conpair/{sample}.concordance.benchmark.txt"
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
        singularity: 
            flexible_container_img(config,['singularity','conpair','sif'],image_url = config['singularity']['conpair']['repo'])
        benchmark:
            "{project}/{genome_version}/results/benchmarks/conpair/{sample}.contamination.benchmark.txt"
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

# samtools flagstat
rule bam_flagstat:
    input:
        bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
    output:
        stat="{project}/{genome_version}/results/stats/{sample_type}/wgs_metrics/{sample}-{group}.bam.flagstat",
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    threads:10
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mapping/{sample_type}/{sample}-{group}.flagstat.benchmark.txt"
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.stat}
        """

##ngs
rule ngs_bit_sample_gender:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        tsv="{project}/{genome_version}/results/qc/conpair/qc/ngs_bit/{sample}-gender.tsv",
    params:
        build_version = get_ngsbit_build
    singularity: 
        flexible_container_img(config,['singularity','ngs_bit','sif'],image_url = config['singularity']['ngs_bit']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/ngs-bit/{sample}.ngsbit_gender.benchmark.txt"
    shell:
        """
        SampleGender -in {input.NC} -method xy  -out {output.tsv}
        """

rule ngs_bit_mapping:
    input:
        bam="{project}/{genome_version}/results/recal/paired/{sample}-{sample_group}.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        txt="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}_{sample_group}_mapqc.txt",
    singularity: 
        flexible_container_img(config,['singularity','ngs_bit','sif'],image_url = config['singularity']['ngs_bit']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/ngs-bit/{sample}.ngsbit_mappingqc_{sample_group}.benchmark.txt"
    shell:
        """
        MappingQC -in {input.bam} \
        -roi {input.bed} \
        -txt \
        -ref {input.ref} \
        -out {output.txt}
        """

# rule discvrseq_variant_qc:
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#     output:
#         txt="{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T_contamination.txt",
#     params:
#         build_version = '' if isinstance(conpair_marker_defalut, bool) and conpair_marker_defalut is True else '-M ' + config['singularity']['conpair'][genome_version]['marker']
#     singularity: 
#         flexible_container_img(config,['singularity','discvrseq','sif'],image_url = config['singularity']['discvrseq']['repo'])
#     benchmark:
#         "{project}/{genome_version}/results/benchmarks/qc/{sample}.discvrseq_variant_qc.benchmark.txt"
#     shell:
#         """
#         java -jar /DISCVRSeq.jar  VariantQC \
#         -R {input.ref} \
#         -V {input.vcf} \
#         -O {output.html} \
#         -rd {output.json}
#         """
