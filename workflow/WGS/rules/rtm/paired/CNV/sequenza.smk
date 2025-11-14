rule sequenza_gc_bins:
    input:
        ref=config['resources'][genome_version]['REFFA'],
    output:
        gc="{project}/{genome_version}/results/cnv/paired/sequenza/{genome_version}.gc50.wig.gz"
    threads: 8
    conda:
        config['softwares']['sequenza']['conda']
    shell:
        """
        sequenza-utils gc_wiggle -w 50 \
        --fasta {input.ref} \
        --output {output.gc}
        """

rule sequenza_bam2seqz:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        gc="{project}/{genome_version}/results/cnv/paired/sequenza/{genome_version}.gc50.wig.gz"
    output:
        seqz="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}.seqz.gz",
    threads: 8
    params:
        gc="{project}/{genome_version}/results/cnv/paired/sequenza/{genome_version}.gc50.wig.gz"
    conda:
        config['softwares']['sequenza']['conda']
    shell:
        """
        sequenza-utils bam2seqz \
        --normal {input.NC} \
        --tumor {input.Tum} \
        --fasta {input.ref} -gc {params.gc} \
        --output {output.seqz}
        """

rule sequenza_seqz_binning:
    input:
        seqz="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}.seqz.gz"
    output:
        bin_seqz="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}.bin50_seqz.gz",
    conda:
        config['softwares']['sequenza']['conda']
    shell:
        """
        sequenza-utils seqz_binning -w 50 --seqz {input.seqz} -o {output.bin_seqz}
        """


rule sequenza_call:
    input:
        bin_seqz="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}.bin50_seqz.gz",
    output:
        segment="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}/{sample}_segments.txt",
    params:
        wd="{project}/{genome_version}/results/cnv/paired/sequenza/{sample}",
    script:
        "../../../../scripts/sequenza.R"