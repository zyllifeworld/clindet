rule sequenza_bam2seqz:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        gc="{project}/{genome_version}/genome/{genome_version}.gc50Base.txt.gz"
    output:
        seqz="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}.seqz.gz",
    threads: 8
    conda:
        config['softwares']['sequenza']['conda']
    shell:
        """
        sequenza-utils bam2seqz \
        --normal {input.NC} \
        --tumor {input.Tum} \
        --fasta {input.ref} -gc {input.gc} \
        --output {output.seqz}
        """

rule sequenza_seqz_binning:
    input:
        seqz="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}.seqz.gz"
    output:
        bin_seqz="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}.bin50_seqz.gz",
    conda:
        config['softwares']['sequenza']['conda']
    shell:
        """
        sequenza-utils seqz_binning -w 50 --seqz {input.seqz} -o {output.bin_seqz}
        """


rule sequenza_call:
    input:
        bin_seqz="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}.bin50_seqz.gz",
    output:
        segment="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}_segments.txt",
    params:
        wd="{project}/{genome_version}/results/cnv/sequenza/{sample}",
    script:
        "../../../../scripts/sequenza.R"