rule bicseq2_samtools_tumor:
    input:
        bam=""
    output:
        log=""
    conda:
        ""
    threads: 20
    params:
        samtools=,
        out_wd=,
        ref_dir=,
        mapfile_dir=,
    shell:
        """
        ../../scripts/bicseq2/bicseq2_norm.sh {params.out_wd} {params.ref_dif} {params.mapfile_dir} {input.bam} {params.samtools} {threads}
        touch {ouput.log}
        """

rule bicseq2_samtools_normal:
    input:
        bam=""
    output:
        log=""
    conda:
        ""
    threads: 20
    params:
        samtools=,
        out_wd=,
        ref_dir=,
        mapfile_dir=,
    shell:
        """
        ../../scripts/bicseq2/bicseq2_norm.sh {params.out_wd} {params.ref_dif} {params.mapfile_dir} {input.bam} {params.samtools} {threads}
        touch {ouput.log}
        """

rule bicseq2_norm:
    input:
        ""
    output:
        ""
    conda:
        ""
    shell:
        """

        """

rule bicseq2_seg:
    input:
        ""
    output:
        ""
    conda:
        ""
    shell:
        """
        
        """