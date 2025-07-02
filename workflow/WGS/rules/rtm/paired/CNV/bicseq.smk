rule bicseq:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
    params:
        # ASCAT should change config file because allelCounter need chr prefix in hg19 version
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        lociprefix=config['softwares']['ascat'][genome_version]['loci_1000'],
        allelesprefix=config['softwares']['ascat'][genome_version]['alleles_1000'],
        GCcontentfile=config['softwares']['ascat'][genome_version]['GCcontentfile'],
        replictimingfile=config['softwares']['ascat'][genome_version]['replictimingfile'],
        # sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        """
        $bam_processing_lib\/bicseq_preprocess.sh $tumor_bam $ref $norm_map cn_norm &>bicseq_preprocess.log &
        """