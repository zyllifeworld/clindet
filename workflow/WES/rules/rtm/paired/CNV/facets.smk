### PoN of factesCH
rule CNA_snp_pileup_nor:
    input:
        # Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        pileup="{project}/{genome_version}/results/cnv/db/pon_pileups/{sample}_pileup.tsv.gz"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        dbsnp=config['resources'][genome_version]['DBSNP'],
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 2
    shell:
        "snp-pileup -g -A -q15 -r20 -Q20 -P100 -v {params.dbsnp} {output.pileup} {input.NC}"

rule pon_facetsCH:
    input:
        expand("{project}/{genome_version}/results/cnv/db/pon_pileups/{MM_sample}_pileup.tsv.gz",MM_sample = paired_samples,project = project,genome_version = genome_version)
    output:
        log='logs/paired/facets_PoNDB_{sample}.log'
    conda: 'clindet'
    params:
        wd="{project}/{genome_version}/results/cnv/db/pon_pileups"
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GenomicsDBImport -R {params.ref} -L {params.bed} \
        --merge-input-intervals true --sites-only-vcf-output true \
        --genomicsdb-workspace-path /public/ClinicalExam/lj_sih/resource/mutFilter/pon_db \
        {params.vcfs}
        touch {output.log}
        """

### call facetsCH
rule CNA_snp_pileup_tum:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        # NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        pileup="{project}/{genome_version}/results/cnv/db/paired_facetsCH/{sample}/{sample}_pileup.tsv.gz",
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        dbsnp=config['resources'][genome_version]['DBSNP'],
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 2
    shell:
        "snp-pileup -g -A -q15 -r20 -Q20 -P100 -v {params.dbsnp} {output.pileup} {input.Tum}"
