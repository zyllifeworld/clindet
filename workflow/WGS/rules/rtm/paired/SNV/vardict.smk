rule vardict_wgs_bed:
    input:
        reference=config['resources'][genome_version]['REFFA']
    output:
        genome_bed="{project}/{genome_version}/results/vcf/paired/{genome_version}.bed"
    shell:
        """
        awk 'BEGIN {FS="\t"}; {print $1 "\t0\t" $2}' {input.reference}.fai > {output.genome_bed}
        """

rule vardict_paired_mode:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions="{project}/{genome_version}/results/vcf/paired/{genome_version}.bed",
        bam="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        normal="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict_raw.vcf"
    params:
        extra="",
        bed_columns="-c 1 -S 2 -E 3 -g 4",  # Optional, default is -c 1 -S 2 -E 3 -g 4
        allele_frequency_threshold="0.01",  # Optional, default is 0.01
        post_scripts="testsomatic.R | var2vcf_paired.pl -N "
    threads: 10
    conda: config['conda']['clindet_main']
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.vardict.benchmark.txt"
    shell:
        """
        vardict-java -G {input.reference} \
        {params.extra} \
        -th {threads} \
        {params.bed_columns} \
        -f {params.allele_frequency_threshold} -N '{wildcards.sample}_T|{wildcards.sample}_NC' \
        -b "{input.bam}|{input.normal}" \
        {input.regions} | {params.post_scripts} '{wildcards.sample}_T|{wildcards.sample}_NC' -f {params.allele_frequency_threshold} > {output.vcf}
        """

rule vardict_filter_somatic:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict_raw.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict.vcf"
    threads: 1
    params:
        caller='vardict'
    conda: config['conda']['clindet_main']
    shell:
        """
        bcftools view -i '(INFO/STATUS~"StrongSomatic" || INFO/STATUS~"LikelySomatic") && FILTER="PASS" && INFO/SSF <= 0.05' {input.vcf} > {output.vcf} 
        """

rule vardict_filter_germline:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict_raw.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf_germline/paired/{sample}/vardict_germline.vcf"
    threads: 1
    params:
        caller='vardict'
    conda: config['conda']['clindet_main']
    shell:
        """
        bcftools view -i 'INFO/STATUS~"Germline" && FILTER="PASS"' {input.vcf} > {output.vcf} 
        """