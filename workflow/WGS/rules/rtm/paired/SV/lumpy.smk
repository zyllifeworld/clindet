rule extract_lumpy_evidence:
    input:
        bam="{project}/{genome_version}/results/recal/paired/{sample}-{group}.bam"
    output:
        discordants=temp(
            "{project}/{genome_version}/results/sv/paired/lumpy/{sample}/evidence/{sample}-{group}.discordants.bam"
        ),
        discordants_bai=temp(
            "{project}/{genome_version}/results/sv/paired/lumpy/{sample}/evidence/{sample}-{group}.discordants.bam.bai"
        ),
        splitters=temp(
            "{project}/{genome_version}/results/sv/paired/lumpy/{sample}/evidence/{sample}-{group}.splitters.bam"
        ),
        splitters_bai=temp(
            "{project}/{genome_version}/results/sv/paired/lumpy/{sample}/evidence/{sample}-{group}.splitters.bam.bai"
        ),
    log:
        "{project}/{genome_version}/logs/sv/paired/lumpy/{sample}/{sample}-{group}.extract_evidence.log"
    benchmark:
        "{project}/{genome_version}/benchmarks/sv/paired/lumpy/{sample}/{sample}-{group}.extract_evidence.tsv"
    threads: 4
    singularity:
        flexible_container_img(config,['singularity','lumpy','sif'],image_url = config['singularity']['lumpy']['repo'])
    shell:
        """
        # Discordant read pairs:
        # -F 1294 removes unmapped, secondary, QC-fail, duplicate etc.
        # -f 1 requires paired reads
        samtools view -@ {threads} -b -F 1294 -f 1 {input.bam} \
          | samtools sort -@ {threads} -o {output.discordants} -

        samtools index -@ {threads} {output.discordants}

        # Split reads:
        # Extract reads with SA tag.
        samtools view -@ {threads} -h {input.bam} \
          | extractSplitReads_BwaMem -i stdin \
          | samtools sort -@ {threads} -o {output.splitters} -

        samtools index -@ {threads} {output.splitters}

        echo "Finished extracting LUMPY evidence for {wildcards.unit}" > {log}
        """

rule lumpy_call_paired:
    input:
        tumour_bam="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        normal_bam="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        tumour_discordants="{project}/{genome_version}/results/sv/paired/lumpy/{sample}/evidence/{sample}-T.discordants.bam",
        tumour_splitters="{project}/{genome_version}/results/sv/paired/lumpy/{sample}/evidence/{sample}-T.splitters.bam",
        normal_discordants="{project}/{genome_version}/results/sv/paired/lumpy/{sample}/evidence/{sample}-NC.discordants.bam",
        normal_splitters="{project}/{genome_version}/results/sv/paired/lumpy/{sample}/evidence/{sample}-NC.splitters.bam",
    output:
        vcf="{project}/{genome_version}/results/sv/paired/lumpy/{sample}/{sample}.lumpy.vcf"
    log:
        "{project}/{genome_version}/logs/sv/paired/lumpy/{sample}/{sample}.lumpy.log"
    benchmark:
        "{project}/{genome_version}/benchmarks/sv/paired/lumpy/{sample}/{sample}.lumpy.tsv"
    threads: 8
    singularity:
        flexible_container_img(config,['singularity','lumpy','sif'],image_url = config['singularity']['lumpy']['repo'])
    shell:
        """
        lumpyexpress \
          -B {input.tumour_bam},{input.normal_bam} \
          -S {input.tumour_splitters},{input.normal_splitters} \
          -D {input.tumour_discordants},{input.normal_discordants} \
          -o {output.vcf} \
          > {log} 2>&1
        """

rule lumpy_svtyper_paired:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/lumpy/{sample}/{sample}.lumpy.vcf",
        tumour_bam="{project}/{genome_version}/results/alignment/paired/{sample}/{sample}-T.dedup.bam",
        tumour_bai="{project}/{genome_version}/results/alignment/paired/{sample}/{sample}-T.dedup.bam.bai",
        normal_bam="{project}/{genome_version}/results/alignment/paired/{sample}/{sample}-NC.dedup.bam",
        normal_bai="{project}/{genome_version}/results/alignment/paired/{sample}/{sample}-NC.dedup.bam.bai",
    output:
        vcf="{project}/{genome_version}/results/sv/paired/lumpy/{sample}/{sample}.lumpy.svtyper.vcf"
    log:
        "{project}/{genome_version}/logs/sv/paired/lumpy/{sample}/{sample}.svtyper.log"
    benchmark:
        "{project}/{genome_version}/benchmarks/sv/paired/lumpy/{sample}/{sample}.svtyper.tsv"
    threads: 4
    singularity:
        flexible_container_img(config,['singularity','lumpy','sif'],image_url = config['singularity']['lumpy']['repo'])
    shell:
        """
        svtyper \
          -i {input.vcf} \
          -B {input.tumour_bam},{input.normal_bam} \
          -o {output.vcf} \
          > {log} 2>&1
        """