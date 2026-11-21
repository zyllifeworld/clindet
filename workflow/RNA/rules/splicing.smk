##
## Alternative-splicing quantification adapted from TCASIA_pipeline
## (workflows/02_as_calling).
##
## Upstream reuse:
##   - STAR_splice_map: dedicated STAR alignment with TCASIA parameters
##   - splicing_salmon_quant: Salmon against the splicing transcriptome index
##   - splicing_detect_strandness: auto-detect from Salmon meta_info.json
##   - splicing_read_length: fixed or samtools-detected read length
##
## Downstream callers (per sample):
##   - rMATS, MAJIQ/Voila, SUPPA2, SplAdder
##
## Outputs:
##   {project}/{genome_version}/results/AS/{sample}/...
##
from pathlib import Path

from snakemake.exceptions import WorkflowError

# The SUPPA2 event file is shared by all samples and has no {sample} wildcard.
# Run it on the head node so the Slurm profile's {wildcards.sample} log path is
# only evaluated for per-sample rules.
localrules: splicing_suppa2_events

_SPLICE_GV = config["project"]["genome_version"]
_SPLICE_CFG = (
    config.get("softwares_params", {})
    .get(_SPLICE_GV, {})
    .get("splicing")
)

if not _SPLICE_CFG:
    raise WorkflowError(
        f"splicing stage requested but "
        f"softwares_params.{_SPLICE_GV}.splicing is missing"
    )

AS_RESOURCE_DIR = (
    "{project}/{genome_version}/results/AS/resources/{genome_version}"
)


def _splice_cfg(wildcards):
    shared = (
        config.get("softwares_params", {})
        .get(wildcards.genome_version, {})
        .get("splicing", {})
    )
    # Per-batch overrides, e.g. use GENCODE v49 for one test batch while the
    # shared config keeps v50.
    run_cfg = config.get("run_params", {}).get("splicing", {}) or {}
    if not run_cfg:
        return shared
    merged = dict(shared)
    merged.update(run_cfg)
    return merged


def _tool_cfg(wildcards, tool):
    return _splice_cfg(wildcards).get(tool) or {}


def get_splice_gtf(wildcards, tool):
    """Tool-specific GTF override, then splicing.gtf, then genome GTF."""
    return (
        _tool_cfg(wildcards, tool).get("gtf")
        or _splice_cfg(wildcards).get("gtf")
        or config["resources"][wildcards.genome_version]["GTF"]
    )


def get_splice_alignment_gtf(wildcards):
    return (
        _splice_cfg(wildcards).get("gtf")
        or config["resources"][wildcards.genome_version]["GTF"]
    )


def get_splice_gff(wildcards, tool="majiq"):
    return (
        _tool_cfg(wildcards, tool).get("gff")
        or _splice_cfg(wildcards).get("gff")
        or config["resources"][wildcards.genome_version]["GFF"]
    )


def get_splice_star_index(wildcards):
    return (
        _splice_cfg(wildcards).get("star_index")
        or config["softwares_params"][wildcards.genome_version]["star"]["index"]
    )


def get_splice_salmon_index(wildcards, tool="suppa2"):
    return (
        _tool_cfg(wildcards, tool).get("salmon_index")
        or _splice_cfg(wildcards).get("salmon_index")
        or config["softwares_params"][wildcards.genome_version]["salmon"]["index"]
    )


def get_splice_strandness_mode(wildcards):
    return _splice_cfg(wildcards).get("strandness", "auto")


def get_splice_read_len_mode(wildcards):
    return _splice_cfg(wildcards).get("read_len_mode", "fixed")


def get_splice_read_len(wildcards):
    return int(_splice_cfg(wildcards).get("read_len", 150))


def get_splice_layout(wildcards):
    r2 = samples_info.loc[wildcards.sample, "R2_file_path"]
    return "PE" if is_valid_file_path(r2) else "SE"


def get_splice_star_gtf_option(wildcards):
    use_gtf = _splice_cfg(wildcards).get("star_use_gtf_at_align", True)
    return "true" if use_gtf else "false"


def get_suppa_events(wildcards):
    configured = (_tool_cfg(wildcards, "suppa2").get("events") or "").strip()
    if configured:
        return configured
    return f"{AS_RESOURCE_DIR}/suppa2_events.ioe".format(**wildcards)


# ---------------------------------------------------------------------------
# Alignment and per-sample inputs
# ---------------------------------------------------------------------------

rule STAR_splice_map:
    input:
        unpack(get_rna_fastq)
    output:
        bam="{project}/{genome_version}/results/mapped/STAR_splice/{sample}/{sample}.sorted.bam",
        bai="{project}/{genome_version}/results/mapped/STAR_splice/{sample}/{sample}.sorted.bam.bai",
    params:
        out_dir=(
            "{project}/{genome_version}/results/mapped/STAR_splice/{sample}/"
        ),
        star_index=lambda wc: get_splice_star_index(wc),
        gtf=lambda wc: get_splice_alignment_gtf(wc),
        use_gtf=lambda wc: get_splice_star_gtf_option(wc),
        read_files_in=get_star_read_files_in,
        layout=get_rna_layout,
        limit_bam_sort_ram=39050942993,
    threads: 10
    log:
        "{project}/{genome_version}/logs/AS/{sample}.star_splice.log"
    conda:
        flexible_conda_env(config, ["conda", "rna"], env_yaml="envs/rsem.yaml")
    shell:
        """
        mkdir -p {params.out_dir} "$(dirname {log})"

        GTF_OPT=""
        if [ "{params.use_gtf}" = "true" ]; then
            GTF_OPT="--sjdbGTFfile {params.gtf}"
        fi

        STAR \
            --runThreadN {threads} \
            --genomeDir "{params.star_index}" \
            --readFilesIn {params.read_files_in} \
            --outFileNamePrefix {params.out_dir}{wildcards.sample}_ \
            --outSAMunmapped Within \
            --readFilesCommand zcat \
            --outFilterMismatchNmax 15 \
            --outSAMattributes All \
            --outStd Log \
            --limitBAMsortRAM {params.limit_bam_sort_ram} \
            --outSAMtype BAM Unsorted \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMstrandField intronMotif \
            --twopassMode Basic \
            --quantMode GeneCounts \
            $GTF_OPT \
            >> {log} 2>&1

        samtools sort \
            -@ {threads} \
            -m 1G \
            -o {output.bam} \
            {params.out_dir}{wildcards.sample}_Aligned.out.bam \
            >> {log} 2>&1
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        rm -f {params.out_dir}{wildcards.sample}_Aligned.out.bam
        """


rule splicing_salmon_quant:
    input:
        unpack(get_rna_fastq)
    output:
        quant="{project}/{genome_version}/results/AS/{sample}/salmon/quant.sf"
    params:
        index=lambda wc: get_splice_salmon_index(wc, "suppa2"),
        reads=get_salmon_reads,
        out_dir="{project}/{genome_version}/results/AS/{sample}/salmon",
    threads: 8
    log:
        "{project}/{genome_version}/logs/AS/{sample}.salmon.log"
    conda:
        flexible_conda_env(config, ["conda", "rna"], env_yaml="envs/rsem.yaml")
    shell:
        """
        mkdir -p "$(dirname {log})"
        rm -rf {params.out_dir}
        mkdir -p {params.out_dir}
        salmon quant \
            -i "{params.index}" \
            -l A \
            --gcBias \
            --validateMappings \
            {params.reads} \
            -p {threads} \
            -o {params.out_dir} \
            > {log} 2>&1
        """


rule splicing_detect_strandness:
    input:
        quant="{project}/{genome_version}/results/AS/{sample}/salmon/quant.sf"
    output:
        txt="{project}/{genome_version}/results/AS/{sample}/strandness.txt"
    params:
        mode=lambda wc: get_splice_strandness_mode(wc),
        default="fr-firststrand",
    conda:
        flexible_conda_env(config, ["conda", "clindet_main"], env_yaml="envs/clindet.yaml")
    script:
        "../scripts/detect_strandness.py"


rule splicing_read_length:
    input:
        bam="{project}/{genome_version}/results/mapped/STAR_splice/{sample}/{sample}.sorted.bam",
    output:
        txt="{project}/{genome_version}/results/AS/{sample}/read_length.txt"
    params:
        mode=lambda wc: get_splice_read_len_mode(wc),
        fixed=lambda wc: get_splice_read_len(wc),
    conda:
        flexible_conda_env(config, ["conda", "clindet_main"], env_yaml="envs/clindet.yaml")
    shell:
        """
        mkdir -p "$(dirname {output.txt})"
        if [ "{params.mode}" = "detect" ]; then
            samtools stats {input.bam} \
                | awk -F '\t' '$1 == "SN" && $2 == "maximum length:" {{print $3}}' \
                > {output.txt}
            [ -s {output.txt} ] || echo {params.fixed} > {output.txt}
        else
            echo {params.fixed} > {output.txt}
        fi
        """


# ---------------------------------------------------------------------------
# rMATS (single-sample PSI, statistics off)
# ---------------------------------------------------------------------------

rule splicing_rmats:
    input:
        bam="{project}/{genome_version}/results/mapped/STAR_splice/{sample}/{sample}.sorted.bam",
        strandness="{project}/{genome_version}/results/AS/{sample}/strandness.txt",
        read_len="{project}/{genome_version}/results/AS/{sample}/read_length.txt",
    output:
        results=directory(
            "{project}/{genome_version}/results/AS/{sample}/rmats"
        ),
    params:
        gtf=lambda wc: get_splice_gtf(wc, "rmats"),
        cstat=lambda wc: float(_tool_cfg(wc, "rmats").get("cstat", 0.0001)),
        extra=lambda wc: _tool_cfg(wc, "rmats").get("extra", ""),
        tmp="{project}/{genome_version}/results/AS/{sample}/work/rmats",
        layout=lambda wc: get_splice_layout(wc),
    threads: 10
    log:
        "{project}/{genome_version}/logs/AS/{sample}.rmats.log"
    singularity:
        flexible_container_img(config, ["singularity", "splicing_rmats", "sif"])
    shell:
        """
        mkdir -p {params.tmp} "$(dirname {log})"
        rm -rf {output.results}
        mkdir -p {output.results}

        printf '%s\n' {input.bam} > {params.tmp}/bam_list.txt
        STRAND=$(cat {input.strandness})
        READLEN=$(cat {input.read_len})
        if [ "{params.layout}" = "PE" ]; then
            RNA_TYPE=paired
        else
            RNA_TYPE=single
        fi

        rmats.py \
            --b1 {params.tmp}/bam_list.txt \
            --gtf "{params.gtf}" \
            --readLength "$READLEN" \
            --nthread {threads} \
            --od {output.results} \
            --tmp {params.tmp} \
            --libType "$STRAND" \
            --cstat {params.cstat} \
            --statoff \
            --variable-read-length \
            -t "$RNA_TYPE" {params.extra} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# MAJIQ / Voila
# ---------------------------------------------------------------------------

rule splicing_majiq_ini:
    input:
        bam="{project}/{genome_version}/results/mapped/STAR_splice/{sample}/{sample}.sorted.bam",
        strandness="{project}/{genome_version}/results/AS/{sample}/strandness.txt",
        read_len="{project}/{genome_version}/results/AS/{sample}/read_length.txt",
    output:
        ini="{project}/{genome_version}/results/AS/{sample}/majiq/build.ini",
    params:
        genome=lambda wc: _tool_cfg(wc, "majiq").get("genome", "hg38"),
        bam_dir=lambda wc: str(Path(
            f"{wc.project}/{wc.genome_version}/results/mapped/STAR_splice/{wc.sample}"
        )),
        bam_stem=lambda wc: f"{wc.sample}.sorted",
    conda:
        flexible_conda_env(config, ["conda", "clindet_main"], env_yaml="envs/clindet.yaml")
    shell:
        """
        mkdir -p "$(dirname {output.ini})"
        READLEN=$(cat {input.read_len})
        STRAND=$(cat {input.strandness})
        case "$STRAND" in
            fr-firststrand)  MAJIQ_STRAND=reverse ;;
            fr-secondstrand) MAJIQ_STRAND=forward ;;
            *)               MAJIQ_STRAND=none ;;
        esac
        printf '[info]\nreadlen=%s\nbamdirs=%s\ngenome=%s\nstrandness=%s\n[experiments]\n%s=%s\n' \
            "$READLEN" {params.bam_dir} {params.genome} "$MAJIQ_STRAND" \
            {wildcards.sample} {params.bam_stem} > {output.ini}
        """


rule splicing_majiq_build:
    input:
        ini="{project}/{genome_version}/results/AS/{sample}/majiq/build.ini",
    output:
        majiq="{project}/{genome_version}/results/AS/{sample}/majiq/{sample}.sorted.majiq",
        sj="{project}/{genome_version}/results/AS/{sample}/majiq/{sample}.sorted.sj",
        sql="{project}/{genome_version}/results/AS/{sample}/majiq/splicegraph.sql",
    params:
        gff=lambda wc: get_splice_gff(wc, "majiq"),
        license=lambda wc: _tool_cfg(wc, "majiq").get("license", ""),
        minreads=lambda wc: int(_tool_cfg(wc, "majiq").get("minreads", 10)),
        output_dir="{project}/{genome_version}/results/AS/{sample}/majiq",
    threads: 8
    log:
        "{project}/{genome_version}/logs/AS/{sample}.majiq_build.log"
    singularity:
        flexible_container_img(config, ["singularity", "splicing_majiq", "sif"])
    shell:
        """
        mkdir -p {params.output_dir} "$(dirname {log})"
        majiq --license "{params.license}" build "{params.gff}" \
            -c "{input.ini}" -j {threads} -o "{params.output_dir}" \
            --minreads {params.minreads} > {log} 2>&1
        """


rule splicing_majiq_psi:
    input:
        majiq="{project}/{genome_version}/results/AS/{sample}/majiq/{sample}.sorted.majiq",
    output:
        psi="{project}/{genome_version}/results/AS/{sample}/majiq/{sample}.psi.tsv",
        voila="{project}/{genome_version}/results/AS/{sample}/majiq/{sample}.psi.voila",
    params:
        license=lambda wc: _tool_cfg(wc, "majiq").get("license", ""),
        output_dir="{project}/{genome_version}/results/AS/{sample}/majiq",
    threads: 8
    log:
        "{project}/{genome_version}/logs/AS/{sample}.majiq_psi.log"
    singularity:
        flexible_container_img(config, ["singularity", "splicing_majiq", "sif"])
    shell:
        """
        mkdir -p {params.output_dir} "$(dirname {log})"
        majiq --license "{params.license}" psi "{input.majiq}" \
            -n {wildcards.sample} -j {threads} -o "{params.output_dir}" \
            > {log} 2>&1
        """


rule splicing_voila_modulize:
    input:
        sql="{project}/{genome_version}/results/AS/{sample}/majiq/splicegraph.sql",
        voila="{project}/{genome_version}/results/AS/{sample}/majiq/{sample}.psi.voila",
    output:
        modules=directory(
            "{project}/{genome_version}/results/AS/{sample}/majiq/modulized"
        ),
    params:
        license=lambda wc: _tool_cfg(wc, "majiq").get("license", ""),
    threads: 4
    log:
        "{project}/{genome_version}/logs/AS/{sample}.voila_modulize.log"
    singularity:
        flexible_container_img(config, ["singularity", "splicing_majiq", "sif"])
    shell:
        """
        mkdir -p "$(dirname {log})"
        rm -rf {output.modules}
        voila --license "{params.license}" modulize -d {output.modules} \
            {input.sql} {input.voila} -j {threads} > {log} 2>&1
        """


rule splicing_voila_tsv:
    input:
        sql="{project}/{genome_version}/results/AS/{sample}/majiq/splicegraph.sql",
        voila="{project}/{genome_version}/results/AS/{sample}/majiq/{sample}.psi.voila",
    output:
        result="{project}/{genome_version}/results/AS/{sample}/majiq/{sample}.voila.tsv",
    params:
        license=lambda wc: _tool_cfg(wc, "majiq").get("license", ""),
    threads: 4
    log:
        "{project}/{genome_version}/logs/AS/{sample}.voila_tsv.log"
    singularity:
        flexible_container_img(config, ["singularity", "splicing_majiq", "sif"])
    shell:
        """
        mkdir -p "$(dirname {output.result})" "$(dirname {log})"
        voila --license "{params.license}" tsv -f {output.result} \
            {input.sql} {input.voila} -j {threads} > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# SUPPA2 (reuse the splicing Salmon quant, generate events from GTF)
# ---------------------------------------------------------------------------

rule splicing_suppa2_events:
    input:
        gtf=lambda wc: get_splice_gtf(wc, "suppa2"),
    output:
        events="{project}/{genome_version}/results/AS/resources/{genome_version}/suppa2_events.ioe",
    params:
        prefix="{project}/{genome_version}/results/AS/resources/{genome_version}/suppa2_events",
    log:
        "{project}/{genome_version}/logs/AS/suppa2_events.{genome_version}.log"
    singularity:
        flexible_container_img(config, ["singularity", "splicing_suppa2", "sif"])
    shell:
        """
        mkdir -p "$(dirname {output.events})" "$(dirname {log})"
        rm -f {params.prefix}_*.ioe
        suppa.py generateEvents \
            -i "{input.gtf}" \
            -o "{params.prefix}" \
            -e SE SS MX RI FL \
            -f ioe > {log} 2>&1
        # Keep exactly one header line; each generated .ioe has its own header.
        first=$(ls {params.prefix}_*.ioe | head -n 1)
        head -n 1 "$first" > {output.events}
        grep -h -v '^seqname' {params.prefix}_*.ioe >> {output.events}
        """


rule splicing_suppa2_tpm:
    input:
        quant="{project}/{genome_version}/results/AS/{sample}/salmon/quant.sf",
    output:
        tpm="{project}/{genome_version}/results/AS/{sample}/suppa2/iso_tpm_formatted.txt",
    conda:
        flexible_conda_env(config, ["conda", "clindet_main"], env_yaml="envs/clindet.yaml")
    script:
        "../scripts/select_suppa_fields.py"


rule splicing_suppa2:
    input:
        tpm="{project}/{genome_version}/results/AS/{sample}/suppa2/iso_tpm_formatted.txt",
        events=lambda wc: get_suppa_events(wc),
    output:
        psi="{project}/{genome_version}/results/AS/{sample}/suppa2/{sample}_event.psi",
    params:
        min_tpm=lambda wc: float(_tool_cfg(wc, "suppa2").get("min_tpm", 1)),
        prefix="{project}/{genome_version}/results/AS/{sample}/suppa2/{sample}_event",
    log:
        "{project}/{genome_version}/logs/AS/{sample}.suppa2.log"
    singularity:
        flexible_container_img(config, ["singularity", "splicing_suppa2", "sif"])
    shell:
        """
        mkdir -p "$(dirname {output.psi})" "$(dirname {log})"
        suppa.py psiPerEvent \
            -i "{input.events}" \
            -f {params.min_tpm} \
            -e "{input.tpm}" \
            -o "{params.prefix}" > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# SplAdder
# ---------------------------------------------------------------------------

rule splicing_spladder:
    input:
        bam="{project}/{genome_version}/results/mapped/STAR_splice/{sample}/{sample}.sorted.bam",
    output:
        results=directory(
            "{project}/{genome_version}/results/AS/{sample}/spladder"
        ),
    params:
        annotation=lambda wc: get_splice_gtf(wc, "spladder"),
        confidence=lambda wc: int(_tool_cfg(wc, "spladder").get("confidence", 3)),
        merge_strategy=lambda wc: _tool_cfg(wc, "spladder").get("merge_strategy", "single"),
        event_types=lambda wc: _tool_cfg(wc, "spladder").get(
            "event_types",
            "exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons",
        ),
    threads: 8
    log:
        "{project}/{genome_version}/logs/AS/{sample}.spladder.log"
    conda:
        flexible_conda_env(config, ["conda", "splicing_spladder"], env_yaml="envs/spladder.yaml")
    shell:
        """
        rm -rf {output.results}
        mkdir -p {output.results} "$(dirname {log})"
        spladder build \
            --bams "{input.bam}" \
            --annotation "{params.annotation}" \
            --outdir {output.results} \
            --parallel {threads} \
            --merge-strat {params.merge_strategy} \
            --event-types {params.event_types} \
            --confidence {params.confidence} -v > {log} 2>&1
        """
