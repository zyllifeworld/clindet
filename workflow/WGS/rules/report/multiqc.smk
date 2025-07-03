import os
from os.path import abspath, join, dirname, basename

# prep purple/conpair/samtools/fastp/gatk/picard results for multiqc
rule prep_multiqc_data:
    input:
        fastp_tum = "{project}/{genome_version}/results/trimmed/fastp/{sample}-T-fastp.json",
        fastp_nc = "{project}/{genome_version}/results/trimmed/fastp/{sample}-NC-fastp.json",
        conpair_concord = "{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T_concordance.txt",
        conpair_contam = "{project}/{genome_version}/results/qc/conpair/paired/{sample}/{sample}-T_contamination.txt",
        purple_stats = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.purity.tsv',
        purple_qc = '{project}/{genome_version}/results/cnv/paired/purple/{sample}/purple/{sample}.purple.qc',
        dedup_metrics = "{project}/{genome_version}/results/qc/dedup/{paired}/{sample}-T.metrics.txt",
        picard = "{project}/{genome_version}/results/stats/{paired}/wes_metrics/{sample}-T.txt"
    output:
        filelist                = '{project}/{genome_version}/results/multiqc_data/{sample}/filelist.txt',
        renamed_file_dir        = directory('{project}/{genome_version}/results/multiqc_data/{sample}'),
    params:
        data_dir                = '{project}/{genome_version}/results/multiqc_data/{sample}'
    group: 'multiqc'
    run:
        if not os.path.exists(output.renamed_file_dir):
            os.makedirs(output.renamed_file_dir)
        qc_files = []
        # qc results
        # change conpair result name
        renamed_conpair_concord    = join(params.data_dir, basename(input.conpair_concord) \
            .replace('-T_', '.'))
        renamed_conpair_contam = join(params.data_dir, basename(input.conpair_contam) \
            .replace('-T_', '.'))
        shell(f'cp {input.conpair_concord} {renamed_conpair_concord}')
        shell(f'cp {input.conpair_contam} {renamed_conpair_contam}')
        qc_files.extend([
            renamed_conpair_concord,
            renamed_conpair_contam,
        ])

        qc_files.extend([
            input.fastp_tum,
            input.fastp_nc,
            input.purple_stats,
            input.purple_qc,
            input.dedup_metrics,
            input.picard,
        ])

        try:
            with open(output.filelist, 'w') as out:
                for fp in qc_files:
                    if fp:
                        print(fp, file=out)
        except OSError as e:
            err(str(e))

### multiqc on conpair/fastp/
rule combined_multiqc_prep_multiqc_data:
    input:
        filelists=expand('{project}/{genome_version}/results/multiqc_data/{sample}/filelist.txt',project = project,genome_version = genome_version,sample = [*unpaired_samples, *paired_samples]),
    output:
        filelist='{project}/{genome_version}/results/multiqc/filelist.txt',
    group: 'combined_multiqc'
    run:
        qc_files = []
        for filelist in input.filelists:
            with open(filelist) as f:
                for l in f:
                    l = l.strip()
                    if l and l not in qc_files:
                        qc_files.append(l)
        # output final filelist
        try:
            with open(output.filelist, 'w') as out:
                for fp in qc_files:
                    if fp:
                        print(fp, file=out)
        except OSError as e:
            err(str(e))


# rule combined_multiqc:
#     input:
#         umccrise_conf_yaml  = join(package_path(), 'multiqc', 'multiqc_config.yaml'),
#         filelist            = 'multiqc/multiqc_data/filelist.txt',
#         generated_conf_yaml = 'multiqc/multiqc_data/generated_conf.yaml',
#     output:
#         html_file           = 'multiqc/multiqc_report.html'
#     group: 'combined_multiqc'
#     run:
#         list_files = input.filelist
#         shell(f'LC_ALL=$LC_ALL LANG=$LANG multiqc -f -o . -l {list_files}'
#                 f' -c {input.umccrise_conf_yaml} -c {input.generated_conf_yaml} --filename {output.html_file}')

