rule esvee_prep:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/sv/paired/esvee/{sample}/SV_gridss_{sample}.vcf",
    params:
        ref=config['resources'][genome_version]['REFFA'],
        wd="{project}/{genome_version}/results/sv/paired/esvee/{sample}",
        known_fusion_bed=config['singularity']['hmftools'][genome_version]['esvee']['known_fusion_bed'],
        ref_version=config['singularity']['hmftools'][genome_version]['esvee']['ref_genome_version']
        # if else like  (exprs ? c1 : c2) in C++/javascript is OK too,eg on next line
        # blacklist=("" if config['singularity']['gridss'][genome_version]['blacklist'] == "" else " -b " + config['singularity']['gridss'][genome_version]['blacklist'])
    singularity: config['singularity']['gridss']['sif']
    threads:16
    shell:
        """
            java -cp ~/softwares/hmftools/esvee_v1.0-rc.5.jar  \
            com.hartwig.hmftools.esvee.prep.PrepApplication \
            -sample '{wildcards.sample},{wildcards.sample}_NC' \
            -bam_file '{input.Tum},{input.NC}' \
            -ref_genome {input.ref} \
            -ref_genome_version {params.ref_version} \
            -known_fusion_bed {params.known_fusion_bed} \
            -bamtool ~/.conda/envs/snake/bin/sambamba \
            -output_dir {params.wd} \
            -threads {threads}
        """

# java  -Xmx120g  -cp ~/softwares/hmftools/esvee_v1.0-rc.5.jar              com.hartwig.hmftools.esvee.prep.PrepApplication             -sample 'COL0829,COL0829_NC'             -bam_file 'WGS/b37/results/recal/paired/COL0829-T.bam,WGS/b37/results/recal/paired/COL0829-NC.bam'             -ref_genome /public/ClinicalExam/lj_sih/projects/project_clindet/reference/b37/Homo_sapiens_assembly19.fasta             -ref_genome_version 37             -known_fusion_bed /public/ClinicalExam/lj_sih/projects/project_clindet/reference/b37/pipeline/v5_34/ref/37/sv/known_fusions.37.bedpe             -bamtool ~/.conda/envs/snake/bin/sambamba             -output_dir WGS/b37/results/sv/paired/esvee/COL0829             -threads 16

# rule esvee_ref_depth:
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#         ref=config['resources'][genome_version]['REFFA'],
#     output:
#         vcf="{project}/{genome_version}/results/sv/paired/esvee/{sample}/SV_gridss_{sample}.vcf",
#     params:
#         ref=config['resources'][genome_version]['REFFA'],
#         wd="{project}/{genome_version}/results/sv/paired/esvee/{sample}",
#         blacklist=get_gridss_blacklist
#         # if else like  (exprs ? c1 : c2) in C++/javascript is OK too,eg on next line
#         # blacklist=("" if config['singularity']['gridss'][genome_version]['blacklist'] == "" else " -b " + config['singularity']['gridss'][genome_version]['blacklist'])
#     singularity: config['singularity']['gridss']['sif']
#     shell:
#         """
#             java -cp esvee.jar com.hartwig.hmftools.esvee.assembly.AssemblyApplication 
#             -tumor TUMOR_SAMPLE_ID 
#             -reference REF_SAMPLE_ID
#             -tumor_bam /sample_data/TUMOR_SAMPLE_ID.bam
#             -reference_bam /sample_data/REF_SAMPLE_ID.bam
#             -junction_file /sample_data/output/TUMOR_SAMPLE_ID.esvee.prep.junction.tsv
#             -ref_genome /path_to_ref_genome_fasta/
#             -ref_genome_version 38
#             -write_types 'JUNC_ASSEMBLY;PHASED_ASSEMBLY;ALIGNMENTS;BREAKEND;VCF'
#             -output_dir /sample_data/output/ 
#             -threads 16
#         """


# rule esvee_call:
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#         ref=config['resources'][genome_version]['REFFA'],
#     output:
#         vcf="{project}/{genome_version}/results/sv/paired/esvee/{sample}/SV_gridss_{sample}.vcf",
#     params:
#         ref=config['resources'][genome_version]['REFFA'],
#         wd="{project}/{genome_version}/results/sv/paired/esvee/{sample}",
#         blacklist=get_gridss_blacklist
#         # if else like  (exprs ? c1 : c2) in C++/javascript is OK too,eg on next line
#         # blacklist=("" if config['singularity']['gridss'][genome_version]['blacklist'] == "" else " -b " + config['singularity']['gridss'][genome_version]['blacklist'])
#     singularity: config['singularity']['gridss']['sif']
#     shell:
#         """
#             java -cp esvee.jar com.hartwig.hmftools.esvee.depth.DepthAnnotator \
#             -sample 'TUMOR_SAMPLE_ID,REF_SAMPLE_ID'
#             -bam_file '/sample_data/TUMOR_SAMPLE_ID.bam,/sample_data/REF_SAMPLE_ID.bam'
#             -input_vcf TUMOR_SAMPLE_ID.esee.raw.vcf.gz
#             -output_vcf TUMOR_SAMPLE_ID.esee.ref_depth.vcf.gz
#             -ref_genome /path_to_ref_genome_fasta/
#             -ref_genome_version 38
#             -threads 16
#         """

# rule esvee_assemb:
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#         ref=config['resources'][genome_version]['REFFA'],
#     output:
#         vcf="{project}/{genome_version}/results/sv/paired/esvee/{sample}/SV_gridss_{sample}.vcf",
#     params:
#         ref=config['resources'][genome_version]['REFFA'],
#         wd="{project}/{genome_version}/results/sv/paired/esvee/{sample}",
#         blacklist=get_gridss_blacklist
#         # if else like  (exprs ? c1 : c2) in C++/javascript is OK too,eg on next line
#         # blacklist=("" if config['singularity']['gridss'][genome_version]['blacklist'] == "" else " -b " + config['singularity']['gridss'][genome_version]['blacklist'])
#     singularity: config['singularity']['gridss']['sif']
#     shell:
#         """
#             java -cp esvee.jar com.hartwig.hmftools.esvee.caller.CallerApplication 
#             -sample TUMOR_SAMPLE_ID
#             -reference REF_SAMPLE_ID
#             -input_vcf /sample_data/output/TUMOR_SAMPLE_ID.esvee.ref_depth.vcf.gz
#             -ref_genome_version 38
#             -pon_sgl_file /ref_data/sgl_pon.38.bed.gz
#             -pon_sv_file /ref_data/sv_pon.38.bedpe.gz
#             -known_hotspot_file /ref_data/known_fusions.38.bedpe
#             -repeat_mask_file /ref_data/repeat_mask_data.37.fa.gz
#             -output_dir /sample_data/output/ 
#         """





