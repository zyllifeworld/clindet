include: "SNV/vardict.smk"
include: "SNV/HaplotypeCaller.smk"
include: "SNV/Mutect2.smk"
if recall_pon or not pre_pon_db:
    # print('recall')
    include: "SNV/Mutect2_pon.smk"
    # include: "SNV/Mutect2.smk"
else:
    # print('not recall')
    include: "SNV/Mutect2.smk"
# include: "SNV/Mutect2_filter_tmp.smk"
include: "SNV/Strelka.smk"
include: "SNV/Muse.smk"
include: "SNV/sentieon.smk"
if genome_version in ['b37','hg38']:
    # print('recall')
    include: "SNV/sage.smk"
    # include: "SNV/Mutect2.smk"

include: "SNV/DeepVariant.smk"
