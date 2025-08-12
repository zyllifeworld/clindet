include: "SNV/vardict.smk"
include: "SNV/HaplotypeCaller.smk"
if recall_pon or not pre_pon_db:
    include: "SNV/Mutect2_pon.smk"
else:
    include: "SNV/Mutect2.smk"

include: "SNV/Strelka.smk"
include: "SNV/Muse.smk"
include: "SNV/sentieon.smk"
if genome_version in ['b37','hg38']:
    include: "SNV/sage.smk"

include: "SNV/DeepVariant.smk"
