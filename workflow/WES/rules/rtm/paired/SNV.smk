## rtm
include: "SNV/vardict.smk"
include: "SNV/HaplotypeCaller.smk"
if (recall_pon or not pre_pon_db) and (not custome_pon_db):
    include: "SNV/Mutect2_pon.smk"
elif custome_pon_db:
    include: "SNV/Mutect2_with_pon.smk"
else:
    include: "SNV/Mutect2.smk"

include: "SNV/Strelka.smk"
include: "SNV/varscan2.smk"
include: "SNV/Muse.smk"
include: "SNV/sage.smk"
include: "SNV/Lofreq.smk"
include: "SNV/DeepVariant.smk"
