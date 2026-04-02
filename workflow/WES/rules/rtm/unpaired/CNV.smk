if genome_version in ['b37','hg38']:
    include:"CNV/purple.smk"
include: "CNV/freec.smk"