ascat_config = config['softwares']['ascat_wgs'].get(genome_version, False)
if ascat_config:
    include:"CNV/ASCAT.smk"
### PoN of factesCH
include:"CNV/facets.smk"
### purple only work for b37 and hg38
if genome_version in ['b37','hg38']:
    include:"CNV/purple.smk"

battenberg_config = config['softwares']['ascat_wgs'].get(genome_version, False)
if battenberg_config:
    include:"CNV/Battenberg.smk"
##### sequenza section
include:"CNV/sequenza.smk"
