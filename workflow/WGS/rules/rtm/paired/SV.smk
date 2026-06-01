#### delly workflow
include: "SV/delly.smk"
#### BRASS work flow
ascat_config = config['softwares_params'][genome_version].get('ascat_wgs', False)
if ascat_config and genome_version in ['b37','hg38']:
    include: "SV/BRASS.smk"
#### gridss workflow
include:"SV/gridss.smk"
#### svaba workflow
include:"SV/svaba.smk"
#### linx workflow
include:"SV/linx.smk"
#### igcaller workflow igcaller for B-cell
igcaller_config = config['singularity'].get('igcaller',{}).get('sif', False)
if igcaller_config:
    include:"SV/igcaller.smk"

include:"SV/jasmine_merge.smk"
