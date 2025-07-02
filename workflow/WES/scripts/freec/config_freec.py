import configparser
from snakemake.script import snakemake
config = configparser.ConfigParser()
config.optionxform = str ## don not covert to lower case!!!


## config temple
config_tem = snakemake.input['ini_template']
## input BAM file,for WES, a target panel bed
tumor_bam = snakemake.input['Tum']
normal_bam = snakemake.input['NC']
## freec params
chrLenFile = snakemake.params['chrLenFile']
outputDir = snakemake.params['outputDir']
chrFiles = snakemake.params['chrFiles']
bed = snakemake.params['bed']
### output config for each sample
out_config = snakemake.output['config']

config.read(config_tem)

config['general']['chrLenFile'] = chrLenFile
config['general']['outputDir'] = outputDir
config['general']['chrFiles'] = chrFiles

config['sample']['mateFile'] = tumor_bam
config['control']['mateFile'] = normal_bam
# config['BAF']['SNPfile'] = 
config['target']['captureRegions'] = bed

with open(out_config, 'w') as configfile:
    config.write(configfile)
