def get_config_value(config_dict, keys, default=None, params=''):
    """
    安全获取嵌套配置值
    Args:
        config_dict: 配置字典
        keys: 键的列表，如 ['mutect2', 'intervals']
        default: 默认值
        param: 参数前缀，如 '-L' 或 '--germline-resource'
    Returns:
        如果 param 不为空: 返回 "param_prefix value"
        否则: 返回 value
    """
    current = config_dict
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            return default

    if current is None or current == '': #参数为空
        return default
    else:
        return f"{params} {current}"

def get_tumor_fastq(wildcards):
    fq_files = list(samples_info.loc[wildcards.sample,['Tumor_R1_file_path','Tumor_R2_file_path']])
    return {
        "R1":fq_files[0],
        "R2":fq_files[1],
    }

def get_normal_fastq(wildcards):
    fq_files = list(samples_info.loc[wildcards.sample,['Normal_R1_file_path','Normal_R2_file_path']])
    return {
        "R1":fq_files[0],
        "R2":fq_files[1],
    }

def get_sample_bed(wildcards):
    bed_file = samples_info.loc[wildcards.sample,'Target_file_bed']
    return(bed_file)

def get_gender(wildcards):
    gender = samples_info.loc[wildcards.sample,'gender']
    return(gender)


def get_vcf_name(wildcards):
    if wildcards.sample in paired_samples:
        if wildcards.caller == 'vardict' or wildcards.caller == 'vardict_filter' or wildcards.caller == 'vardict_germline':          
            name =  "--tumor-id " + wildcards.sample + '_T' + " --normal-id " + wildcards.sample + '_NC'
        elif wildcards.caller == 'sage' or wildcards.caller == 'pave':
            name =  "--tumor-id " + wildcards.sample + " --normal-id " + wildcards.sample + '_NC'
        elif wildcards.caller == 'cgppindel':
            name = "--tumor-id TUMOUR --normal-id NORMAL"
        elif wildcards.caller == 'caveman':
            name = "--tumor-id TUMOUR --normal-id NORMAL"
        elif wildcards.caller == 'muse':
            name = "--tumor-id TUMOR --normal-id NORMAL"
        elif wildcards.caller == 'varscan2':
            name = "--tumor-id TUMOR --normal-id NORMAL"
        else: 
            name =  "--tumor-id " + wildcards.sample + '_T' + " --normal-id " + wildcards.sample + '_NC'
    else:
        name =  "--tumor-id " + wildcards.sample + '_T' + " --normal-id " + wildcards.sample + '_NORMAL'
    return(name)

def vcf2vcf_name(wildcards):
    if wildcards.sample in paired_samples:
        if wildcards.caller == 'vardict' or wildcards.caller == 'vardict_filter' or wildcards.caller == 'vardict_germline':          
            name =  "--vcf-tumor-id " + wildcards.sample + '_T' + " --vcf-normal-id " + wildcards.sample + '_NC'
        elif wildcards.caller == 'sage' or wildcards.caller == 'pave':
            name =  "--vcf-tumor-id " + wildcards.sample + " --vcf-normal-id " + wildcards.sample + '_NC'
        elif wildcards.caller == 'cgppindel':
            name = "--vcf-tumor-id TUMOUR --vcf-normal-id NORMAL"
        elif wildcards.caller == 'caveman':
            name = "--vcf-tumor-id TUMOUR --vcf-normal-id NORMAL"
        elif wildcards.caller == 'muse':
            name = "--vcf-tumor-id TUMOR --vcf-normal-id NORMAL"
        elif wildcards.caller == 'varscan2':
            name = "--vcf-tumor-id TUMOR --vcf-normal-id NORMAL"
        else: 
            name =  "--vcf-tumor-id " + wildcards.sample + '_T' + " --vcf-normal-id " + wildcards.sample + '_NC'
    else:
        name =  "--vcf-tumor-id " + wildcards.sample + '_T' + " --vcf-normal-id " + wildcards.sample + '_NORMAL'
    return(name)

def get_vcf_file(wildcards):
    if wildcards.sample in paired_samples:
        file_name =  f"{project}/{genome_version}/results/vcf/" + 'paired' + "/{sample}/{caller}.vcf"
    else:
        file_name =  f"{project}/{genome_version}/results/vcf/" + 'unpaired' + "/{sample}/{caller}.vcf"
    return(file_name)

def choose_vcf2maf(wildcards):
    if vcf2maf == 'raw':
        vcf = '{project}/{genome_version}/results/vcf/paired/{sample}/{caller}.vcf'
    else:
        vcf = '{project}/{genome_version}/results/vcf_norm/paired/{sample}/{caller}.vcf'

def get_bcftools_filter_rules(wildcards):
    caller = wildcards.caller.lower()

    rule_map = {
        "mutect2": {
            "include": 'FILTER="PASS"',
            "exclude": ''
        },
        "vardict": {
            "include": 'FILTER="PASS" && QUAL >= 30',
            "exclude": ''
        },
        "strelka2": {
            "include": 'FILTER="PASS"',
            "exclude": ''
        },
        "caveman": {
            "include": 'FILTER="PASS"',
            "exclude": ''
        },
        "vardict": {
            "include": 'FILTER="PASS"',
            "exclude": ''
        },
        "haplotypecaller": {
            "include": 'FILTER="PASS" && QUAL >= 30',
            "exclude": 'MQ < 40'
        },
        "muse": {
            "include": '',
            "exclude": 'FILTER~"Tier5"'
        },
        "sage": {
            "include": '',
            "exclude": 'FILTER~"Tier5"'
        },
        "cgppindel": {
            "include": '',
            "exclude": 'FILTER~"FF010|FF012"'
        }
    }

    return rule_map.get(caller, {
        "include": 'FILTER="PASS"',
        "exclude": ''
    })


def get_bcftools_include_expr(wildcards):
    return get_bcftools_filter_rules(wildcards).get("include", "")


def get_bcftools_exclude_expr(wildcards):
    return get_bcftools_filter_rules(wildcards).get("exclude", "")


def build_bcftools_filter_cmd(wildcards):
    include_expr = get_bcftools_include_expr(wildcards)
    exclude_expr = get_bcftools_exclude_expr(wildcards)

    cmd_parts = []
    if include_expr:
        cmd_parts.append(f"-i '{include_expr}'")
    if exclude_expr:
        cmd_parts.append(f"-e '{exclude_expr}'")

    return " ".join(cmd_parts)

def flexible_conda_env(config_dict, keys, default=None):
    current = config_dict
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            return default

    if current is None or current == '': #参数为空
        return default
    else:
        return current
