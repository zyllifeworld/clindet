rule make_region_bed_list:
    output:
        bed_list = "{project}/{genome_version}/results/bed_list.tsv"
    params:
        bed_dict = lambda wildcards: config["softwares_params"][wildcards.genome_version]['annotate_beds']
    run:
        import os
        import pandas as pd

        os.makedirs(os.path.dirname(output.bed_list), exist_ok=True)

        bed_dict = params.bed_dict

        df = pd.DataFrame({
            "name": list(bed_dict.keys()),
            "path": [resolve_bed_path(p) for p in bed_dict.values()]
        })

        df.to_csv(output.bed_list, sep="\t", index=False)

rule flag_mutation_pairead_maf:
    input: 
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf",
        bed_list="{project}/{genome_version}/results/bed_list.tsv"
    output: 
        maf="{project}/{genome_version}/results/maf/paired/{sample}/anno/{sample}.maf"
    conda:
        flexible_conda_env(config,['conda','clindet_mut'],env_yaml = 'envs/clindet_mut.yaml')
    log:
        "{project}/{genome_version}/results/maf/paired/{sample}/anno/{sample}.log"
    params:
        genome_version=lambda wildcards: wildcards.genome_version,
        chr_col="Chromosome",
        start_col="Start_Position",
        end_col="End_Position",
        keep_chr_prefix=False,
        ignore_strand=True,
        max_region_labels=10
    script:
        "../scripts/flag_mutation_maf.R"