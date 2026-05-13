### caveman_normal_pane
rule battenberg_call: 
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam"
    output:
        # out_dir=directory('results/vcf/paired/{sample}/caveman'),
        log='logs/paired/cgpbattenberg_{sample}.log'
    threads: 20
    params:
        ref=config['resources'][genome_version]['REFFA'],
        igbed=config['softwares_params'][genome_version]['cgpbattenberg']['ignorebed'],
        gc=config['softwares_params'][genome_version]['cgpbattenberg']['gc'],
        impute=config['softwares_params'][genome_version]['cgpbattenberg']['impute'],
        prob_loci=config['softwares_params'][genome_version]['cgpbattenberg']['prob_loci'],
        loci_1000=config['softwares_params'][genome_version]['cgpbattenberg']['loci_1000'],
        out_dir='results/vcf/paired/{sample}/battenberg',
        log='results/vcf/paired/{sample}/battenberg/log'
        # gender=get_gender
    singularity:
        flexible_container_img(config,['singularity','cgpbattenberg','sif'],image_url = config['softwares_params'][genome_version]['cgpbattenberg']['repo'])
    shell:
        """
        battenberg.pl -outdir {params.out_dir} \
        -r {params.ref}.fai \
        -tb {input.Tum} -prob-loci {params.prob_loci} \
        -nb {input.NC} \
        -ge XY \
        -impute-info {params.impute} \
        -thousand-genomes-loc  {params.loci_1000}\
        -ignore-contigs-file {params.igbed} \
        -gc-correction-loc {params.gc} \
        -species Human \
        -assembly 37 \
        -t {threads}
        touch {output.log}
        """