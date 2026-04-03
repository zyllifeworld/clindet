## 
rule maf2vcf:
    input:
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.maf",
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}.vcf"
    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        out_dir="{project}/{genome_version}/results/maf/paired/{sample}/merge"
    shell:
        """
        maf2vcf.pl --input-maf {input.maf} \
        --output-dir {params.out_dir} --ref-fasta {input.ref} \
        --output-vcf {output.vcf}  
        sed -i 's/ID=AD,Number=R/ID=AD,Number=2/g' {output.vcf} 
        """

rule pcgr_anno:
    input:
        vcf="{project}/{genome_version}/results/vcf2vcf/paired/{sample}/merge/{sample}.vcf",
        ref=config['resources'][genome_version]['REFFA']
    output:

    conda:
        config['softwares']['vcf2maf']['conda']
    params:
        out_dir="{project}/{genome_version}/results/maf/paired/{sample}/merge"
    shell:
        """
        maf2vcf.pl --input-maf {input.maf} \
        --output-dir {params.out_dir} --ref-fasta {input.ref} \
        --output-vcf {output.vcf}  
        sed -i 's/ID=AD,Number=R/ID=AD,Number=2/g' {output.vcf} 
        """