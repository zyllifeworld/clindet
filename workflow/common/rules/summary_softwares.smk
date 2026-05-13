rule clindet_main_version:
    output:
        ""
    conda:
        
    shell:
        """
        conda list | grep -v ^# | awk -v var=clindet_main '{{ print var, $0 }}' 
        """

rule clindet_vep_version:
    output:
        ""
    conda:
        
    shell:
        """
        conda list | grep -v ^# | awk -v var=clindet_main '{{ print var, $0 }}' 
        """

rule clindet_hmftool_version:
    output:
        ""
    conda:
        
    shell:
        """
        conda list | grep -v ^# | awk -v var=clindet_main '{{ print var, $0 }}' 
        """

rule clindet_rsem_version:
    output:
        ""
    conda:
        
    shell:
        """
        conda list | grep -v ^# | awk -v var=clindet_main '{{ print var, $0 }}' 
        """


#### singularity software version