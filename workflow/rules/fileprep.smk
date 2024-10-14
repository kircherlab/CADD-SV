rule prep_files:
    input:
        bed="input/id_{set}.{format}",
        genome="ucsc/hg38.fa.bgz"
    conda:
        "../envs/preprocessing.yml"
    output:
        CB="beds/{set}/{set}{format}_CBinput.bed",
        SB="beds/{set}/{set}{format}_SBinput.bed"
    shell:
        """
        mkdir -p beds/{wildcards.set}
        
        python workflow/scripts/retrieve_seq_for_nt.py <(bash workflow/scripts/input_{wildcards.format}_preprocessing.sh {input.bed}) {input.genome} > {output.SB}
        
        cut -f1,2,3,4 {output.SB} > {output.CB}
        """

rule flanks:
    input:
        CB="beds/{set}/{set}{format}_CBinput.bed",
        genome="ucsc/hg38.fa.sorted.genome",
    conda:
        "../envs/SV.yml"
    params:
        tmpup="beds/{set}/{set}{format}_CBinput.bedup_tmp",
        flanksize=config["flank"],
    output:
        up="beds/{set}/{set}{format}_CBinput.bedup",
        down="beds/{set}/{set}{format}_CBinput.beddown",
        flank="beds/{set}/{set}{format}_CBinput.bedflank"
    shell:
        """
        cat {input.CB} | awk 'BEGIN{{OFS = "\t"}}{{if ($2 == 0) $2+=1 ; print $0}}' > {params.tmpup}
        bedtools flank -i {params.tmpup} -g {input.genome} -l {params.flanksize} -r 0| awk 'BEGIN{{OFS = "\t"}}{{if ($2 == 0) $2+=1 ; print}}' | sort -k1,1 -k2,2n - > {output.up}
        bedtools flank -i {input.CB} -g {input.genome} -r {params.flanksize} -l 0 | sort -k1,1 -k2,2n - > {output.down}
        bedtools slop -i {input.CB} -g {input.genome} -b {params.flanksize} | sort -k1,1 -k2,2n - > {output.flank}
        
        """

rule split_bed_SB:
    input:
        SB="beds/{set}/{set}{format}_SBinput.bed"
    conda:
        "../envs/preprocessing.yml"
    output:
        flag="beds/{set}/{set}{format}_SB.flag"
    params:
        split_dir="beds/{set}/split_for_SB"
    shell:
        """
        mkdir -p {params.split_dir}
        split -l 5 {input.SB} {params.split_dir}/part_
        touch {output.flag}
        """

rule run_nt_script:
    input:
        "beds/{set}/{set}{format}_SB.flag"
    conda:
        "../envs/NT.yml"
    params:
        split_dir="beds/{set}/split_for_SB"
    output:
        flag="beds/{set}/{set}{format}_NT.flag"
    shell:
        """
        python workflow/scripts/run_segmentNT.py {params.split_dir}
        cat {input} > {output.flag}
        """

rule SB_features_generation:
    input:
        "beds/{set}/{set}{format}_NT.flag"
    conda:
        "../envs/NT.yml"
    params:
        split_dir="beds/{set}/split_for_SB"
    output:
         "{set}/{set}{format}_SB_features.bed"
    shell:
        """
        python workflow/scripts/SB_features_df_generator.py {params.split_dir} > {output}
        """
