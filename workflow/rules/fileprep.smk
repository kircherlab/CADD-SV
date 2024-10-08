rule prep_files:
    input:
        "input/id_{set}.{format}",
    conda:
        "../envs/preprocessing.yml"
    output:
        CB="beds/{set}/{set}{format}_CBinput.bed",
    shell:
        """
        mkdir -p beds/{wildcards.set}
        bash workflow/scripts/input_{wildcards.format}_preprocessing.sh {input} | cut -f1,2,3,4 > {output.CB}
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
