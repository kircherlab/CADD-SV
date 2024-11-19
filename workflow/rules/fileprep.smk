import os
gpu_available = bool(os.environ.get("CUDA_VISIBLE_DEVICES"))
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
        up="beds/{set}/{set}{format}_CBinput.bedup{flanksize}",
        down="beds/{set}/{set}{format}_CBinput.beddown{flanksize}",
    shell:
        """
        cat {input.CB} | awk 'BEGIN{{OFS = "\t"}}{{if ($2 == 0) $2+=1 ; print $0}}' > {params.tmpup}
        bedtools flank -i {params.tmpup} -g {input.genome} -l {params.flanksize} -r 0| awk 'BEGIN{{OFS = "\t"}}{{if ($2 == 0) $2+=1 ; print}}' | sort -k1,1 -k2,2n - > {output.up}
        bedtools flank -i {input.CB} -g {input.genome} -r {params.flanksize} -l 0 | sort -k1,1 -k2,2n - > {output.down}
        
        """


rule run_nt_script:
    input:
        "beds/{set}/{set}{format}_SBinput.bed"
    conda:
        "../envs/NT.yml"
    resources:
        gpu=1 if gpu_available else 0,
        mem_gb=200
    output:
        "beds/{set}/{set}{format}_SBprobabilities.h5"
    shell:
        """
        python workflow/scripts/run_segmentNT.py {input} {output}
        """

rule SB_features_generation:
    input:
        coordinates="beds/{set}/{set}{format}_SBinput.bed",
        probabilities="beds/{set}/{set}{format}_SBprobabilities.h5"
    conda:
        "../envs/NT.yml"
    output:
         "{set}/{set}{format}_SBfeatures.bed"
    shell:
        """
        python workflow/scripts/SB_features_df_generator.py {input.coordinates} {input.probabilities} > {output}
        """
