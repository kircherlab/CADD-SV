import os
gpu_available = bool(os.environ.get("CUDA_VISIBLE_DEVICES"))
rule prep_files:
    input:
        bed="input/id_{set}.{format}",
        genome=ANNOT_DIR + "/ucsc/hg38.fa.bgz"
    conda:
        "../envs/preprocessing.yml"
    output:
        CB="beds/{set}/{set}{format}_CBinput.bed",
        SB="beds/{set}/{set}{format}_SBinput.bed",
        SB2="beds/{set}/{set}{format}_SBrefinput.bed"
    shell:
        """
        mkdir -p beds/{wildcards.set}
        
        python {workflow.basedir}/scripts/retrieve_seq_for_nt.py <(bash {workflow.basedir}/scripts/input_{wildcards.format}_preprocessing.sh {input.bed}) {input.genome} > {output.SB}
        python {workflow.basedir}/scripts/retrieve_seq_for_ref.py <(bash {workflow.basedir}/scripts/input_{wildcards.format}_preprocessing.sh {input.bed}) {input.genome} > {output.SB2}
        cut -f1,2,3,4 {output.SB} > {output.CB}
        """

rule flanks:
    input:
        CB="beds/{set}/{set}{format}_CBinput.bed",
        genome=ANNOT_DIR + "/ucsc/hg38.fa.sorted.genome",
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
if config["sequence_model"]:
    rule run_nt_script_ref:
        input:
            "beds/{set}/{set}{format}_SBrefinput.bed"
        conda:
            "../envs/NT.yml"
        resources:
            gpu=1 if gpu_available else 0,
            mem_gb=200
        output:
            "beds/{set}/{set}{format}_SBrefprobabilities.h5"
        shell:
            """
            python {workflow.basedir}/scripts/run_segmentNT.py {input} {output}
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
            python {workflow.basedir}/scripts/run_segmentNT.py {input} {output}
            """


    rule SB_features_generation:
        input:
            coordinates="beds/{set}/{set}{format}_SBinput.bed",
            probabilities="beds/{set}/{set}{format}_SBprobabilities.h5",
            refprobabilities="beds/{set}/{set}{format}_SBrefprobabilities.h5"
        conda:
            "../envs/NT.yml"
        output:
            SB = "beds/{set}/{set}{format}_SBfeatures.bed",
            SBref = "beds/{set}/{set}{format}_SBreffeatures.bed"
        shell:
            """
            python {workflow.basedir}/scripts/SB_features_df_generator.py {input.coordinates} {input.probabilities} > {output.SB}
            python {workflow.basedir}/scripts/SB_features_df_generator.py {input.coordinates} {input.refprobabilities} > {output.SBref}
            """


    rule DB_features_generation:
        input:
            coordinates="beds/{set}/{set}{format}_SBinput.bed",
            probabilities="beds/{set}/{set}{format}_SBprobabilities.h5",
            coordinates_ref="beds/{set}/{set}{format}_SBrefinput.bed",
            probabilities_ref="beds/{set}/{set}{format}_SBrefprobabilities.h5",
            SB = "beds/{set}/{set}{format}_SBfeatures.bed",
            SBref = "beds/{set}/{set}{format}_SBreffeatures.bed"
        conda:
            "../envs/NT.yml"
        output:
            DB = "beds/{set}/{set}{format}_DBfeatures.bed"
        shell:
            """
            python {workflow.basedir}/scripts/DB_features_df_generator.py \
        {input.coordinates} {input.probabilities} \
        {input.coordinates_ref} {input.probabilities_ref} > {output.DB}


            """
