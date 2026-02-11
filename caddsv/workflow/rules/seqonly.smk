# Sequence-only mode rules for CADD-SV
# Processes input with 2 sequences (sample and reference) without coordinate-based annotations

import os

gpu_available = bool(os.environ.get("CUDA_VISIBLE_DEVICES"))


rule seqonly_preprocess:
    """
    Preprocess sequence-only input: shrink middle if > 1008bp, transform N runs.
    Input format: TSV with columns REF, ALT, [TYPE], [ID]
    """
    input:
        tsv="input/id_{set}.tsv"
    output:
        sample="beds/{set}/{set}_seqonly_sample.bed",
        ref="beds/{set}/{set}_seqonly_ref.bed"
    conda:
        "../envs/preprocessing.yml"
    shell:
        """
        mkdir -p beds/{wildcards.set}
        python {workflow.basedir}/scripts/preprocess_sequences.py {input.tsv} {output.sample} {output.ref}
        """


rule seqonly_run_nt_sample:
    """Run SegmentNT on sample sequences."""
    input:
        "beds/{set}/{set}_seqonly_sample.bed"
    conda:
        "../envs/NT.yml"
    resources:
        gpu=1 if gpu_available else 0,
        mem_gb=200
    output:
        "beds/{set}/{set}_seqonly_sample.h5"
    shell:
        """
        python {workflow.basedir}/scripts/run_segmentNT.py {input} {output}
        """


rule seqonly_run_nt_ref:
    """Run SegmentNT on reference sequences."""
    input:
        "beds/{set}/{set}_seqonly_ref.bed"
    conda:
        "../envs/NT.yml"
    resources:
        gpu=1 if gpu_available else 0,
        mem_gb=200
    output:
        "beds/{set}/{set}_seqonly_ref.h5"
    shell:
        """
        python {workflow.basedir}/scripts/run_segmentNT.py {input} {output}
        """


rule seqonly_SB_features:
    """Generate SB features from sample SegmentNT output."""
    input:
        coordinates="beds/{set}/{set}_seqonly_sample.bed",
        probabilities="beds/{set}/{set}_seqonly_sample.h5"
    conda:
        "../envs/NT.yml"
    output:
        "beds/{set}/{set}_seqonly_SB.bed"
    shell:
        """
        python {workflow.basedir}/scripts/SB_features_df_generator.py {input.coordinates} {input.probabilities} > {output}
        """


rule seqonly_SBref_features:
    """Generate SBref features from reference SegmentNT output."""
    input:
        coordinates="beds/{set}/{set}_seqonly_ref.bed",
        probabilities="beds/{set}/{set}_seqonly_ref.h5"
    conda:
        "../envs/NT.yml"
    output:
        "beds/{set}/{set}_seqonly_SBref.bed"
    shell:
        """
        python {workflow.basedir}/scripts/SB_features_df_generator.py {input.coordinates} {input.probabilities} > {output}
        """


rule seqonly_DB_features:
    """Generate DB features (delta between sample and reference)."""
    input:
        coordinates="beds/{set}/{set}_seqonly_sample.bed",
        probabilities="beds/{set}/{set}_seqonly_sample.h5",
        coordinates_ref="beds/{set}/{set}_seqonly_ref.bed",
        probabilities_ref="beds/{set}/{set}_seqonly_ref.h5",
        SB="beds/{set}/{set}_seqonly_SB.bed",
        SBref="beds/{set}/{set}_seqonly_SBref.bed"
    conda:
        "../envs/NT.yml"
    output:
        "beds/{set}/{set}_seqonly_DB.bed"
    shell:
        """
        python {workflow.basedir}/scripts/DB_features_df_generator.py \
            {input.coordinates} {input.probabilities} \
            {input.coordinates_ref} {input.probabilities_ref} > {output}
        """


rule seqonly_scoring:
    """Score variants using the sequence-only model."""
    input:
        SB="beds/{set}/{set}_seqonly_SB.bed",
        SBref="beds/{set}/{set}_seqonly_SBref.bed",
        DB="beds/{set}/{set}_seqonly_DB.bed"
    conda:
        "../envs/training.yml"
    output:
        "beds/{set}/output/{set}_seqonly_score.tsv"
    shell:
        """
        mkdir -p beds/{wildcards.set}/output
        python {workflow.basedir}/scripts/scoring_seqonly.py {input.SB} {input.SBref} {input.DB} {output}
        python {workflow.basedir}/scripts/add_phred.py {output}
        """
