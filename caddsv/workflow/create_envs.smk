# Minimal input-free workflow used by `caddsv get envs --use-conda`.
# Snakemake creates these environments and exits before running the shell jobs.

rule all_envs:
    input:
        ".caddsv-envs/preprocessing.ready",
        ".caddsv-envs/sv.ready",
        ".caddsv-envs/nt.ready",
        ".caddsv-envs/training.ready"


rule coordinate_based_envs:
    input:
        ".caddsv-envs/preprocessing.ready",
        ".caddsv-envs/sv.ready",
        ".caddsv-envs/training.ready"


rule preprocessing_env:
    output:
        ".caddsv-envs/preprocessing.ready"
    conda:
        "envs/preprocessing.yml"
    shell:
        "touch {output}"


rule sv_env:
    output:
        ".caddsv-envs/sv.ready"
    conda:
        "envs/SV.yml"
    shell:
        "touch {output}"


rule nt_env:
    output:
        ".caddsv-envs/nt.ready"
    conda:
        "envs/NT.yml"
    shell:
        "touch {output}"


rule training_env:
    output:
        ".caddsv-envs/training.ready"
    conda:
        "envs/training.yml"
    shell:
        "touch {output}"
