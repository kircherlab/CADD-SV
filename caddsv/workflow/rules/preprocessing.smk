rule prep_chr1:
    input:
        "beds/{set}/{set}{format}_CBinput.{bedflanks}",
    conda:
        "../envs/SV.yml"
    output:
        "beds/{set}/{set}{format}_wchr.{bedflanks}",
    shell:
        """
        cut -f1,2,3 {input} | sort -k1,1 -k2,2n > {output}
        """


rule prep_chr2:
    input:
        "beds/{set}/{set}{format}_wchr.{bedflanks}",
    conda:
        "../envs/SV.yml"
    output:
        "beds/{set}/{set}{format}_nochr.{bedflanks}",
    shell:
        """
        sed 's/^chr\\|%$//g' {input} > {output}
        """


rule prep_merg1:
    input:
        nochr="beds/{set}/{set}{format}_nochr.{bedflanks}",
        wchr="beds/{set}/{set}{format}_wchr.{bedflanks}",
    conda:
        "../envs/SV.yml"
    output:
        nochr="beds/{set}/{set}{format}_nochr_merged.{bedflanks}",
        wchr="beds/{set}/{set}{format}_merged.{bedflanks}",
    shell:
        """
        bedops --merge {input.nochr} > {output.nochr}
        bedops --merge {input.wchr} > {output.wchr}

        """
