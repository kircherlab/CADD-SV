rule prep_chr1:
    input:
        "beds/{set}/{set}{format}_CBinput.{bedflanks}",
    conda:
        "../envs/SV.yml"
    container:
        container_for("sv")
    output:
        "beds/{set}/{set}{format}_wchr.{bedflanks}",
    shell:
        """
        # Annotation inputs must retain CB input order.
        cut -f1,2,3 {input} > {output}
        """


rule prep_chr2:
    input:
        "beds/{set}/{set}{format}_wchr.{bedflanks}",
    conda:
        "../envs/SV.yml"
    container:
        container_for("sv")
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
    container:
        container_for("sv")
    output:
        nochr="beds/{set}/{set}{format}_nochr_merged.{bedflanks}",
        wchr="beds/{set}/{set}{format}_merged.{bedflanks}",
    shell:
        """
        sort-bed {input.nochr} | bedops --merge - > {output.nochr}
        sort-bed {input.wchr} | bedops --merge - > {output.wchr}

        """
