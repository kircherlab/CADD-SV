# CTCF Peaks from wang_et_al
rule CTCF:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/CTCF/Wangetal_hg38.bed", # should I move it?
    conda:
        "../envs/SV.yml"
    output:
        temp("{set}/{set}{format}_ctcf.{bedflanks}"),
    shell:
        """
        bedtools coverage -b {input.anno} -a {input.bed} > {output}
        """

rule ultraconserved:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/ultraconserved/ultraconserved_hg38_muliz120M_sort.bed",
    conda:
        "../envs/SV.yml"
    output:
        temp("{set}/{set}{format}_ultraconserved.{bedflanks}"),
    shell:
        """
        bedtools coverage -b {input.anno} -a {input.bed} > {output}
        """

# gc content from ucsc
rule GC:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/GC/gc5Base.bedGraph.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        #t1=temp("{set}/{set}{format}_gc.{bedflanks}"),
        t1="{set}/{set}{format}_gc.{bedflanks}",
    shell:
        """
        python workflow/scripts/mean_map.py {input.bed} {input.anno}  > {output.t1}
        """

# ensemble genome build temporary genenames and attribute extraction
rule gene_model_tmp:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        merg="beds/{set}/{set}{format}_nochr_merged.{bedflanks}",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz",
    conda:
        "../envs/SV.yml"
    output:
        #temp("{set}/{set}{format}_gm_tmp.{bedflanks}"),
        "{set}/{set}{format}_gm_tmp.{bedflanks}",
    shell:
        """
        tabix {input.anno} -R {input.merg} | awk '{{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \\"\\";"; }}' | gtf2bed - | cut -f1,2,3,8,10  > {output}
        """

# genemodels
rule gene_model:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="{set}/{set}{format}_gm_tmp.{bedflanks}",
    conda:
        "../envs/SV.yml"
    output:
        #temp("{set}/{set}{format}_genemodel.{bedflanks}"),
        "{set}/{set}{format}_genemodel.{bedflanks}",
    shell:
        """
        paste <( bedtools coverage -b <(grep "exon" {input.anno}) -a {input.bed} | cut -f 1,2,3,5 ) <(bedtools coverage -b <(grep "transcript" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "gene" {input.anno} ) -a {input.bed} | cut -f 5) <( bedtools coverage -b <( grep "start_codon" {input.anno} ) -a {input.bed} | cut -f 5 ) <( bedtools coverage -b <(grep "stop_codon" {input.anno}) -a {input.bed} | cut -f 5) <(bedtools coverage -b <(grep "three_prime_utr" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "five_prime_utr" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "CDS" {input.anno}) -a {input.bed} | cut -f 5 ) > {output}
        """

rule gene_model_dist:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed",
    conda:
        "../envs/SV.yml"
    output:
        #temp("{set}/{set}{format}_genetic_dist.{bedflanks}"),
        "{set}/{set}{format}_genetic_dist.{bedflanks}",
    shell:
        """
        grep "exon" {input.anno} | bedtools closest -d -t first -b stdin -a {input.bed} |cut -f 1,2,3,9 |paste - <(grep "gene" {input.anno} | bedtools closest -d -t first -b stdin -a {input.bed} | cut -f 9)|paste - <(grep "start_codon" {input.anno} |bedtools closest -d -t first -b stdin -a {input.bed} | cut -f 9 ) >> {output} 
        
        """

# gene names necessary for PlI extraction
rule gene_names:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.GENES.bed",
    conda:
        "../envs/SV.yml"
    output:
        "{set}/{set}{format}_genenames.{bedflanks}",
    shell:
        """
        bedtools map -b {input.anno} -a {input.bed} -c 4 -o distinct > {output}
        """

# PLi
rule pli:
    input:
        gn="{set}/{set}{format}_genenames.{bedflanks}",
        pli="annotations/gnomad/pli_exac.csv",
    conda:
        "../envs/SV.yml"
    output:
        #temp("{set}/{set}{format}_pli.{bedflanks}"),
        "{set}/{set}{format}_pli.{bedflanks}",
    shell:
        """
        python workflow/scripts/PLIextract.py {input.gn} {input.pli} {output}
        
        """

# retrieving mean and max CADD scores from cadd_10score.bed.gz
rule cadd_PC_phylop:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/PhastCons/CADD_PC_PhyloP_scores.bed.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        #temp("{set}/{set}{format}_CADD_PC_PhyloP_maxsum.{bedflanks}"),
        "{set}/{set}{format}_CADD_PC_PhyloP_maxsum.{bedflanks}",
    shell:
        """
        python workflow/scripts/rule_cadd_PC_phylop.py {input.bed} {input.anno}  > {output}
        """

rule cadd2:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/CADD/CADD_GRCh38-v1.5.bedGraph_90q_12.bed.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        #temp("{set}/{set}{format}_cadd2_count.{bedflanks}"),
        "{set}/{set}{format}_cadd2_count.{bedflanks}",
    shell:
        """
        python workflow/scripts/count_map.py {input.bed} {input.anno}  > {output}
        """

rule gerp:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        #temp("{set}/{set}{format}_gerp_max.{bedflanks}"),
        "{set}/{set}{format}_gerp_max.{bedflanks}",
    shell:
        """
        python workflow/scripts/max_map.py {input.bed} {input.anno}  > {output}
        """

rule gerp2:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        #temp("{set}/{set}{format}_gerp2_count.{bedflanks}"),
        "{set}/{set}{format}_gerp2_count.{bedflanks}",
    shell:
        """
        python workflow/scripts/count_map.py {input.bed} {input.anno}  > {output}
        """


rule LINSIGHT:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/linsight/LINSIGHT_hg38_sort.bed.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        #temp("{set}/{set}{format}_linsight_sum.{bedflanks}"),
        "{set}/{set}{format}_linsight_sum.{bedflanks}",
    shell:
        """
         python workflow/scripts/sum_map.py {input.bed} {input.anno}  > {output}
        """

# TODO: cambiare Rscript in un python script (l'obiettivo è eliminare tutti i pacchetti R dall'env SV...)
# enhancer promotor links
rule EP:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        ep="annotations/enhancer-promoter-links/sorted_encode.bed",
    conda:
        "../envs/SV.yml"
    output:
        o1=temp("{set}/{set}{format}_EP.{bedflanks}"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.ep} | python workflow/scripts/annotateHIC.py /dev/stdin {output.o1}
        """

# frequently interacting regulatory elements
rule fire:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="annotations/FIRE/fire_{celllines}.bed",
    conda:
        "../envs/SV.yml"
    output:
        temp("{set}/{set}{format}_fire_{celllines}.{bedflanks}"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4,4 -o max,min > {output}
        """


rule fire2:
    input:
        expand("{{set}}/{{set}}{{format}}_fire_{fire}.{{bedflanks}}", fire=CL),
    output:
        temp("{set}/{set}{format}_fire.{bedflanks}"),
    shell:
        "paste {input} > {output}"


# hic data hESC hg18
rule HIC_hESC:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        hic="annotations/hic/hESC/combined/sorted.total.combined.domain",
    conda:
        "../envs/SV.yml"
    output:
        o1=temp("{set}/{set}{format}_HIC_hESC.{bedflanks}"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | python workflow/scripts/annotateHIC.py /dev/stdin {output.o1}
        """

# hic data from encode
rule HIC_encode:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        hic="annotations/Encode-HIC/{cells}/sorted_{tad}.bed.gz",
    conda:
        "../envs/SV.yml"
    output:
        o1=temp("{set}/{set}{format}_encode_{cells}_{tad}_hic.{bedflanks}"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | python workflow/scripts/annotateHIC.py /dev/stdin {output.o1}
        """


# merging encode HIC data
rule mergetad:
    input:
        expand("{{set}}/{{set}}{{format}}_encode_{cells}_{tad}_hic.{{bedflanks}}", cells=CELLS, tad=TAD),
    output:
        temp("{set}/{set}{format}_HIC.{bedflanks}"),
    shell:
        "paste {input} > {output}"

# micro syntenic regions
rule microsynteny:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        hic="annotations/synteny/microsynteny.bed",
    conda:
        "../envs/SV.yml"
    output:
        o1=temp("{set}/{set}{format}_microsynteny.{bedflanks}"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | python workflow/scripts/annotateHIC.py /dev/stdin {output.o1}
        """


rule ccr:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="annotations/ccr/ccrs.all.bed.gz",
    conda:
        "../envs/SV.yml"
    output:
        #temp("{set}/{set}{format}_ccr_mean.{bedflanks}"),
        "{set}/{set}{format}_ccr_mean.{bedflanks}",
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4 -o max > {output}
        """

# directionality index of HIC data from genomegitar database
rule genomegitar1:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/genomegitar/{gg}/DI_sort.bed",
    conda:
        "../envs/SV.yml"
    output:
        temp("{set}/{set}{format}_genomegitar_{gg}.{bedflanks}"),
    shell:
        """
        bedtools map -a {input.bed} \
        -b {input.anno} -c 4,4 -o max,min > {output}
        """

rule genomegitar2:
    input:
        expand("{{set}}/{{set}}{{format}}_genomegitar_{gg}.{{bedflanks}}", gg=genomegitars),
    output:
        temp("{set}/{set}{format}_DI.{bedflanks}"),
    shell:
        "paste {input} > {output}"

rule genomegitar3:
    input:
        "{set}/{set}{format}_DI.{bedflanks}",
    output:
        maxgg=temp("{set}/{set}{format}_DI_max.{bedflanks}"),
        mingg=temp("{set}/{set}{format}_DI_min.{bedflanks}"),
    shell:
        """
        cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i<m)m=$i;print m}}' > {output.mingg}
        cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i;print m}}' > {output.maxgg}"""


rule MPC:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="annotations/MPC/transcript_constraints_hg38liftover.bg.gz",
    conda:
        "../envs/SV.yml"
    output:
        temp("{set}/{set}{format}_MPC_mean.{bedflanks}"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4 -o mean > {output}
        """

rule RemapTF:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="annotations/ReMap/ReMap2_overlapTF_hg38.bg.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        #temp("{set}/{set}{format}_remapTF_mean.{bedflanks}"),
        "{set}/{set}{format}_remapTF_mean.{bedflanks}",
    shell:
        """
        python workflow/scripts/mean_map2.py {input.bed} {input.anno}  > {output}
        """

# various regulatory and epigenetic features from ENCODE (see ENCODE list)
rule encode:
    input:
        merg="beds/{set}/{set}{format}_nochr_merged.{bedflanks}",
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno= "annotations/encode/{encodes}/{encodes}_merged_90quant.bed.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        temp("{set}/{set}{format}_encode_{encodes}_mean.{bedflanks}"),
    shell:
        """
        python workflow/scripts/rule_encode.py {input.bed} {input.anno}  > {output}
        """

rule encode2:
    input:
        expand("{{set}}/{{set}}{{format}}_encode_{encodes}_mean.{{bedflanks}}", encodes=ENCODES),
    output:
        temp("{set}/{set}{format}_encode.{bedflanks}"),
    shell:
        "paste {input} > {output}"

# genome states infered from chromHMM
rule chromHMM_MAX:
    input:
        bed="beds/{set}/{set}{format}_nochr.{bedflanks}",
        anno="annotations/chromhmm/chromHMM_GRCh38.bg.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        #temp("{set}/{set}{format}_chromHMM_max.{bedflanks}"),
        "{set}/{set}{format}_chromHMM_max.{bedflanks}",
    shell:
        """
        python workflow/scripts/rule_chromHMM.py {input.bed} {input.anno}  > {output}
        """


rule Fantom5_counts:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/fantom5/F5.hg38.enhancers_sort.bed",
    conda:
        "../envs/SV.yml"
    output:
        temp("{set}/{set}{format}_f5_counts.{bedflanks}"),
    shell:
        """
        bedtools coverage -a {input.bed} -b {input.anno} -counts  > {output}
        """


rule HI:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/DDD_HI/hg38_HI_Predictions_version3_sort.bed",
    conda:
        "../envs/SV.yml"
    output:
        temp("{set}/{set}{format}_dddhi.{bedflanks}"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 5 -o max > {output}
        """

rule deepc:
    input:
        bed="beds/{set}/{set}{format}_wchr.{bedflanks}",
        anno="annotations/deepc/saliencies_merged_gm12878_5kb_10bp.bed.gz",
    conda:
        "../envs/preprocessing.yml"
    output:
        #temp("{set}/{set}{format}_deepc.{bedflanks}"),
        "{set}/{set}{format}_deepc.{bedflanks}",
    shell:
        """
        python workflow/scripts/max_map.py {input.bed} {input.anno}  > {output}
        """


rule complete_script:
    input:
        cadd="{set}/{set}{format}_CADD_PC_PhyloP_maxsum.{bedflanks}",
        cadd2="{set}/{set}{format}_cadd2_count.{bedflanks}",
        ccr="{set}/{set}{format}_ccr_mean.{bedflanks}",
        chromHMM="{set}/{set}{format}_chromHMM_max.{bedflanks}",
        ctcf="{set}/{set}{format}_ctcf.{bedflanks}",
        di_min="{set}/{set}{format}_DI_min.{bedflanks}",
        di_max="{set}/{set}{format}_DI_max.{bedflanks}",
        encode="{set}/{set}{format}_encode.{bedflanks}",
        ep="{set}/{set}{format}_EP.{bedflanks}",
        fire="{set}/{set}{format}_fire.{bedflanks}",
        gc="{set}/{set}{format}_gc.{bedflanks}",
        gm="{set}/{set}{format}_genemodel.{bedflanks}",
        gerp="{set}/{set}{format}_gerp_max.{bedflanks}",
        gerp2="{set}/{set}{format}_gerp2_count.{bedflanks}",
        hic="{set}/{set}{format}_HIC.{bedflanks}",
        hesc="{set}/{set}{format}_HIC_hESC.{bedflanks}",
        microsyn="{set}/{set}{format}_microsynteny.{bedflanks}",
        mpc="{set}/{set}{format}_MPC_mean.{bedflanks}",
        pli="{set}/{set}{format}_pli.{bedflanks}",
        remapTF="{set}/{set}{format}_remapTF_mean.{bedflanks}",
        f5="{set}/{set}{format}_f5_counts.{bedflanks}",
        hi="{set}/{set}{format}_dddhi.{bedflanks}",
        deepc="{set}/{set}{format}_deepc.{bedflanks}",
        ultrac="{set}/{set}{format}_ultraconserved.{bedflanks}",
        g_dist="{set}/{set}{format}_genetic_dist.{bedflanks}",
        linsight="{set}/{set}{format}_linsight_sum.{bedflanks}",
        header="annotations/header.txt",
    output:
        "{set}/{set}{format}_matrix.{bedflanks}",
    conda:
        "../envs/SV.yml"
    shell:
        """
        paste <(cut -f1-11 {input.cadd}) <(cut -f4 {input.cadd2}) <(cut -f4 {input.ccr}) <(cut -f4-28 {input.chromHMM}) <(cut -f4,5,7 {input.ctcf}) <(cut -f1 {input.di_min}) <(cut -f1 {input.di_max}) <(cut -f4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65 {input.encode}) <(cut -f4-7 {input.ep}) <(cut -f4,5,9,10,14,15,19,20,24,25 {input.fire}) <(cut -f4 {input.gc}) <(cut -f4-11 {input.gm}) <(cut -f4 {input.gerp}) <(cut -f4 {input.gerp2}) <(cut -f4,5,6,7,11,12,13,14,18,19,20,21,25,26,27,28 {input.hic}) <(cut -f4-7 {input.hesc}) <(cut -f4-7 {input.microsyn}) <(cut -f4 {input.mpc}) <(cut -f4 {input.pli}) <(cut -f4,5,6 {input.g_dist}) <(cut -f4 {input.remapTF}) <(cut -f4 {input.f5}) <(cut -f4 {input.hi}) <(cut -f4 {input.deepc}) <(cut -f4,5,7 {input.ultrac}) <(cut -f4 {input.linsight})| cat {input.header} - > {output}
        """
