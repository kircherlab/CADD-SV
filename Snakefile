# Snakefile
# Annotation script for feature retrieval
# CADD-SV

# this container defines the underlying OS for each job when using the workflow
# # with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


CELLS = ["A549", "Caki2"]
TAD = ["nested", "tad"]

# set of encode datasets in annotation/
ENCODES = [
    "DNase-seq",
    "H2AFZ",
    "H3K27ac",
    "H3K27me3",
    "H3K36me3",
    "H3K4me1",
    "H3K4me2",
    "H3K4me3",
    "H3K79me2",
    "H3K9ac",
    "H3K9me3",
    "H4K20me1",
    "totalRNA-seq",
]

# set of genomegitar datasets for HIC directionality index
genomegitars = [
    "GSM1055800_DI",
    "GSM1055805_DI",
    "GSM1081530_DI",
    "GSM1267196_DI",
    "GSM1267200_DI",
    "GSM1294038_DI",
    "GSM1294039_DI",
    "GSM1551599_DI",
    "GSM1551629_DI",
    "GSM1608505_DI",
    "GSM1718021_DI",
    "GSM1906332_DI",
    "GSM1906333_DI",
    "GSM1906334_DI",
    "GSM1909121_DI",
    "GSM455133_DI",
    "GSM862723_DI",
    "GSM862724_DI",
    "GSM927075_DI",
]

# Cell lines from Encode for HIC datasets
CL = ["gm12878", "msc", "mes", "imr90", "h1"]


rule all:
    input:
      matrix=expand("output/{sets}_score.bed", sets=config["dataset"])


rule prep_chr1:
    input:
        "input/id_{set}.bed",
    conda:
        "envs/SV.yml"
    output:
        "beds/{set}/{set}_wchr.bed",
    shell:
        """
        cut -f1,2,3 {input} > {output}
        """


rule prep_chr2:
    input:
        "beds/{set}/{set}_wchr.bed",
    conda:
        "envs/SV.yml"
    output:
        "beds/{set}/{set}_nochr.bed",
    shell:
        """
        sed 's/^chr\\|%$//g' {input} > {output}
        """


rule prep_merg1:
    input:
        nochr="beds/{set}/{set}_nochr.bed",
        wchr="beds/{set}/{set}_wchr.bed",
    conda:
        "envs/SV.yml"
    output:
        nochr="beds/{set}/{set}_nochr_merged.bed",
        wchr="beds/{set}/{set}_merged.bed",
    shell:
        """
        bedtools merge -i {input.nochr} > {output.nochr}
        bedtools merge -i {input.wchr} > {output.wchr}

        """


# CTCF Peaks from wang_et_al
rule CTCF:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/CTCF/Wangetal_hg38.bed",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/{set}_ctcf.bed"),
    shell:
        """
        bedtools coverage -b {input.anno} -a {input.bed} > {output.t1}
        """


rule ultraconserved:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/ultraconserved/ultraconserved_hg38_muliz120M_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/{set}_ultraconserved.bed"),
    shell:
        """
        bedtools coverage -b {input.anno} -a {input.bed} > {output.t1}
        """


# gc content from ucsc
rule GC:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/GC/gc5Base.bedGraph.gz",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/{set}_gc.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o mean; done < {input.bed}) > {output.t1}
        """


# ensemble genome build temporary genenames and attribute extraction
rule gene_model_tmp:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        merg="beds/{set}/{set}_nochr_merged.bed",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/gm_tmp.bed"),
    shell:
        """
        tabix {input.anno} -R {input.merg} | awk '{{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \\"\\";"; }}' | gtf2bed - | cut -f1,2,3,8,10  > {output.t1}
        """


# genemodels
rule gene_model:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="{set}/gm_tmp.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_genemodel.bed"),
    shell:
        """
        paste <( bedtools coverage -b <(grep "exon" {input.anno}) -a {input.bed} | cut -f 1,2,3,5 ) <(bedtools coverage -b <(grep "transcript" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "gene" {input.anno} ) -a {input.bed} | cut -f 5) <( bedtools coverage -b <( grep "start_codon" {input.anno} ) -a {input.bed} | cut -f 5 ) <( bedtools coverage -b <(grep "stop_codon" {input.anno}) -a {input.bed} | cut -f 5) <(bedtools coverage -b <(grep "three_prime_utr" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "five_prime_utr" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "CDS" {input.anno}) -a {input.bed} | cut -f 5 ) > {output}
        """


# grep "exon" {input.anno} | bedtools coverage -b stdin -a {input.bed} |cut -f 1,2,3,5 |paste - <(grep "transcript" {input.anno} |bedtools coverage -b stdin -a {input.bed} |cut -f 5 ) |paste - <(grep "gene" {input.anno} | bedtools coverage -b stdin -a {input.bed} | cut -f 5)|paste - <(grep "start_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "stop_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5) |paste - <(grep "three_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "five_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "CDS" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 ) >> {output}   """


rule gene_model_dist:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_genetic_dist.bed"),
    shell:
        """
        grep "exon" {input.anno} | bedtools closest -d -t first -b stdin -a {input.bed} |cut -f 1,2,3,9 |paste - <(grep "gene" {input.anno} | bedtools closest -d -t first -b stdin -a {input.bed} | cut -f 9)|paste - <(grep "start_codon" {input.anno} |bedtools closest -d -t first -b stdin -a {input.bed} | cut -f 9 ) >> {output}   """


# gene names necessary for PlI extraction
rule gene_names:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.GENES.bed",
    conda:
        "envs/SV.yml"
    output:
        "{set}/{set}_genenames.bed",
    shell:
        """
        bedtools map -b {input.anno} -a {input.bed} -c 4 -o distinct > {output}
        """


# PLi
rule pli:
    input:
        gn="{set}/{set}_genenames.bed",
        pli="annotations/gnomad/pli_exac.csv",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_pli.bed"),
    shell:
        """
        Rscript --vanilla scripts/PLIextract.R {input.gn} {input.pli} {output}
        """


# retrieving mean and max CADD scores from cadd_10score.bed.gz
rule cadd_PC_phylop:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/PhastCons/CADD_PC_PhyloP_scores.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_CADD_PC_PhyloP_maxsum.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}')| cat annotations/dummy8.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4,4,5,5,6,6,7,7,8,8 -o max,sum,max,sum,max,sum,max,sum,max,sum; done < {input.bed}) > {output}
        """


rule cadd2:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/CADD/CADD_GRCh38-v1.5.bedGraph_90q_12.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_cadd2_count.bed"),
    shell:
        """
        (while read -r line; do bedtools coverage -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}')| cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -counts; done < {input.bed}) > {output}
        """


# tabix {input.anno} -R {input.merg} | bedtools coverage -a {input.bed} -b stdin -counts > {output}
# merg="beds/{set}/{set}_merged.bed",


rule gerp:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_gerp_mean.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}')| cat annotations/dummy5_nochr.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o max; done < {input.bed}) > {output}
        """


rule gerp2:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_gerp2_count.bed"),
    shell:
        """
        (while read -r line; do bedtools coverage -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy5_nochr.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -counts; done < {input.bed}) > {output}
        """


# see cadd2
rule LINSIGHT:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/linsight/LINSIGHT_hg38_sort.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_linsight_sum.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o sum; done < {input.bed}) > {output}

        """


# tabix {input.anno} -R {input.merg} | awk '!seen[$0]++' | bedtools map -a {input.bed} -b stdin -c 4 -o sum > {output}


# enhancer promotor links
rule EP:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        ep="annotations/enhancer-promoter-links/sorted_encode.bed",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_EP.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.ep} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


# frequently interacting regulatory elements
rule fire:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/FIRE/fire_{celllines}.bed",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/{set}_fire_{celllines}.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4,4 -o max,min > {output.t1}
        """


rule fire2:
    input:
        expand("{{set}}/{{set}}_fire_{fire}.bed", fire=CL),
    output:
        temp("{set}/{set}_fire.bed"),
    shell:
        "paste {input} > {output}"


# hic data hESC hg18
rule HIC_hESC:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        hic="annotations/hic/hESC/combined/sorted.total.combined.domain",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_HIC_hESC.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


# hic data from encode
rule HIC_encode:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        hic="annotations/Encode-HIC/{cells}/sorted_{tad}.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_encode_{cells}_{tad}_hic.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


# merging encode HIC data
rule mergetad:
    input:
        expand("{{set}}/{{set}}_encode_{cells}_{tad}_hic.bed", cells=CELLS, tad=TAD),
    output:
        temp("{set}/{set}_HIC.bed"),
    shell:
        "paste {input} > {output}"


# micro syntenic regions
rule microsynteny:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        hic="annotations/synteny/microsynteny.bed",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_microsynteny.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


rule ccr:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/ccr/ccrs.all.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_ccr_mean.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4 -o max > {output}
        """


# directionality index of HIC data from genomegitar database
rule genomegitar1:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/genomegitar/{gg}/DI_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/{set}_genomegitar_{gg}.bed"),
    shell:
        """
        bedtools map -a {input.bed} \
        -b {input.anno} -c 4,4 -o max,min > {output.t1}
        """


rule genomegitar2:
    input:
        expand("{{set}}/{{set}}_genomegitar_{gg}.bed", gg=genomegitars),
    output:
        temp("{set}/{set}_DI.bed"),
    shell:
        "paste {input} > {output}"


rule genomegitar3:
    input:
        "{set}/{set}_DI.bed",
    output:
        maxgg=temp("{set}/{set}_DI_max.bed"),
        mingg=temp("{set}/{set}_DI_min.bed"),
    shell:
        """
        cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i<m)m=$i;print m}}' > {output.mingg}
        cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i;print m}}' > {output.maxgg}"""


rule MPC:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/MPC/transcript_constraints_hg38liftover.bg.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_MPC_mean.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4 -o mean > {output}
        """


rule RemapTF:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/ReMap/ReMap2_overlapTF_hg38.bg.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_remapTF_mean.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4_nochr.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o mean; done < {input.bed}) > {output}
        """


# various regulatory and epigenetic features from ENCODE (see ENCODE list)
rule encode:
    input:
        merg="beds/{set}/{set}_nochr_merged.bed",
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/encode/{encodes}/{encodes}_merged_90quant.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_encode_{encodes}_mean.bed"),
    shell:
        "tabix {input.anno} -R {input.merg} | cat annotations/dummy5_nochr.bed - | bedtools map -a {input.bed} -b - -c 4,4 -o max,sum > {output}"


rule encode2:
    input:
        expand("{{set}}/{{set}}_encode_{encodes}_mean.bed", encodes=ENCODES),
    output:
        temp("{set}/{set}_encode.bed"),
    shell:
        "paste {input} > {output}"


# genome states infered from chromHMM
rule chromHMM_MAX:
    input:
        bed="beds/{set}/{set}_nochr.bed",
        anno="annotations/chromhmm/chromHMM_GRCh38.bg.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_chromHMM_max.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(cat annotations/dummy_chrhmm.bed; tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 -o max; done < {input.bed}) > {output}
        """


#  bedtools map -a {input.bed} -b stdin -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 -o max > {output}


rule Fantom5_counts:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/fantom5/F5.hg38.enhancers_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_f5_counts.bed"),
    shell:
        """
        bedtools coverage -a {input.bed} -b {input.anno} -counts  > {output}
        """


rule HI:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/DDD_HI/hg38_HI_Predictions_version3_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_dddhi.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 5 -o max > {output}
        """


rule deepc:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        anno="annotations/deepc/saliencies_merged_gm12878_5kb_10bp.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_deepc.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o max; done < {input.bed}) > {output}
        """


rule complete_script:
    input:
        cadd="{set}/{set}_CADD_PC_PhyloP_maxsum.bed",
        cadd2="{set}/{set}_cadd2_count.bed",
        ccr="{set}/{set}_ccr_mean.bed",
        chromHMM="{set}/{set}_chromHMM_max.bed",
        ctcf="{set}/{set}_ctcf.bed",
        di_min="{set}/{set}_DI_min.bed",
        di_max="{set}/{set}_DI_max.bed",
        encode="{set}/{set}_encode.bed",
        ep="{set}/{set}_EP.bed",
        fire="{set}/{set}_fire.bed",
        gc="{set}/{set}_gc.bed",
        gm="{set}/{set}_genemodel.bed",
        gerp="{set}/{set}_gerp_mean.bed",
        gerp2="{set}/{set}_gerp2_count.bed",
        hic="{set}/{set}_HIC.bed",
        hesc="{set}/{set}_HIC_hESC.bed",
        microsyn="{set}/{set}_microsynteny.bed",
        mpc="{set}/{set}_MPC_mean.bed",
        pli="{set}/{set}_pli.bed",
        remapTF="{set}/{set}_remapTF_mean.bed",
        f5="{set}/{set}_f5_counts.bed",
        hi="{set}/{set}_dddhi.bed",
        deepc="{set}/{set}_deepc.bed",
        ultrac="{set}/{set}_ultraconserved.bed",
        g_dist="{set}/{set}_genetic_dist.bed",
        linsight="{set}/{set}_linsight_sum.bed",
        header="annotations/header.txt",
    output:
        "{set}/matrix.bed",
    conda:
        "envs/SV.yml"
    shell:
        """
        paste <(cut -f1-11 {input.cadd}) <(cut -f4 {input.cadd2}) <(cut -f4 {input.ccr}) <(cut -f4-28 {input.chromHMM}) <(cut -f4,5,7 {input.ctcf}) <(cut -f1 {input.di_min}) <(cut -f1 {input.di_max}) <(cut -f4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65 {input.encode}) <(cut -f4-7 {input.ep}) <(cut -f4,5,9,10,14,15,19,20,24,25 {input.fire}) <(cut -f4 {input.gc}) <(cut -f4-11 {input.gm}) <(cut -f4 {input.gerp}) <(cut -f4 {input.gerp2}) <(cut -f4,5,6,7,11,12,13,14,18,19,20,21,25,26,27,28 {input.hic}) <(cut -f4-7 {input.hesc}) <(cut -f4-7 {input.microsyn}) <(cut -f4 {input.mpc}) <(cut -f4 {input.pli}) <(cut -f4,5,6 {input.g_dist}) <(cut -f4 {input.remapTF}) <(cut -f4 {input.f5}) <(cut -f4 {input.hi}) <(cut -f4 {input.deepc}) <(cut -f4,5,7 {input.ultrac}) <(cut -f4 {input.linsight})| cat {input.header} - > {output}
        """


###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


rule prep_chr_100bpup:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        genome="annotations/hs38.fa.genome",
    output:
        wchr="beds/{set}/{set}_100bpup.bed",
        nchr="beds/{set}/{set}_100bpup_nchr.bed",
    conda:
        "envs/SV.yml"
    shell:
        """
        bedtools flank -i {input.bed} -g {input.genome} -l 100 -r 0 > {output.wchr}
        bedtools flank -i {input.bed} -g {input.genome} -l 100 -r 0 | sed 's/^chr\|%$//g' > {output.nchr}

        """


rule prep_merg1_100bpup:
    input:
        nchr="beds/{set}/{set}_100bpup_nchr.bed",
        wchr="beds/{set}/{set}_100bpup.bed",
    conda:
        "envs/SV.yml"
    output:
        nchr="beds/{set}/{set}_100bpup_nchr_merged.bed",
        wchr="beds/{set}/{set}_100bpup_merged.bed",
    shell:
        """
        bedtools merge -i {input.nchr} > {output.nchr}
        bedtools merge -i {input.wchr} > {output.wchr}

        """


# CTCF Peaks from wang_et_al
rule CTCF_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/CTCF/Wangetal_hg38.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_ctcf.bed"),
    shell:
        """
        bedtools coverage -b {input.anno} -a {input.bed} > {output}
        """


rule ultraconserved_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/ultraconserved/ultraconserved_hg38_muliz120M_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_ultraconserved.bed"),
    shell:
        """
        bedtools coverage -b {input.anno} -a {input.bed} > {output}
        """


# gc content from ucsc
rule GC_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/GC/gc5Base.bedGraph.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_gc.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o mean; done < {input.bed}) > {output}
        """


# ensemble genome build temporary genenames and attribute extraction
rule gene_model_tmp_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        merg="beds/{set}/{set}_100bpup_nchr_merged.bed",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/gm_tmp_100bpup.bed"),
    shell:
        """
        tabix {input.anno} -R {input.merg} | awk '{{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \\"\\";"; }}' | gtf2bed - | cut -f1,2,3,8,10  > {output.t1}
        """


# genemodels
rule gene_model_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="{set}/gm_tmp_100bpup.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_genemodel.bed"),
    shell:
        """
        paste <( bedtools coverage -b <(grep "exon" {input.anno}) -a {input.bed} | cut -f 1,2,3,5 ) <(bedtools coverage -b <(grep "transcript" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "gene" {input.anno} ) -a {input.bed} | cut -f 5) <( bedtools coverage -b <( grep "start_codon" {input.anno} ) -a {input.bed} | cut -f 5 ) <( bedtools coverage -b <(grep "stop_codon" {input.anno}) -a {input.bed} | cut -f 5) <(bedtools coverage -b <(grep "three_prime_utr" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "five_prime_utr" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "CDS" {input.anno}) -a {input.bed} | cut -f 5 ) > {output}
        """


# grep "exon" {input.anno} | bedtools coverage -b stdin -a {input.bed} |cut -f 1,2,3,5 |paste - <(grep "transcript" {input.anno} |bedtools coverage -b stdin -a {input.bed} |cut -f 5 ) |paste - <(grep "gene" {input.anno} | bedtools coverage -b stdin -a {input.bed} | cut -f 5)|paste - <(grep "start_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "stop_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5) |paste - <(grep "three_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "five_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "CDS" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 ) >> {output}   """


rule gene_model_dist_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_genetic_dist.bed"),
    shell:
        """
        grep "exon" {input.anno} | bedtools closest -d -t first -b stdin -a {input.bed} |cut -f 1,2,3,9 |paste - <(grep "gene" {input.anno} | bedtools closest -d -t first -b stdin -a {input.bed} | cut -f 9)|paste - <(grep "start_codon" {input.anno} |bedtools closest -d -t first -b stdin -a {input.bed} | cut -f 9 ) >> {output}   """


# gene names necessary for PlI extraction
rule gene_names_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.GENES.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_genenames.bed"),
    shell:
        """
        bedtools map -b {input.anno} -a {input.bed} -c 4 -o distinct > {output}
        """


# PLi
rule pli_100bpup:
    input:
        gn="{set}/{set}_100bpup_genenames.bed",
        pli="annotations/gnomad/pli_exac.csv",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_pli.bed"),
    shell:
        """
        Rscript --vanilla scripts/PLIextract.R {input.gn} {input.pli} {output}
        """


# retrieving mean and max CADD scores from cadd_10score.bed.gz
rule cadd_PC_phylop_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/PhastCons/CADD_PC_PhyloP_scores.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_CADD_PC_PhyloP_maxsum.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}')| cat annotations/dummy8.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4,4,5,5,6,6,7,7,8,8 -o max,sum,max,sum,max,sum,max,sum,max,sum; done < {input.bed}) > {output}
        """


rule cadd2_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/CADD/CADD_GRCh38-v1.5.bedGraph_90q_12.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_cadd2_count.bed"),
    shell:
        """
        (while read -r line; do bedtools coverage -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}')| cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -counts; done < {input.bed}) > {output}
        """


# tabix {input.anno} -R {input.merg} | bedtools coverage -a {input.bed} -b stdin -counts > {output}
# merg="beds/{set}/{set}_merged.bed",


rule gerp_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_gerp_mean.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}')| cat annotations/dummy5_nchr.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o max; done < {input.bed}) > {output}
        """


rule gerp2_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_gerp2_count.bed"),
    shell:
        """
        (while read -r line; do bedtools coverage -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy5_nchr.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -counts; done < {input.bed}) > {output}
        """


# see cadd2
rule LINSIGHT_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/linsight/LINSIGHT_hg38_sort.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_linsight_sum.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o sum; done < {input.bed}) > {output}

        """


# tabix {input.anno} -R {input.merg} | awk '!seen[$0]++' | bedtools map -a {input.bed} -b stdin -c 4 -o sum > {output}


# enhancer promotor links
rule EP_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        ep="annotations/enhancer-promoter-links/sorted_encode.bed",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_100bpup_EP.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.ep} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


# frequently interacting regulatory elements
rule fire_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/FIRE/fire_{celllines}.bed",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/{set}_100bpup_fire_{celllines}.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4,4 -o max,min > {output.t1}
        """


rule fire2_100bpup:
    input:
        expand("{{set}}/{{set}}_100bpup_fire_{fire}.bed", fire=CL),
    output:
        temp("{set}/{set}_100bpup_fire.bed"),
    shell:
        "paste {input} > {output}"


# hic data hESC hg18
rule HIC_hESC_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        hic="annotations/hic/hESC/combined/sorted.total.combined.domain",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_100bpup_HIC_hESC.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


# hic data from encode
rule HIC_encode_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        hic="annotations/Encode-HIC/{cells}/sorted_{tad}.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_100bpup_encode_{cells}_{tad}_hic.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


# merging encode HIC data
rule mergetad_100bpup:
    input:
        expand(
            "{{set}}/{{set}}_100bpup_encode_{cells}_{tad}_hic.bed", cells=CELLS, tad=TAD
        ),
    output:
        temp("{set}/{set}_100bpup_HIC.bed"),
    shell:
        "paste {input} > {output}"


# micro syntenic regions
rule microsynteny_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        hic="annotations/synteny/microsynteny.bed",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_100bpup_microsynteny.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


rule ccr_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/ccr/ccrs.all.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_ccr_mean.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4 -o max > {output}
        """


# directionality index of HIC data from genomegitar database
rule genomegitar1_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/genomegitar/{gg}/DI_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/{set}_100bpup_genomegitar_{gg}.bed"),
    shell:
        """
        bedtools map -a {input.bed} \
        -b {input.anno} -c 4,4 -o max,min > {output.t1}
        """


rule genomegitar2_100bpup:
    input:
        expand("{{set}}/{{set}}_100bpup_genomegitar_{gg}.bed", gg=genomegitars),
    output:
        temp("{set}/{set}_100bpup_DI.bed"),
    shell:
        "paste {input} > {output}"


rule genomegitar3_100bpup:
    input:
        "{set}/{set}_100bpup_DI.bed",
    output:
        maxgg=temp("{set}/{set}_100bpup_DI_max.bed"),
        mingg=temp("{set}/{set}_100bpup_DI_min.bed"),
    shell:
        """
        cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i<m)m=$i;print m}}' > {output.mingg}
        cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i;print m}}' > {output.maxgg}
        """


rule MPC_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/MPC/transcript_constraints_hg38liftover.bg.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_MPC_mean.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4 -o mean > {output}
        """


rule RemapTF_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/ReMap/ReMap2_overlapTF_hg38.bg.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_remapTF_mean.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4_nchr.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o mean; done < {input.bed}) > {output}
        """


# various regulatory and epigenetic features from ENCODE (see ENCODE list)
rule encode_100bpup:
    input:
        merg="beds/{set}/{set}_100bpup_nchr_merged.bed",
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/encode/{encodes}/{encodes}_merged_90quant.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_encode_{encodes}_mean.bed"),
    shell:
        "tabix {input.anno} -R {input.merg} | cat annotations/dummy5_nchr.bed - | bedtools map -a {input.bed} -b - -c 4,4 -o max,sum > {output}"


rule encode2_100bpup:
    input:
        expand("{{set}}/{{set}}_100bpup_encode_{encodes}_mean.bed", encodes=ENCODES),
    output:
        temp("{set}/{set}_100bpup_encode.bed"),
    shell:
        "paste {input} > {output}"


# genome states infered from chromHMM
rule chromHMM_MAX_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup_nchr.bed",
        anno="annotations/chromhmm/chromHMM_GRCh38.bg.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_chromHMM_max.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(cat annotations/dummy_chrhmm.bed; tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 -o max; done < {input.bed}) > {output}
        """


#  bedtools map -a {input.bed} -b stdin -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 -o max > {output}


rule Fantom5_counts_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/fantom5/F5.hg38.enhancers_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_f5_counts.bed"),
    shell:
        """
        bedtools coverage -a {input.bed} -b {input.anno} -counts  > {output}
        """


rule HI_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/DDD_HI/hg38_HI_Predictions_version3_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_dddhi.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 5 -o max > {output}
        """


rule deepc_100bpup:
    input:
        bed="beds/{set}/{set}_100bpup.bed",
        anno="annotations/deepc/saliencies_merged_gm12878_5kb_10bp.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpup_deepc.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o max; done < {input.bed}) > {output}
        """


rule complete_script_100bpup:
    input:
        cadd="{set}/{set}_100bpup_CADD_PC_PhyloP_maxsum.bed",
        cadd2="{set}/{set}_100bpup_cadd2_count.bed",
        ccr="{set}/{set}_100bpup_ccr_mean.bed",
        chromHMM="{set}/{set}_100bpup_chromHMM_max.bed",
        ctcf="{set}/{set}_100bpup_ctcf.bed",
        di_min="{set}/{set}_100bpup_DI_min.bed",
        di_max="{set}/{set}_100bpup_DI_max.bed",
        encode="{set}/{set}_100bpup_encode.bed",
        ep="{set}/{set}_100bpup_EP.bed",
        fire="{set}/{set}_100bpup_fire.bed",
        gc="{set}/{set}_100bpup_gc.bed",
        gm="{set}/{set}_100bpup_genemodel.bed",
        gerp="{set}/{set}_100bpup_gerp_mean.bed",
        gerp2="{set}/{set}_100bpup_gerp2_count.bed",
        hic="{set}/{set}_100bpup_HIC.bed",
        hesc="{set}/{set}_100bpup_HIC_hESC.bed",
        microsyn="{set}/{set}_100bpup_microsynteny.bed",
        mpc="{set}/{set}_100bpup_MPC_mean.bed",
        pli="{set}/{set}_100bpup_pli.bed",
        remapTF="{set}/{set}_100bpup_remapTF_mean.bed",
        f5="{set}/{set}_100bpup_f5_counts.bed",
        hi="{set}/{set}_100bpup_dddhi.bed",
        deepc="{set}/{set}_100bpup_deepc.bed",
        ultrac="{set}/{set}_100bpup_ultraconserved.bed",
        g_dist="{set}/{set}_100bpup_genetic_dist.bed",
        linsight="{set}/{set}_100bpup_linsight_sum.bed",
        header="annotations/header.txt",
    output:
        "{set}/matrix_100bpup.bed",
    conda:
        "envs/SV.yml"
    shell:
        """
        paste <(cut -f1-11 {input.cadd}) <(cut -f4 {input.cadd2}) <(cut -f4 {input.ccr}) <(cut -f4-28 {input.chromHMM}) <(cut -f4,5,7 {input.ctcf}) <(cut -f1 {input.di_min}) <(cut -f1 {input.di_max}) <(cut -f4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65 {input.encode}) <(cut -f4-7 {input.ep}) <(cut -f4,5,9,10,14,15,19,20,24,25 {input.fire}) <(cut -f4 {input.gc}) <(cut -f4-11 {input.gm}) <(cut -f4 {input.gerp}) <(cut -f4 {input.gerp2}) <(cut -f4,5,6,7,11,12,13,14,18,19,20,21,25,26,27,28 {input.hic}) <(cut -f4-7 {input.hesc}) <(cut -f4-7 {input.microsyn}) <(cut -f4 {input.mpc}) <(cut -f4 {input.pli}) <(cut -f4,5,6 {input.g_dist}) <(cut -f4 {input.remapTF}) <(cut -f4 {input.f5}) <(cut -f4 {input.hi}) <(cut -f4 {input.deepc}) <(cut -f4,5,7 {input.ultrac}) <(cut -f4 {input.linsight})| cat {input.header} - > {output}
        """


###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


rule prep_chr_100bpdown:
    input:
        bed="beds/{set}/{set}_wchr.bed",
        genome="annotations/hs38.fa.genome",
    output:
        wchr="beds/{set}/{set}_100bpdown.bed",
        nchr="beds/{set}/{set}_100bpdown_nchr.bed",
    conda:
        "envs/SV.yml"
    shell:
        """
        bedtools flank -i {input.bed} -g {input.genome} -l 0 -r 100 | bedtools sort > {output.wchr}
        bedtools flank -i {input.bed} -g {input.genome} -l 0 -r 100 | bedtools sort | sed 's/^chr\|%$//g' > {output.nchr}

        """


rule prep_merg1_100bpdown:
    input:
        nchr="beds/{set}/{set}_100bpdown_nchr.bed",
        wchr="beds/{set}/{set}_100bpdown.bed",
    conda:
        "envs/SV.yml"
    output:
        nchr="beds/{set}/{set}_100bpdown_nchr_merged.bed",
        wchr="beds/{set}/{set}_100bpdown_merged.bed",
    shell:
        """
        bedtools merge -i {input.nchr} > {output.nchr}
        bedtools merge -i {input.wchr} > {output.wchr}

        """


# CTCF Peaks from wang_et_al
rule CTCF_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/CTCF/Wangetal_hg38.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_ctcf.bed"),
    shell:
        """
        bedtools coverage -b {input.anno} -a {input.bed} > {output}
        """


rule ultraconserved_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/ultraconserved/ultraconserved_hg38_muliz120M_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_ultraconserved.bed"),
    shell:
        """
        bedtools coverage -b {input.anno} -a {input.bed} > {output}
        """


# gc content from ucsc
rule GC_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/GC/gc5Base.bedGraph.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_gc.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o mean; done < {input.bed}) > {output}
        """


# ensemble genome build temporary genenames and attribute extraction
rule gene_model_tmp_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        merg="beds/{set}/{set}_100bpdown_nchr_merged.bed",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/gm_tmp_100bpdown.bed"),
    shell:
        """
        tabix {input.anno} -R {input.merg} | awk '{{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \\"\\";"; }}' | gtf2bed - | cut -f1,2,3,8,10  > {output.t1}
        """


# genemodels
rule gene_model_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="{set}/gm_tmp_100bpdown.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_genemodel.bed"),
    shell:
        """
        paste <( bedtools coverage -b <(grep "exon" {input.anno}) -a {input.bed} | cut -f 1,2,3,5 ) <(bedtools coverage -b <(grep "transcript" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "gene" {input.anno} ) -a {input.bed} | cut -f 5) <( bedtools coverage -b <( grep "start_codon" {input.anno} ) -a {input.bed} | cut -f 5 ) <( bedtools coverage -b <(grep "stop_codon" {input.anno}) -a {input.bed} | cut -f 5) <(bedtools coverage -b <(grep "three_prime_utr" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "five_prime_utr" {input.anno}) -a {input.bed} | cut -f 5 ) <(bedtools coverage -b <(grep "CDS" {input.anno}) -a {input.bed} | cut -f 5 ) > {output}
        """


# grep "exon" {input.anno} | bedtools coverage -b stdin -a {input.bed} |cut -f 1,2,3,5 |paste - <(grep "transcript" {input.anno} |bedtools coverage -b stdin -a {input.bed} |cut -f 5 ) |paste - <(grep "gene" {input.anno} | bedtools coverage -b stdin -a {input.bed} | cut -f 5)|paste - <(grep "start_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "stop_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5) |paste - <(grep "three_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "five_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "CDS" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 ) >> {output}   """


rule gene_model_dist_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_genetic_dist.bed"),
    shell:
        """
        grep "exon" {input.anno} | bedtools closest -d -t first -b stdin -a {input.bed} |cut -f 1,2,3,9 |paste - <(grep "gene" {input.anno} | bedtools closest -d -t first -b stdin -a {input.bed} | cut -f 9)|paste - <(grep "start_codon" {input.anno} |bedtools closest -d -t first -b stdin -a {input.bed} | cut -f 9 ) >> {output}   """


# gene names necessary for PlI extraction
rule gene_names_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.GENES.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_genenames.bed"),
    shell:
        """
        bedtools map -b {input.anno} -a {input.bed} -c 4 -o distinct > {output}
        """


# PLi
rule pli_100bpdown:
    input:
        gn="{set}/{set}_100bpdown_genenames.bed",
        pli="annotations/gnomad/pli_exac.csv",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_pli.bed"),
    shell:
        """
        Rscript --vanilla scripts/PLIextract.R {input.gn} {input.pli} {output}
        """


# retrieving mean and max CADD scores from cadd_10score.bed.gz
rule cadd_PC_phylop_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/PhastCons/CADD_PC_PhyloP_scores.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_CADD_PC_PhyloP_maxsum.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}')| cat annotations/dummy8.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4,4,5,5,6,6,7,7,8,8 -o max,sum,max,sum,max,sum,max,sum,max,sum; done < {input.bed}) > {output}
        """


rule cadd2_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/CADD/CADD_GRCh38-v1.5.bedGraph_90q_12.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_cadd2_count.bed"),
    shell:
        """
        (while read -r line; do bedtools coverage -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}')| cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -counts; done < {input.bed}) > {output}
        """


# tabix {input.anno} -R {input.merg} | bedtools coverage -a {input.bed} -b stdin -counts > {output}
# merg="beds/{set}/{set}_merged.bed",


rule gerp_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_gerp_mean.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}')| cat annotations/dummy5_nchr.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o max; done < {input.bed}) > {output}
        """


rule gerp2_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_gerp2_count.bed"),
    shell:
        """
        (while read -r line; do bedtools coverage -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy5_nchr.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -counts; done < {input.bed}) > {output}
        """


# see cadd2
rule LINSIGHT_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/linsight/LINSIGHT_hg38_sort.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_linsight_sum.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o sum; done < {input.bed}) > {output}

        """


# tabix {input.anno} -R {input.merg} | awk '!seen[$0]++' | bedtools map -a {input.bed} -b stdin -c 4 -o sum > {output}


# enhancer promotor links
rule EP_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        ep="annotations/enhancer-promoter-links/sorted_encode.bed",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_100bpdown_EP.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.ep} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


# frequently interacting regulatory elements
rule fire_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/FIRE/fire_{celllines}.bed",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/{set}_100bpdown_fire_{celllines}.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4,4 -o max,min > {output.t1}
        """


rule fire2_100bpdown:
    input:
        expand("{{set}}/{{set}}_100bpdown_fire_{fire}.bed", fire=CL),
    output:
        temp("{set}/{set}_100bpdown_fire.bed"),
    shell:
        "paste {input} > {output}"


# hic data hESC hg18
rule HIC_hESC_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        hic="annotations/hic/hESC/combined/sorted.total.combined.domain",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_100bpdown_HIC_hESC.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


# hic data from encode
rule HIC_encode_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        hic="annotations/Encode-HIC/{cells}/sorted_{tad}.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_100bpdown_encode_{cells}_{tad}_hic.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


# merging encode HIC data
rule mergetad_100bpdown:
    input:
        expand(
            "{{set}}/{{set}}_100bpdown_encode_{cells}_{tad}_hic.bed",
            cells=CELLS,
            tad=TAD,
        ),
    output:
        temp("{set}/{set}_100bpdown_HIC.bed"),
    shell:
        "paste {input} > {output}"


# micro syntenic regions
rule microsynteny_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        hic="annotations/synteny/microsynteny.bed",
    conda:
        "envs/SV.yml"
    output:
        o1=temp("{set}/{set}_100bpdown_microsynteny.bed"),
    shell:
        """
        bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla scripts/annotateHIC.R stdin {output.o1}
        """


rule ccr_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/ccr/ccrs.all.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_ccr_mean.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4 -o max > {output}
        """


# directionality index of HIC data from genomegitar database
rule genomegitar1_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/genomegitar/{gg}/DI_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        t1=temp("{set}/{set}_100bpdown_genomegitar_{gg}.bed"),
    shell:
        """
        bedtools map -a {input.bed} \
        -b {input.anno} -c 4,4 -o max,min > {output.t1}
        """


rule genomegitar2_100bpdown:
    input:
        expand("{{set}}/{{set}}_100bpdown_genomegitar_{gg}.bed", gg=genomegitars),
    output:
        temp("{set}/{set}_100bpdown_DI.bed"),
    shell:
        "paste {input} > {output}"


rule genomegitar3_100bpdown:
    input:
        "{set}/{set}_100bpdown_DI.bed",
    output:
        maxgg=temp("{set}/{set}_100bpdown_DI_max.bed"),
        mingg=temp("{set}/{set}_100bpdown_DI_min.bed"),
    shell:
        """
        cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i<m)m=$i;print m}}' > {output.mingg}
        cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i;print m}}' > {output.maxgg}"""


rule MPC_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/MPC/transcript_constraints_hg38liftover.bg.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_MPC_mean.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 4 -o mean > {output}
        """


rule RemapTF_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/ReMap/ReMap2_overlapTF_hg38.bg.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_remapTF_mean.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4_nchr.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o mean; done < {input.bed}) > {output}
        """


# various regulatory and epigenetic features from ENCODE (see ENCODE list)
rule encode_100bpdown:
    input:
        merg="beds/{set}/{set}_100bpdown_nchr_merged.bed",
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/encode/{encodes}/{encodes}_merged_90quant.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_encode_{encodes}_mean.bed"),
    shell:
        "tabix {input.anno} -R {input.merg} | cat annotations/dummy5_nchr.bed - | bedtools map -a {input.bed} -b - -c 4,4 -o max,sum > {output}"


rule encode2_100bpdown:
    input:
        expand("{{set}}/{{set}}_100bpdown_encode_{encodes}_mean.bed", encodes=ENCODES),
    output:
        temp("{set}/{set}_100bpdown_encode.bed"),
    shell:
        "paste {input} > {output}"


# genome states infered from chromHMM
rule chromHMM_MAX_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown_nchr.bed",
        anno="annotations/chromhmm/chromHMM_GRCh38.bg.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_chromHMM_max.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(cat annotations/dummy_chrhmm.bed; tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 -o max; done < {input.bed}) > {output}
        """


#  bedtools map -a {input.bed} -b stdin -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 -o max > {output}


rule Fantom5_counts_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/fantom5/F5.hg38.enhancers_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_f5_counts.bed"),
    shell:
        """
        bedtools coverage -a {input.bed} -b {input.anno} -counts  > {output}
        """


rule HI_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/DDD_HI/hg38_HI_Predictions_version3_sort.bed",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_dddhi.bed"),
    shell:
        """
        bedtools map -a {input.bed} -b {input.anno} -c 5 -o max > {output}
        """


rule deepc_100bpdown:
    input:
        bed="beds/{set}/{set}_100bpdown.bed",
        anno="annotations/deepc/saliencies_merged_gm12878_5kb_10bp.bed.gz",
    conda:
        "envs/SV.yml"
    output:
        temp("{set}/{set}_100bpdown_deepc.bed"),
    shell:
        """
        (while read -r line; do bedtools map -b <(tabix {input.anno} $(echo $line | awk '{{ print $1":"$2"-"$3+1}}') | cat annotations/dummy4.bed - ) -a <(echo $line | awk 'BEGIN{{ OFS="\t" }}{{ print $1,$2,$3}}') -c 4 -o max; done < {input.bed}) > {output}
        """


rule complete_script_100bpdown:
    input:
        cadd="{set}/{set}_100bpdown_CADD_PC_PhyloP_maxsum.bed",
        cadd2="{set}/{set}_100bpdown_cadd2_count.bed",
        ccr="{set}/{set}_100bpdown_ccr_mean.bed",
        chromHMM="{set}/{set}_100bpdown_chromHMM_max.bed",
        ctcf="{set}/{set}_100bpdown_ctcf.bed",
        di_min="{set}/{set}_100bpdown_DI_min.bed",
        di_max="{set}/{set}_100bpdown_DI_max.bed",
        encode="{set}/{set}_100bpdown_encode.bed",
        ep="{set}/{set}_100bpdown_EP.bed",
        fire="{set}/{set}_100bpdown_fire.bed",
        gc="{set}/{set}_100bpdown_gc.bed",
        gm="{set}/{set}_100bpdown_genemodel.bed",
        gerp="{set}/{set}_100bpdown_gerp_mean.bed",
        gerp2="{set}/{set}_100bpdown_gerp2_count.bed",
        hic="{set}/{set}_100bpdown_HIC.bed",
        hesc="{set}/{set}_100bpdown_HIC_hESC.bed",
        microsyn="{set}/{set}_100bpdown_microsynteny.bed",
        mpc="{set}/{set}_100bpdown_MPC_mean.bed",
        pli="{set}/{set}_100bpdown_pli.bed",
        remapTF="{set}/{set}_100bpdown_remapTF_mean.bed",
        f5="{set}/{set}_100bpdown_f5_counts.bed",
        hi="{set}/{set}_100bpdown_dddhi.bed",
        deepc="{set}/{set}_100bpdown_deepc.bed",
        ultrac="{set}/{set}_100bpdown_ultraconserved.bed",
        g_dist="{set}/{set}_100bpdown_genetic_dist.bed",
        linsight="{set}/{set}_100bpdown_linsight_sum.bed",
        header="annotations/header.txt",
    output:
        "{set}/matrix_100bpdown.bed",
    conda:
        "envs/SV.yml"
    shell:
        """
        paste <(cut -f1-11 {input.cadd}) <(cut -f4 {input.cadd2}) <(cut -f4 {input.ccr}) <(cut -f4-28 {input.chromHMM}) <(cut -f4,5,7 {input.ctcf}) <(cut -f1 {input.di_min}) <(cut -f1 {input.di_max}) <(cut -f4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65 {input.encode}) <(cut -f4-7 {input.ep}) <(cut -f4,5,9,10,14,15,19,20,24,25 {input.fire}) <(cut -f4 {input.gc}) <(cut -f4-11 {input.gm}) <(cut -f4 {input.gerp}) <(cut -f4 {input.gerp2}) <(cut -f4,5,6,7,11,12,13,14,18,19,20,21,25,26,27,28 {input.hic}) <(cut -f4-7 {input.hesc}) <(cut -f4-7 {input.microsyn}) <(cut -f4 {input.mpc}) <(cut -f4 {input.pli}) <(cut -f4,5,6 {input.g_dist}) <(cut -f4 {input.remapTF}) <(cut -f4 {input.f5}) <(cut -f4 {input.hi}) <(cut -f4 {input.deepc}) <(cut -f4,5,7 {input.ultrac}) <(cut -f4 {input.linsight})| cat {input.header} - > {output}
        """


###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


rule scoring:
    input:
        span="{set}/matrix.bed",
        flank_up="{set}/matrix_100bpdown.bed",
        flank_down="{set}/matrix_100bpup.bed",
    conda:
        "envs/SV.yml"
    output:
      temp("output/{set}.score"),
    params:
        name="{set}"
    shell:
        """
        Rscript --vanilla scripts/scoring.R {params.name} {input.span} {input.flank_up} {input.flank_down} {output}
        """



rule sort:
    input:
      score="output/{set}.score",
      header="annotations/header_final.txt",
    conda:
        "envs/SV.yml"
    output:
        "output/{set}_score.bed"
    shell:
        """
        bedtools sort -i {input.score} | cat {input.header} - > {output}
        """





