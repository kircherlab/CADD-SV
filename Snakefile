# Snakefile
# Annotation script for feature retrieval 
# CADD-SV

SETS=["benigndel","sim_cdel","sim_cdel1","sim_cdel2","sim_cdel3","sim_cdel4","sim_cdel5","sim_cdel6","sim_cdel7","sim_cdel8","sim_cdel9","pathodel","cdel"] #"DDD4","DDD3","sim_cdel1","sim_cdel2","sim_cdel3","sim_cdel4","sim_cdel5","sim_cdel6","sim_cdel7","sim_cdel8","sim_cdel9","cdel","sim_cdel","hins","benigndel","pathodel_smaller1mb","hdel500","cins500","sim_cdel","sim_hins","sim_cins500","sim_hdel500","pathodel_smaller1mb_downstream","pathodel_smaller1mb_upstream","gnomad-sv_DEL_1","gnomad-sv_DEL_1_upstream","gnomad-sv_DEL_1_downstream","benigndel_upstream","benigndel_downstream","test","pathoins500","benignins500","pathodup_smaller1mb","benigndup_smaller50kb"]    #"gnomad-sv_DEL","benigndel","sim_cdel","cdel","pathodel","gnomad-sv_DEL_1", ,"hins","hdel","sim_cins","sim_hdel","cins","cdel","sim_hins","benigndel","pathodel","gnomad-sv_DEL_1"

CELLS=["A549","Caki2"]
TAD=["nested","tad"]

#set of encode datasets in annotation/
ENCODES=["DNase-seq","H2AFZ","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me2","H3K4me3","H3K79me2","H3K9ac","H3K9me3","H4K20me1","totalRNA-seq"] 

#set of genomegitar datasets for HIC directionality index
genomegitars=["GSM1055800_DI","GSM1055805_DI","GSM1081530_DI","GSM1267196_DI","GSM1267200_DI","GSM1294038_DI","GSM1294039_DI","GSM1551599_DI","GSM1551629_DI","GSM1608505_DI","GSM1718021_DI","GSM1906332_DI","GSM1906333_DI","GSM1906334_DI","GSM1909121_DI","GSM455133_DI","GSM862723_DI","GSM862724_DI","GSM927075_DI"]

#Cell lines from Encode for HIC datasets  
CL=["gm12878","msc","mes","imr90","h1"]

rule all:
  input:
    matrix=expand("{sets}/score.bed",sets=SETS)
     
#CTCF Peaks from wang_et_al
rule CTCF:
  input:
    bed="../beds/{set}.bed",
    anno="../dependencies/CTCF/Wangetal_hg38.bed"
  conda: "envs/SV.yml"
  output: "{set}/{set}_ctcf.bed"
  shell: """
    bedtools coverage -b {input.anno} -a {input.bed} > {output}
    """

#gc content from ucsc
rule GC:
  input:
    merg="../beds/{set}_merged.bed",
    bed="../beds/{set}.bed",
    anno="../dependencies/GC/gc5Base.bedGraph.gz"
  conda: "envs/SV.yml"
  output: "{set}/{set}_gc.bed"
  shell: """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o sum | awk '{{print $1 "\t" ($2) "\t" ($2) "\t" ($4*5/($3-$2)) }}' > {output}
     """

#ensemble genome build temporary genenames and attribute extraction 
rule gene_model_tmp:
  input:
    bed="../beds/{set}_nochr.bed",
    merg="../beds/{set}_nochr_merged.bed",
    anno="../dependencies/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz"
  conda: "envs/SV.yml"
  output: 
      t1=temp("{set}/gm_tmp.bed")
  shell: """
    tabix {input.anno} -R {input.merg} | awk '{{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \\"\\";"; }}' | gtf2bed - | cut -f1,2,3,8,10  > {output.t1}
    """
# genemodels 
rule gene_model:
  input:
    bed="../beds/{set}_nochr.bed",
    anno="{set}/gm_tmp.bed"
  conda: "envs/SV.yml"
  output: "{set}/{set}_genemodel.bed"
  shell: """
    grep "exon" {input.anno} | bedtools coverage -b stdin -a {input.bed} |cut -f 1,2,3,5 |paste - <(grep "transcript" {input.anno} |bedtools coverage -b stdin -a {input.bed} |cut -f 5 ) |paste - <(grep "gene" {input.anno} | bedtools coverage -b stdin -a {input.bed} | cut -f 5)|paste - <(grep "start_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "stop_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5) |paste - <(grep "three_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "five_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "CDS" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 ) >> {output}   """
  
# gene names necessary for PlI extraction
rule gene_names:
  input:  
    bed="../beds/{set}_nochr.bed",
    anno="{set}/gm_tmp.bed"
  conda: "envs/SV.yml"
  output: "{set}/{set}_genenames.bed"
  shell: """
    cut -f 1,2,3 {input.anno} | paste - <(cut -f 5 {input.anno} | tr ";" "\n" | grep -E 'gene_name' | cut -f 1 | tr " " "\t" | cut -f 3 ) | bedtools map -b stdin -a {input.bed} -c 4 -o distinct >{output}
    """

#PLi
rule pli:
  input:  
    gn="{set}/{set}_genenames.bed",
    pli="../dependencies/gnomad/pli_exac.csv"
  conda: "envs/SV.yml"
  output: "{set}/{set}_pli.bed"
  shell:  """
    Rscript --vanilla ../scripts/PLIextract.R {input.gn} {input.pli} {output}
    """

#retrieving mean and max CADD scores from cadd_10score.bed.gz
rule cadd_PC_phylop:
  input:  
    bed="../beds/{set}.bed",
    merg="../beds/{set}_merged.bed",
    anno="../dependencies/PhastCons/CADD_PC_PhyloP_scores.bed.gz"
  conda: "envs/SV.yml"
  output: "{set}/{set}_CADD_PC_PhyloP_maxsum.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} -b stdin -c 4,4,5,5,6,6,7,7,8,8 -o max,sum,max,sum,max,sum,max,sum,max,sum > {output}
    """
   
rule gerp:
  input:  
    merg="../beds/{set}_nochr_merged.bed",
    bed="../beds/{set}_nochr.bed",
    anno="../dependencies/gerp/gerp_score2_hg38_MAM.bg.gz"
  conda: "envs/SV.yml"
  output: "{set}/{set}_gerp_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} -b stdin -c 4 -o max > {output}"


#enhancer promotor links
rule EP:
  input:  
    bed="../beds/{set}.bed",
    ep="../dependencies/enhancer-promoter-links/sorted_encode.bed"
  conda: "envs/SV.yml"
  output: 
    o1="{set}/{set}_EP.bed"
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.ep} | Rscript --vanilla ../scripts/annotateHIC.R stdin {output.o1}
    """

#frequently interacting regulatory elements
rule fire:
  input:  
    bed="../beds/{sets}_nochr.bed",
    anno="../dependencies/FIRE/fire_{celllines}.bed"
  conda: "envs/SV.yml"
  output: t1=temp("{sets}/{sets}_fire_{celllines}.bed")
  shell:  """
    bedtools map -a {input.bed} -b {input.anno} -c 4,4 -o max,min > {output.t1}
    """
  
def mergeFire(wc):
  return(expand("{set}/{set}_fire_{fire}.bed",set=wc.set,fire=CL))

rule fire2:
  input: mergeFire
  output: "{set}/{set}_fire.bed"
  shell: "paste {input} > {output}"
  
  
#hic data hESC hg18
rule HIC_hESC:
  input:  
    bed="../beds/{set}.bed",
    hic="../dependencies/hic/hESC/combined/sorted.total.combined.domain"
  conda: "envs/SV.yml"
  output: 
    o1="{set}/{set}_HIC_hESC.bed"
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla ../scripts/annotateHIC.R stdin {output.o1}
     """

#hic data from encode     
rule HIC_encode:
  input:  
    bed="../beds/{sets}.bed",
    hic="../dependencies/Encode-HIC/{cells}/sorted_{tad}.bed.gz"
  conda: "envs/SV.yml"
  output: 
    o1=temp("{sets}/{sets}_encode_{cells}_{tad}_hic.bed")
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla ../scripts/annotateHIC.R stdin {output.o1}
    """

#merging encode HIC data
def mergeHic(wc):
  return(expand("{set}/{set}_encode_{cells}_{tad}_hic.bed",set=wc.set,cells=CELLS,tad=TAD))

rule mergetad:
  input: mergeHic
  output: "{set}/{set}_HIC.bed"
  shell: "paste {input} > {output}"

#micro syntenic regions
rule microsynteny:
  input:  
    bed="../beds/{sets}_nochr.bed",
    hic="../dependencies/synteny/microsynteny.bed"
  conda: "envs/SV.yml"
  output: 
    o1="{sets}/{sets}_microsynteny.bed"
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla ../scripts/annotateHIC.R stdin {output.o1}
    """

rule ccr:
  input:  
    merg="../beds/{set}_nochr_merged.bed",
    bed="../beds/{set}_nochr.bed",
    anno="../dependencies/ccr/ccrs.all.bed.gz"
  conda: "envs/SV.yml"
  output: "{set}/{set}_ccr_mean.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o mean > {output}
    """

# directionality index of HIC data from genomegitar database  
rule genomegitar1:
  input:  
    bed="../beds/{set}.bed",
    anno="../dependencies/genomegitar/{gg}/DI_sort.bed"
  conda: "envs/SV.yml"
  output: 
    t1=temp("{set}/{set}_genomegitar_{gg}.bed")
  shell:  """
    bedtools map -a {input.bed} \
    -b {input.anno} -c 4,4 -o max,min > {output.t1}
    """
    
def mergeGgitars(wc):
  return(expand("{set}/{set}_genomegitar_{gg}.bed",set=wc.set,gg=genomegitars))

rule genomegitar2:
  input: mergeGgitars
  output: temp("{set}/{set}_DI.bed")
  shell: "paste {input} > {output}"

rule genomegitar3:
  input: "{set}/{set}_DI.bed"
  output: 
    maxgg="{set}/{set}_DI_max.bed",
    mingg="{set}/{set}_DI_min.bed"
  shell: """
    cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i<m)m=$i;print m}}' > {output.mingg}
    cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i;print m}}' > {output.maxgg}"""


rule MPC:
  input:  
    merg="../beds/{set}_nochr_merged.bed",
    bed="../beds/{set}_nochr.bed",
    anno="../dependencies/MPC/transcript_constraints_hg38liftover.bg.gz"
  conda: "envs/SV.yml"
  output: "{set}/{set}_MPC_mean.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o mean > {output}
    """
 
rule RemapTF:
  input:  
    merg="../beds/{set}_nochr_merged.bed",
    bed="../beds/{set}_nochr.bed",
    anno="../dependencies/ReMap/ReMap2_overlapTF_hg38.bg.gz"
  conda: "envs/SV.yml"
  output: "{set}/{set}_remapTF_mean.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o mean > {output}
    """

#various regulatory and epigenetic features from ENCODE (see ENCODE list)
rule encode:
  input:  
    merg="../beds/{set}_nochr_merged.bed",
    bed="../beds/{set}_nochr.bed",
    anno="../dependencies/encode/{encodes}/{encodes}_merged_90quant.bed.gz"
  conda: "envs/SV.yml"
  output: temp("{set}/{set}_encode_{encodes}_mean.bed")
  shell:  "tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} -b stdin -c 4 -o mean > {output}"

def mergeEncode(wc):
  return(expand("{set}/{set}_encode_{encodes}_mean.bed",set=wc.set,encodes=ENCODES))

rule encode2:
  input: mergeEncode
  output: "{set}/{set}_encode.bed"
  shell: "paste {input} > {output}"

#genome states infered from chromHMM    
rule chromHMM_MAX:
  input:  
    merg="../beds/{set}_nochr_merged.bed",
    bed="../beds/{set}_nochr.bed",
    anno="../dependencies/chromhmm/chromHMM_GRCh38.bg.gz"
  conda: "envs/SV.yml"
  output: "{set}/{set}_chromHMM_max.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 -o max > {output}
    """

rule Fantom5_counts:
  input:  
    bed="../beds/{set}.bed",
    anno="../dependencies/fantom5/F5.hg38.enhancers_sort.bed"
  conda: "envs/SV.yml"
  output: "{set}/{set}_f5_counts.bed"
  shell:  """
    bedtools coverage -a {input.bed} -b {input.anno} -counts  > {output}
    """

rule HI:
  input:
    bed="../beds/{set}.bed",
    anno="../dependencies/DDD_HI/HI_Predictions_Version3_sort.bed"
  conda: "envs/SV.yml"
  output: "{set}/{set}_dddhi.bed"
  shell: """
    bedtools map -a {input.bed} -b {input.anno} -c 5 -o max > {output}
     """

rule deepc:
  input:
    merg="../beds/{set}_merged.bed",
    bed="../beds/{set}.bed",
    anno="../dependencies/deepc/saliencies_merged_gm12878_5kb_10bp.bed.gz"
  conda: "envs/SV.yml"
  output: "{set}/{set}_deepc.bed"
  shell: """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o max  > {output}
     """


rule complete_script:
  input:
    cadd="{set}/{set}_CADD_PC_PhyloP_maxsum.bed",
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
    hic="{set}/{set}_HIC.bed",
    hesc="{set}/{set}_HIC_hESC.bed",
    microsyn="{set}/{set}_microsynteny.bed",
    mpc="{set}/{set}_MPC_mean.bed",
    pli="{set}/{set}_pli.bed",
    remapTF="{set}/{set}_remapTF_mean.bed",
    f5="{set}/{set}_f5_counts.bed",
    hi="{set}/{set}_dddhi.bed",
    deepc="{set}/{set}_deepc.bed"
  output:"{set}/matrix.bed"
  conda: "envs/SV.yml"
  shell:  """
    paste <(cut -f1,2,3,4,5,6,7,8,9,10,11 {input.cadd}) <(cut -f4 {input.ccr}) <(cut -f4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 {input.chromHMM}) <(cut -f4,5,6,7 {input.ctcf}) <(cut -f1 {input.di_min}) <(cut -f1 {input.di_max}) <(cut -f4,8,12,16,20,24,28,32,36,40,44,48,52 {input.encode}) <(cut -f4,5,6,7 {input.ep}) <(cut -f4,5,9,10,14,15,19,20,24,25 {input.fire}) <(cut -f4 {input.gc}) <(cut -f4,5,6,7,8,9,10,11 {input.gm}) <(cut -f4 {input.gerp}) <(cut -f4,5,6,7,11,12,13,14,18,19,20,21,25,26,27,28 {input.hic}) <(cut -f4,5,6,7 {input.hesc}) <(cut -f4,5,6,7 {input.microsyn}) <(cut -f4 {input.mpc}) <(cut -f4 {input.pli})  <(cut -f4 {input.remapTF}) <(cut -f4 {input.f5}) <(cut -f4 {input.hi}) <(cut -f4 {input.deepc})| cat ../dependencies/header.txt - > {output}
    """
  
rule scoring:
  input: "{sets}/matrix.bed"
  conda: "envs/SV.yml"
  output: "{sets}/score.bed"
  shell:  """
    Rscript --vanilla  ../scripts/scoring.R {input} {output}
    """
    
    
    
