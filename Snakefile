# Snakefile
# Annotation script for feature retrieval and scoring using the CADD-SV model
# CADD-SV

SETS=["test"]#,"sim_cdel2","sim_cdel3","sim_cdel4","sim_cdel5","sim_cdel6","sim_cdel7","sim_cdel8","sim_cdel9","sim_cdel"] #"DDD4","DDD3","sim_cdel1","sim_cdel2","sim_cdel3","sim_cdel4","sim_cdel5","sim_cdel6","sim_cdel7","sim_cdel8","sim_cdel9","cdel","sim_cdel","hins","benigndel","pathodel_smaller1mb","hdel500","cins500","sim_cdel","sim_hins","sim_cins500","sim_hdel500","pathodel_smaller1mb_downstream","pathodel_smaller1mb_upstream","gnomad-sv_DEL_1","gnomad-sv_DEL_1_upstream","gnomad-sv_DEL_1_downstream","benigndel_upstream","benigndel_downstream","test","pathoins500","benignins500","pathodup_smaller1mb","benigndup_smaller50kb"]    #"gnomad-sv_DEL","benigndel","sim_cdel","cdel","pathodel","gnomad-sv_DEL_1", ,"hins","hdel","sim_cins","sim_hdel","cins","cdel","sim_hins","benigndel","pathodel","gnomad-sv_DEL_1"

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
    matrix=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/score.bed",sets=SETS)
     
#CTCF Peaks from wang_et_al
rule CTCF:
  input:
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="../annotations/CTCF/Wangetal_hg38.bed"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_ctcf.bed"
  shell: """
    bedtools coverage -b {input.anno} -a {input.bed} > {output}
    """

#gc content from ucsc
rule GC:
  input:
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="../annotations/GC/gc5Base.bedGraph.gz"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_gc.bed"
  shell: """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o sum | awk '{{print $1 "\t" ($2) "\t" ($2) "\t" ($4*5/($3-$2)) }}' > {output}
     """

#ensemble genome build temporary genenames and attribute extraction 
rule gene_model_tmp:
  input:
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    anno="../annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz"
  conda: "envs/SV.yml"
  output: 
      t1=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/gm_tmp.bed")
  shell: """
    tabix {input.anno} -R {input.merg} | awk '{{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \\"\\";"; }}' | gtf2bed - | cut -f1,2,3,8,10  > {output.t1}
    """
# genemodels 
rule gene_model:
  input:
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/gm_tmp.bed"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genemodel.bed"
  shell: """
    grep "exon" {input.anno} | bedtools coverage -b stdin -a {input.bed} |cut -f 1,2,3,5 |paste - <(grep "transcript" {input.anno} |bedtools coverage -b stdin -a {input.bed} |cut -f 5 ) |paste - <(grep "gene" {input.anno} | bedtools coverage -b stdin -a {input.bed} | cut -f 5)|paste - <(grep "start_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "stop_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5) |paste - <(grep "three_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "five_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "CDS" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 ) >> {output}   """
  
# gene names necessary for PlI extraction
rule gene_names:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/gm_tmp.bed"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genenames.bed"
  shell: """
    cut -f 1,2,3 {input.anno} | paste - <(cut -f 5 {input.anno} | tr ";" "\n" | grep -E 'gene_name' | cut -f 1 | tr " " "\t" | cut -f 3 ) | bedtools map -b stdin -a {input.bed} -c 4 -o distinct >{output}
    """

#PLi
rule pli:
  input:  
    gn="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genenames.bed",
    pli="../annotations/gnomad/pli_exac.csv"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_pli.bed"
  shell:  """
    Rscript --vanilla ../scripts/PLIextract.R {input.gn} {input.pli} {output}
    """

#retrieving mean and max CADD scores from cadd_10score.bed.gz
rule cadd1:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    anno="../annotations/CADD/cadd_10score.bed.gz"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_cadd_meanmax.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} -b stdin -c 4,5,4,5 -o max,max,mean,mean > {output}
    """

#enhancer promotor links
rule EP:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    ep="../annotations/enhancer-promoter-links/sorted_encode.bed"
  conda: "envs/SV.yml"
  output: 
    o1="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_EP.bed"
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.ep} | Rscript --vanilla ../scripts/annotateHIC.R stdin {output.o1}
    """

#frequently interacting regulatory elements
rule fire:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{sets}_nochr.bed",
    anno="../annotations/FIRE/fire_{celllines}.bed"
  conda: "envs/SV.yml"
  output: t1=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_fire_{celllines}.bed")
  shell:  """
    bedtools map -a {input.bed} -b {input.anno} -c 4,4 -o max,min > {output.t1}
    """
  
def mergeFire(wc):
  return(expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_fire_{fire}.bed",set=wc.set,fire=CL))

rule fire2:
  input: mergeFire
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_fire.bed"
  conda: "envs/rules.yml"
  shell: "paste {input} > {output}"
  
  
#hic data hESC hg18
rule HIC_hESC:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    hic="../annotations/hic/hESC/combined/sorted.total.combined.domain"
  conda: "envs/SV.yml"
  output: 
    o1="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_HIC_hESC.bed"
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla ../scripts/annotateHIC.R stdin {output.o1}
     """

#hic data from encode     
rule HIC_encode:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{sets}.bed",
    hic="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/Encode-HIC/{cells}/sorted_{tad}.bed.gz"
  conda: "envs/SV.yml"
  output: 
    o1=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_{cells}_{tad}.bed")
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla ../scripts/annotateHIC.R stdin {output.o1}
    """

#merging encode HIC data
def mergeHic(wc):
  return(expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_{cells}_{tad}.bed",set=wc.set,cells=CELLS,tad=TAD))

rule mergetad:
  input: mergeHic
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_HIC.bed"
  conda: "envs/rules.yml"
  shell: "paste {input} > {output}"

#micro syntenic regions
rule microsynteny:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{sets}_nochr.bed",
    hic="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/synteny/microsynteny.bed"
  conda: "envs/SV.yml"
  output: 
    o1="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_microsynteny.bed"
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} | Rscript --vanilla ../scripts/annotateHIC.R stdin {output.o1}
    """

#evolutionary conservation scores
rule phastcons:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/PhastCons/hg38.phastCons{way}way.bedGraph.gz"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_PhastCons{way}_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg}  | bedtools map -a {input.bed} -b stdin -c 4,4 -o max,mean > {output}"
    
rule phyloP:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/PhastCons/hg38.phyloP20way.bedGraph.gz"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_phyloP20_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg}  | bedtools map -a {input.bed} -b stdin -c 4,4 -o max,mean > {output}"
   
rule gerp:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="../annotations/gerp/gerp_score2_hg38_MAM.bg.gz"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_gerp_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} -b stdin -c 4 -o mean > {output}"

rule ccr:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="../annotations/ccr/ccrs.all.bed.gz"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_ccr_mean.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o mean > {output}
    """

# directionality index of HIC data from genomegitar database  
rule genomegitar1:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/genomegitar/{gg}/DI_sort.bed"
  conda: "envs/SV.yml"
  output: 
    t1=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genomegitar_{gg}.bed")
  shell:  """
    bedtools map -a {input.bed} \
    -b {input.anno} -c 4,4 -o max,min > {output.t1}
    """
    
def mergeGgitars(wc):
  return(expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genomegitar_{gg}.bed",set=wc.set,gg=genomegitars))

rule genomegitar2:
  input: mergeGgitars
  output: temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_DI.bed")
  conda: "envs/rules.yml"
  shell: "paste {input} > {output}"

rule genomegitar3:
  input: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_DI.bed"
  output: 
    maxgg="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_DI_max.bed",
    mingg="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_DI_min.bed"
  conda: "envs/rules.yml"
  shell: """
    cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i<m)m=$i;print m}}' > {output.mingg}
    cut -f 4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 {input} | awk '{{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i;print m}}' > {output.maxgg}"""


rule MPC:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="../annotations/MPC/transcript_constraints_hg38liftover.bg.gz"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_MPC_mean.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o mean > {output}
    """
 
rule RemapTF:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="..annotations/Remap/ReMap2_overlapTF_hg38.bg.gz"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_remapTF_mean.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o mean > {output}
    """

#various regulatory and epigenetic features from ENCODE (see ENCODE list)
rule encode:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="../annotations/encode/{encodes}/{encodes}_merged.bg.gz"
  conda: "envs/SV.yml"
  output: temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_encode_{encodes}_mean.bed")
  shell:  "tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} -b stdin -c 4 -o mean > {output}"

def mergeEncode(wc):
  return(expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_encode_{encodes}_mean.bed",set=wc.set,encodes=ENCODES))

rule encode2:
  input: mergeEncode
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_encode.bed"
  conda: "envs/rules.yml"
  shell: "paste {input} > {output}"

#genome states infered from chromHMM    
rule chromHMM_MAX:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="../annotations/chromhmm/chromHMM_GRCh38.bg.gz"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_chromHMM_max.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 -o max > {output}
    """

rule complete_script:
  input:
    cadd="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_cadd_meanmax.bed",
    ccr="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_ccr_mean.bed",
    chromHMM="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_chromHMM_max.bed",
    ctcf="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_ctcf.bed",
    di_min="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_DI_min.bed",
    di_max="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_DI_max.bed",
    encode="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_encode.bed",
    ep="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_EP.bed",
    fire="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_fire.bed",
    gc="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_gc.bed",
    gm="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genemodel.bed",
    gerp="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_gerp_mean.bed",
    hic="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_HIC.bed",
    hesc="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_HIC_hESC.bed",
    microsyn="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_microsynteny.bed",
    mpc="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_MPC_mean.bed",
    pc1="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_PhastCons20_mean.bed",
    pc2="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_PhastCons30_mean.bed",
    pc3="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_PhastCons100_mean.bed",
    phyloP="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_phyloP20_mean.bed",
    pli="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_pli.bed",
    remapTF="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_remapTF_mean.bed"
  output:"/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/matrix.bed"
  conda: "envs/SV.yml"
  shell:  """
    paste <(cut -f1,2,3,4,5,6,7 {input.cadd}) <(cut -f4 {input.ccr}) <(cut -f4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 {input.chromHMM}) <(cut -f4,5,6,7 {input.ctcf}) <(cut -f1 {input.di_min}) <(cut -f1 {input.di_max}) <(cut -f4,8,12,16,20,24,28,32,36,40,44,48,52 {input.encode}) <(cut -f4,5,6,7 {input.ep}) <(cut -f4,5,9,10,14,15,19,20,24,25 {input.fire}) <(cut -f4 {input.gc}) <(cut -f4,5,6,7,8,9,10,11 {input.gm}) <(cut -f4 {input.gerp}) <(cut -f4,5,6,7,11,12,13,14,18,19,20,21,25,26,27,28 {input.hic}) <(cut -f4,5,6,7 {input.hesc}) <(cut -f4,5,6,7 {input.microsyn}) <(cut -f4 {input.mpc}) <(cut -f4,5 {input.pc1}) <(cut -f4,5 {input.pc2}) <(cut -f4,5 {input.pc3}) <(cut -f4,5 {input.phyloP}) <(cut -f4 {input.pli})  <(cut -f4 {input.remapTF}) | cat ../annotations/header.txt - > {output}
    """
  
rule scoring:
  input: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/matrix.bed"
  conda: "envs/SV.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/score.bed"
  shell:  """
    Rscript --vanilla  ../scripts/scoring.R {input} {output}
    """
    
    
    
    
