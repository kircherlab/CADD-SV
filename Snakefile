

SETS=["sim_cdel1","sim_cdel2","sim_cdel3","sim_cdel4","sim_cdel5","sim_cdel6","sim_cdel7","sim_cdel8","sim_cdel9"]
#["yribg","cdel","sim_cdel","hins","benigndel","pathodel_smaller1mb","hdel500","cins500","sim_cdel","sim_hins","sim_cins500","sim_hdel500","pathodel_smaller1mb_downstream","pathodel_smaller1mb_upstream","gnomad-sv_DEL_1","gnomad-sv_DEL_1_upstream","gnomad-sv_DEL_1_downstream","benigndel_upstream","benigndel_downstream","test","pathoins500","benignins500","pathodup_smaller1mb","benigndup_smaller50kb"]    #"gnomad-sv_DEL","benigndel","sim_cdel","cdel","pathodel","gnomad-sv_DEL_1", ,"hins","hdel","sim_cins","sim_hdel","cins","cdel","sim_hins","benigndel","pathodel","gnomad-sv_DEL_1"

CELLS=["A549","Caki2"]
TAD=["nested","tad"]

ENCODES=["DNase-seq","H2AFZ","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me2","H3K4me3","H3K79me2","H3K9ac","H3K9me3","H4K20me1","totalRNA-seq"] #PhastCons #,"RemapTF","cadd","MPC","mirTarget","ccr","gerp"

genomegitars=["GSM1055800_DI","GSM1055805_DI","GSM1081530_DI","GSM1267196_DI","GSM1267200_DI","GSM1294038_DI","GSM1294039_DI","GSM1551599_DI","GSM1551629_DI","GSM1608505_DI","GSM1718021_DI","GSM1906332_DI","GSM1906333_DI","GSM1906334_DI","GSM1909121_DI","GSM455133_DI","GSM862723_DI","GSM862724_DI","GSM927075_DI"]
  
CL=["gm12878","msc","mes","imr90","h1"]

rule all:
  input:
    encode_hic=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_{cells}nested.bed",sets=SETS,cells=CELLS),
    encode_tad=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_{cells}tad.bed",sets=SETS,cells=CELLS),
    encode=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_encode_{encodes}_mean.bed",sets=SETS,encodes=ENCODES),
    encode2=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_encode.bed",sets=SETS),
    chromHMM=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_chromHMM_max.bed",sets=SETS),
    pc1=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_PhastCons20_mean.bed",sets=SETS),
    pc2=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_PhastCons30_mean.bed",sets=SETS),
    pc3=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_PhastCons100_mean.bed",sets=SETS),
    phyloP=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_phyloP20_mean.bed",sets=SETS),
    #gnomadsv=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_gnomad_sv.bed",sets=SETS),
    hesc=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_HIC_hESC.bed",sets=SETS),
    ep=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_EP.bed",sets=SETS),
    gm=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_genemodel.bed",sets=SETS),
    gn=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_genenames.bed",sets=SETS),
    pli=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_pli.bed",sets=SETS),
    microsyn=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_microsynteny.bed",sets=SETS),
    cadd1=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_cadd_meanmax.bed",sets=SETS),
    genomegitar1=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_genomegitar_{gg}.bed",sets=SETS,gg=genomegitars),
    ctcf=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_ctcf.bed",sets=SETS),
    fire=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_fire_{celllines}.bed",sets=SETS,celllines=CL),
    gc=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_gc.bed",sets=SETS),
    genomegitar2=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_DI.bed",sets=SETS),
    gerp=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_gerp_mean.bed",sets=SETS),
    ccr=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_ccr_mean.bed",sets=SETS),
    remaptf=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_remapTF_mean.bed",sets=SETS),
    #mirtarget=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_mirTarget_mean.bed",sets=SETS),
    mpc=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_MPC_mean.bed",sets=SETS),
    merge1=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_fire.bed",sets=SETS),
    merge2=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_HIC.bed",sets=SETS),
    mingg=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_DI_max.bed",sets=SETS),
    maxgg=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_DI_min.bed",sets=SETS),
    matrix=expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/matrix.bed",sets=SETS)
 
#rule gnomad_sv:
#  input:  
#    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
#    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
#    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/gnomad-sv/gnomad_v2_sv.sites.bed.gz"
#  conda: "envs/rules2.yml"
#  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_gnomad_sv.bed"
#  shell:  """
#    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
#    -b stdin -c 5,5,6,6,30,35,39,42,52,62,72 -o collapse,count,max,mean,max,max,max,max,max,max,max > {output}
#    """
    
    #CTCF Peaks from ucsc
rule CTCF:
  input:
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="../annotations/CTCF/hg38_CTCF.bed.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_ctcf.bed"
  shell: """
    bedtools coverage -b stdin -a {input.bed} > {output}
    """

#gc content from ucsc
rule GC:
  input:
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="../annotations/GC/gc5Base.bedGraph.gz"
  conda: "envs/rules2.yml"
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
  conda: "envs/rules2.yml"
  output: 
      t1=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/test.bed")
  shell: """
    tabix {input.anno} -R {input.merg} | awk '{{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \\"\\";"; }}' | gtf2bed - | cut -f1,2,3,8,10  > {output.t1}
    """
    # genemodels 
rule gene_model:
  input:
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/test.bed"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genemodel.bed"
  shell: """
    grep "exon" {input.anno} | bedtools coverage -b stdin -a {input.bed} |cut -f 1,2,3,5 |paste - <(grep "transcript" {input.anno} |bedtools coverage -b stdin -a {input.bed} |cut -f 5 ) |paste - <(grep "gene" {input.anno} | bedtools coverage -b stdin -a {input.bed} | cut -f 5)|paste - <(grep "start_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "stop_codon" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5) |paste - <(grep "three_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "five_prime_utr" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 )|paste - <(grep "CDS" {input.anno} |bedtools coverage -b stdin -a {input.bed} | cut -f 5 ) >> {output}   """
  
  # gene names necessary for PlI extraction
rule gene_names:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/test.bed"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genenames.bed"
  shell: """
    cut -f 1,2,3 {input.anno} | paste - <(cut -f 5 {input.anno} | tr ";" "\n" | grep -E 'gene_name' | cut -f 1 | tr " " "\t" | cut -f 3 ) | bedtools map -b stdin -a {input.bed} -c 4 -o distinct >{output}
    """

rule pli:
  input:  
    gn="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genenames.bed",
    pli="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/gnomad/pli_exac.csv"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_pli.bed"
  shell:  """
    Rscript --vanilla ../scripts/PLIextract.R {input.gn} {input.pli} {output}
    """

rule cadd1:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/CADD/whole_genome_SNVs.tsv.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_cadd_meanmax.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | awk '{{print $1 "\\t" $2 "\\t" ($2 +1) "\\t" $3 "\\t" $4 "\\t" $5 "\\t" $6 }}' | bedtools map -a {input.bed} -b stdin -c 6,7,6,7 -o max,max,mean,mean > {output}
    """
    #enhancer promotor links
rule EP:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    ep="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/enhancer-promoter-links/sorted_encode.bed"
  conda: "envs/rules2.yml"
  output: 
    o1="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_EP.bed",
    o2=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/EP_tmp.bed")
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.ep} > {output.o2}
    Rscript --vanilla ../scripts/annotateHIC.R {output.o2} {output.o1}
    """

#frequently interacting regulatory elements
rule fire:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{sets}_nochr.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/FIRE/fire_{celllines}.bed"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_fire_{celllines}.bed"
  shell:  """
    bedtools map -a {input.bed} -b {input.anno} -c 4,4 -o max,min > {output}
    """
  
def mergeFire(wc):
  return(expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_fire_{fire}.bed",set=wc.set,fire=CL))

rule fire2:
  input: mergeFire
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_fire.bed"
  conda: "envs/rules.yml"
  shell: "paste {input} > {output}"
  
  
  #hic data from encode 
rule HIC_hESC:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    hic="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/hic/hESC/combined/sorted.total.combined.domain"
  conda: "envs/rules2.yml"
  output: 
    o1="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_HIC_hESC.bed",
    o2=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/HIC_tmp.bed")
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} > {output.o2}
    Rscript --vanilla ../annotateHIC.R {output.o2} {output.o1}
     """
    
rule HIC_encode:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{sets}.bed",
    hic="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/Encode-HIC/{cells}/sorted_nested.bed.gz"
  conda: "envs/rules2.yml"
  output: 
    o1="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_{cells}nested.bed",
    o2=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{cells}_nested_tmp.bed")
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} > {output.o2}
    Rscript --vanilla ../annotateHIC.R {output.o2} {output.o1}
    """

rule HIC_tad:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{sets}.bed",
    hic="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/Encode-HIC/{cells}/sorted_tad.bed.gz"
  conda: "envs/rules2.yml"
  output: 
    o1="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_{cells}tad.bed",
    o2=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{cells}tad_tmp.bed")
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} > {output.o2}
    Rscript --vanilla ../annotateHIC.R {output.o2} {output.o1}
    """

def mergeHic(wc):
  return(expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_{cells}{tad}.bed",set=wc.set,cells=CELLS,tad=TAD))

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
  conda: "envs/rules2.yml"
  output: 
    o1="/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/{sets}_microsynteny.bed",
    o2=temp("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{sets}/microsyn_tmp.bed")
  shell:  """
    bedtools closest -d -t first -a {input.bed} -b {input.hic} > {output.o2}
    Rscript --vanilla ../annotateHIC.R {output.o2} {output.o1}
    """


rule phastcons20:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/PhastCons/hg38.phastCons20way.bedGraph.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_PhastCons20_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg}  | bedtools map -a {input.bed} -b stdin -c 4,4 -o max,mean > {output}"
  
rule phastcons30:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/PhastCons/hg38.phastCons30way.bedGraph.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_PhastCons30_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg}  | bedtools map -a {input.bed} -b stdin -c 4,4 -o max,mean > {output}"

rule phastcons100:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/PhastCons/hg38.phastCons100way.bedGraph.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_PhastCons100_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg}  | bedtools map -a {input.bed} -b stdin -c 4,4 -o max,mean > {output}"
  
rule phyloP:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/PhastCons/hg38.phyloP20way.bedGraph.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_phyloP20_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg}  | bedtools map -a {input.bed} -b stdin -c 4,4 -o max,mean > {output}"
   
   # directionality index of HIC data from genomegitar database  
rule genomegitar1:
  input:  
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
    anno="/fast/groups/ag_kircher/work/SV-scoring/SV/annotations/genomegitar/{gg}/DI_sort.bed"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genomegitar_{gg}.bed"
  shell:  """
    bedtools map -a {input.bed} \
    -b {input.anno} -c 4,4 -o max,min > {output}
    """
    
def mergeGgitars(wc):
  return(expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_genomegitar_{gg}.bed",set=wc.set,gg=genomegitars))

rule genomegitar2:
  input: mergeGgitars
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_DI.bed"
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


#rule mirTarget:
#  input:  
#    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_merged.bed",
#    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}.bed",
#    anno="/fast/groups/ag_kircher/CADD/dependencies/annotations/mirTargetScan/targetScan_hg38_v71.bg.gz"
#  conda: "envs/rules2.yml"
#  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_mirTarget_mean.bed"
#  shell:  """
#    tabix {input.anno} -R {input.merg}  | bedtools map -a {input.bed} \
#    -b stdin -c 4 -o mean > {output}
#   """

rule MPC:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="/fast/groups/ag_kircher/CADD/dependencies/annotations/MPC/transcript_constraints_hg38liftover.bg.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_MPC_mean.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o mean > {output}
    """
 
rule RemapTF:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="/fast/groups/ag_kircher/CADD/dependencies/annotations/Remap/ReMap2_overlapTF_hg38.bg.gz"
  conda: "envs/rules2.yml"
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
    anno="/fast/groups/ag_kircher/CADD/dependencies/annotations/encode/{encodes}/{encodes}_merged.bg.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_encode_{encodes}_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} -b stdin -c 4 -o mean > {output}"

def mergeEncode(wc):
  return(expand("/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_encode_{encodes}_mean.bed",set=wc.set,encodes=ENCODES))

rule encode2:
  input: mergeEncode
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_encode.bed"
  conda: "envs/rules.yml"
  shell: "paste {input} > {output}"


rule gerp:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="/fast/groups/ag_kircher/CADD/dependencies/annotations/gerp/gerp_score2_hg38_MAM.bg.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_gerp_mean.bed"
  shell:  "tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} -b stdin -c 4 -o mean > {output}"


rule ccr:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="/fast/groups/ag_kircher/CADD/dependencies/annotations/ccr/ccrs.all.bed.gz"
  conda: "envs/rules2.yml"
  output: "/fast/groups/ag_kircher/work/SV-scoring/SV/run/{set}/{set}_ccr_mean.bed"
  shell:  """
    tabix {input.anno} -R {input.merg} | bedtools map -a {input.bed} \
    -b stdin -c 4 -o mean > {output}
    """
    
rule chromHMM_MAX:
  input:  
    merg="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr_merged.bed",
    bed="/fast/groups/ag_kircher/work/SV-scoring/SV/beds/{set}_nochr.bed",
    anno="/fast/groups/ag_kircher/CADD/dependencies/annotations/chromhmm/chromHMM_GRCh38.bg.gz"
  conda: "envs/rules2.yml"
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
  conda: "envs/rules2.yml"
  shell:  """
    paste <(cut -f1,2,3,4,5,6,7 {input.cadd}) <(cut -f4 {input.ccr}) <(cut -f4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 {input.chromHMM}) <(cut -f4,5,6,7 {input.ctcf}) <(cut -f1 {input.di_min}) <(cut -f1 {input.di_max}) <(cut -f4,8,12,16,20,24,28,32,36,40,44,48,52 {input.encode}) <(cut -f4,5,6,7 {input.ep}) <(cut -f4,5,9,10,14,15,19,20,24,25 {input.fire}) <(cut -f4 {input.gc}) <(cut -f4,5,6,7,8,9,10,11 {input.gm}) <(cut -f4 {input.gerp}) <(cut -f4,5,6,7,11,12,13,14,18,19,20,21,25,26,27,28 {input.hic}) <(cut -f4,5,6,7 {input.hesc}) <(cut -f4,5,6,7 {input.microsyn}) <(cut -f4 {input.mpc}) <(cut -f4,5 {input.pc1}) <(cut -f4,5 {input.pc2}) <(cut -f4,5 {input.pc3}) <(cut -f4,5 {input.phyloP}) <(cut -f4 {input.pli})  <(cut -f4 {input.remapTF}) > {output}
    """
  
  # chr start end cadd1_max cadd2_max cadd1_mean cadd2_mean ccr cHMM_1 cHMM_2 cHMM_3 cHMM_4 cHMM_5 cHMM_6 cHMM_7 cHMM_8 cHMM_9 cHMM_10 cHMM_11 cHMM_12 cHMM_13 cHMM_14 cHMM_15 cHMM_16 cHMM_17 cHMM_18 cHMM_19 cHMM_20 cHMM_21 cHMM_22 nr_ctcf_BS ctcf_nr_bases ctcf_length ctcf_fraction DI_min DI_max DNase-seq H2AFZ H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K79me2 H3K9ac H3K9me3 H4K20me1 totalRNA-seq EP_intra EP_boundary EP_ext EP_distance gm12878_max gm12878_min msc_max msc_min mes_max mes_min imr90_max imr90_min h1_max h1_min perc_gc exon transcript gene start_codon stop_codon 3utr 5utr cds gerp A549_nested_intra A549_nested_boundary A549_nested_ext A549_nested_dist A549_tad_intra A549_tad_boundary A549_tad_ext A549_tad_dist caki2_nested_intra caki2_nested_boundary caki2_nested_ext caki2_nested_dist caki2_tad_intra caki2_tad_boundary caki2_tad_ext caki2_tad_dist escTAD_intra escTAD_boundary escTAD_ext escTAD_distance microsyn_intra miscrosyn_boundary microsyn_ext microsyn_distance MPC pc20_max pc20_mean pc30_max pc30_mean pc100_max pc100_mean phyoP_max phyloP_mean pli remapTF
    
    
    
    
