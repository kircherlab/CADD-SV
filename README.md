# CADD-SV 

## CADD-SV – a framework to score the effect of structural variants 

Here we describe CADD-SV, a method to retrieve a wide set of annotations in the range and vicinity of a SV. Our tool computes summary statistics and uses a trained linear model to differentiate deleterious from neutral variants. We use human and chimpanzee derived alleles as proxy-neutral and contrast them with matched simulated variants as proxy-pathogenic, an approach that has proven powerful in the interpretation of SNVs (CADD). We show that CADD-SV-scores correlate with known pathogenic variants in individual genomes and allelic diversity. The ability of this model to prioritize functionally relevant, pathogenic variants is unmatched by existing methods.
philip.kleinert@bihealth.de


## Pre-requirements

### Conda

The pipeline depends on [Snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system that wraps up all scripts and runs them highly automated, in various environments (workstations, clusters, grid, or cloud). Further, we use Conda as software/dependency managment tool. Conda can install snakemake and all neccessary software with its dependencies automatically. Conda installation guidlines can be found here:

https://conda.io/projects/conda/en/latest/user-guide/install/index.html

### Snakemake

After installing Conda, you install Snakemake using Conda and the `environment.yaml` provided in this repository. For this purpose, please clone or download and uncompress the repository first. Then change into the root folder of the local repository. 

```bash
git clone https://github.com/kircherlab/CADD-SV
cd CADD-SV
```

We will now initiate the Conda environment, which we will need for getting the Snakemake workflow invoked. Using this environment (`run.caddsv`) snakemake will be installed

```bash
conda env create -n run.caddsv --file environment.yaml
```

The conda environment (`envs/SV.yml`) containing all packages and tools to make CADD-SV run will be installed automatically during the first run only. This can take some time.

### Annotations

CADD-SV depends on various annotations to provide the model with its necessary input features. CADD-SV automatically retrieves and transforms these annotations (see Snakefile) and combines them in bed-format at `/desired-sv-set/matrix.bed`

Annotations can be downloaded and expanded individually. However, to run CADD-SV as desired and to minimize runtime and memory failures use the annotation sets as stored at https://kircherlab.bihealth.org/download/CADD-SV/

```bash
wget https://kircherlab.bihealth.org/download/CADD-SV/v1.0/dependencies.tar.gz
tar -xf dependencies.tar.gz
```

These annotations include:

##### Integrated Scores
CADD (https://krishna.gs.washington.edu/download/CADD/bigWig/) \
LINSIGHT (compgen.cshl.edu/LINSIGHT/LINSIGHT.bw) 

##### Species conservation and constraint metrics
PhastCons (http://hgdownload.cse.ucsc.edu/goldenpath/hg38/) \
Syntenic regions (http://webclu.bio.wzw.tum.de/cgi-bin/syntenymapper/get-species-list.py) \
gerp score (http://mendel.stanford.edu/SidowLab/downloads/gerp/) 

##### Population and disease constraint metrics
Pli score (ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/manuscript_data/forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz) \
conserved coding regions (CCR) (https://www.nature.com/articles/s41588-018-0294-6?WT.feed_name=subjects_population-genetics) \
DDD Happloinsufficiency (https://decipher.sanger.ac.uk/files/downloads/HI_Predictions_Version3.bed.gz) 

##### Epigenetic and regulatory activity
Encode Features such as Histon Modifications and DNase and RNase-seq (https://www.encodeproject.org/help/batch-download/) \
GC content (http://hgdownload.cse.ucsc.edu/gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw) \
chromhmm states of encode celllines (http://compbio.mit.edu/ChromHMM/) 

##### 3D genome organization
CTCF (http://genome.cshlp.org/content/suppl/2012/08/28/22.9.1680.DC1/Table_S2_Location_of_ChIP-seq_binding_positions_in_19_cell_lines.txt) \
Encode-HIC data (https://www.encodeproject.org/search/?type=Experiment&assay_term_name=HiC&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens&status=released) \
enhancer-promoter-links from FOCS (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5930446/) \
FIREs frequently interacting regulatory elements (https://www.sciencedirect.com/science/article/pii/S2211124716314814) \
Directionality index for HIC data from various datasets (https://www.genomegitar.org/processed-data.html) \
deepc saliencies score (http://userweb.molbiol.ox.ac.uk/public/rschwess/container_for_ucsc/data/deepC/saliency_scores/saliencies_merged_gm12878_5kb.bw) 

##### Gene and element annotations
ensembl-gff3 genebuild (ftp://ftp.ensembl.org/pub/release-96/gff3/homo_sapiens/Homo_sapiens.GRCh38.96.chr.gff3.gz) \
Fantom5 enhancers (https://zenodo.org/record/556775#.Xkz3G0oo-70) 


## Config

Almost ready to go. After you prepared the files above, you may need to adjust locations and names of these files in the `config.yml`. 

## List of required input files

- CADD-SV scores SV in bedformat on the GRCh38 genomebuild. The type of SV needs tobe contained in the 4th column. We recomment to split files containing more than 10.000 SVs into smaller files.

  If you want to score SVs in a VCF format or your SVs are not in GRCh38 genomebuild coordinates:
  We provide a environment to handle this.
  
  ```bash
  conda env create -n prepBED --file envs/prepBED.yml
  ```
  
  To convert your VCF into BED format run:
  ```
  conda activate prepBED
  SURVIVOR vcftobed input.vcf 0 -1 output.bed
  cut -f1,2,6,11 output.bed > beds/set_id.bed
  
  ```
  
  To lift hg19 coordinates to GRCh38 apply following steps:
  
  ```
  conda activate prepBED
  liftOver beds/setname_hg19_id.bed /dependencies/hg19ToHg38.over.chain.gz beds/setname_id.bed beds/setname_unlifted.bed
  ```
   
- Annotations in the /dependencies folder
- Models and scripts as cloned from this GIT repository


## Run pipeline

Ready to go! If you run the pipeline on a cluster see the `cluster.json` for an estimate of minimum requirements for the individual jobs. Note that this depends on your dataset size so you may have to adjust this.
To start the pipeline:

```bash
conda activate run.caddsv
# dry run to see if everything works
snakemake  --use-conda --configfile config.yml -n
# run the pipeline
snakemake  --use-conda --configfile config.yml
```

## Output files

The pipeline outputs your SV input containing all annotations in bed format as /SV-to-score/matrix.bed
Further a file only containing the coordinates and the cadd-sv score is outputed at /SV-to-score/score.bed


