
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Analysis pipeline for RNA-seq data

This pipeline is specifically designed for the analysis of RNA-seq data.
It’s written in
[Snakemake](https://snakemake.readthedocs.io/en/stable/). Upon
execution, the pipeline will produce outputs detailing gene read counts,
mutations, fusions, and chromosomal level copy number variations
(gains/losses) derived from the RNA-seq data.

The workflow of this pipeline:
<br>

<img src="pics/RNAseq_diagflow.JPG" align="center" width="60%" height="60%"/>
<br>
<br>

### Dependencies

#### Please install all the needed softwares before running this pipeline:

perl

[Snakemake](https://snakemake.readthedocs.io/en/stable/)

[STAR](https://github.com/alexdobin/STAR)

[samtools](http://www.htslib.org/)

[GATK](https://gatk.broadinstitute.org/hc/en-us)

[HTSeq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html)

[FusionCatcher](https://github.com/ndaniel/fusioncatcher)

[RNApeg](https://github.com/stjude/RNApeg)

[Cicero](https://github.com/stjude/CICERO)

[RNAseqCNV](https://github.com/honzee/RNAseqCNV)

### Configuration

Users can edit the run\_rnaseq.smk file for configurations. Two
parameter panels included: pipeline and project. Parameters in pipeline:
<strong>ref\_fa</strong>, <strong>gtf</strong>,
“0.ref/genome\_anno/Homo\_sapiens.GRCh38.V102.withChr.gtf”
<strong>bed\_DUX4</strong>,
“0.ref/DUX4patch/Homo\_sapiens.GRCh38.V102.withChr.DUX4.bed”
<strong>ref\_star</strong>, “0.ref/index\_star”
<strong>ref\_fusioncatcher</strong>, “0.ref/fusioncatcher/human\_v102”
<strong>ref\_cicero</strong>,
“/ref\_genomes/CICERO/human/GRCh38/reference”
<strong>ref\_RNApeg\_fa</strong> ref\_fa
<strong>ref\_RNApeg\_flat</strong>, “0.ref/genome/refFlat.txt”

<strong>set\_chr<strong>, \[‘chr{}’.format(x) for x in list(range(1,23))
+ \[‘X’, ‘Y’\]\]

<strong>cores\_star</strong>, <strong>cores\_samtools</strong>,
<strong>cores\_fusioncatcher</strong>,<strong>cores\_RNApeg</strong>,<strong>cores\_cicero</strong>,

### Run the pipeline

To run the pipeline with 16 cores:

``` bash
snakemake -s rnaseq.smk -p -j16
```
