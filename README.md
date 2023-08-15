
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

Users can edit the run\_rnaseq.smk file for configurations.

Parameters:

<strong>ref\_fa</strong>, the fasta file of reference genome of human
GRCh38. Users need to download it.

<strong>gtf</strong>, gtf annotation file of the reference genome. Users
need to download it.

<strong>bed\_DUX4</strong>, bed file of DUX4 genes. This file is used in
the read counts patching process for DUX4 genes. Already included in the
0.ref directory.

<strong>ref\_star</strong>, the directory of reference used by STAR to
do alignment. Users will get it after the installation of STAR.

<strong>ref\_fusioncatcher</strong>, the directory of reference used by
FusionCatcher to call gene fusions Users will get it after the
installation of FusionCatcher.

<strong>ref\_cicero</strong>, the directory of reference used by Cicero
to call gene fusions Users will get it after the installation of Cicero.

<strong>ref\_RNApeg\_flat</strong>, the refFlat file used by RNApeg.
Already included in the 0.ref directory.

<strong>cores\_star</strong>, <strong>cores\_samtools</strong>,
<strong>cores\_fusioncatcher</strong>, <strong>cores\_RNApeg</strong>
and <strong>cores\_cicero</strong> are the number of threads used by the
name of software indicated.

<strong>dir\_in</strong>,

<strong>dir\_out</strong>,

<strong>samplelist</strong>,

### Run the pipeline

To run the pipeline with 16 cores:

``` bash
snakemake -s rnaseq.smk -p -j16
```
