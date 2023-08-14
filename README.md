
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Analysis pipeline for RNA-seq data

### The workflow of this pipeline:

<br>

<img src="pics/RNAseq_diagflow.JPG" align="center" width="60%" height="60%"/>
<br> <br>

### Run the pipeline

To run the pipeline with 16 cores:

``` bash
snakemake -s rnaseq.smk -p -j16
```
