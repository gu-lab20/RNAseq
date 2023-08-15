
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

## Dependencies

#### Please install all the needed softwares before running this pipeline:

[perl-v5.26.2](https://www.perl.org/)

[Snakemake-7.16.0](https://snakemake.readthedocs.io/en/stable/)

[STAR-2.7.2b](https://github.com/alexdobin/STAR)

[samtools-1.12](http://www.htslib.org/)

[GATK-4.3.0.0](https://gatk.broadinstitute.org/hc/en-us)

[HTSeq-count-0.13.5](https://htseq.readthedocs.io/en/release_0.11.1/count.html)

[FusionCatcher-1.33](https://github.com/ndaniel/fusioncatcher)

[RNApeg](https://github.com/stjude/RNApeg)

[Cicero-0.3.0p2](https://github.com/stjude/CICERO)

[RNAseqCNV-1.2.2](https://github.com/honzee/RNAseqCNV)

## Configuration

Users can edit the run\_rnaseq.smk file for configurations.

Parameters:

‘<strong>ref\_fa</strong>’, the fasta file of reference genome of human
GRCh38. Users need to download it.

‘<strong>gtf</strong>’, gtf annotation file of the reference genome.
Users need to download it.

‘<strong>bed\_DUX4</strong>’, bed file of DUX4 genes. This file is used
in the read counts patching process for DUX4 genes. Already included in
the 0.ref directory.

‘<strong>ref\_star</strong>’, the directory of reference used by STAR to
do alignment. Users will get it after the installation of STAR.

‘<strong>ref\_fusioncatcher</strong>’, the directory of reference used
by FusionCatcher to call gene fusions. Users will get it after the
installation of FusionCatcher.

‘<strong>ref\_cicero</strong>’, the directory of reference used by
Cicero to call gene fusions. Users will get it after the installation of
Cicero.

‘<strong>ref\_RNApeg\_flat</strong>’, the refFlat file used by RNApeg.
Already included in the 0.ref directory.

‘<strong>cores\_star</strong>’, ‘<strong>cores\_samtools</strong>’,
‘<strong>cores\_fusioncatcher</strong>’,
‘<strong>cores\_RNApeg</strong>’ and ‘<strong>cores\_cicero</strong>’
are the number of threads used by the the corresponding software.

‘<strong>dir\_in</strong>’, the directory of input fastq files. Only gz
compressed paired-end fastq files are supported currently. The file
names should follow the pattern {sample}.R1.fq.gz and {sample}.R2.fq.gz.
If a sample id is COH000456\_D1, then the fastq file names should be
COH000456\_D1.R1.fq.gz and COH000456\_D1.R2.fq.gz.

‘<strong>dir\_out</strong>’, the output directory. Results will be
stored in sub-directories within this folder, each named according to
the respective sample ID.

‘<strong>samplelist</strong>’, the sample ID list that will be processed
for analysis. The corresponding fastq files need to be stored in the
directory ‘<strong>dir\_in</strong>’.

## Run the pipeline

After installation of all the required softwares and propoably
configured the parameters, a dry-run could be executed with:

``` bash
snakemake -s run_rnaseq.smk -pn
```

To run the pipeline with 16 cores:

``` bash
snakemake -s run_rnaseq.smk -p -j16
```
