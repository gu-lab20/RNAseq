################################################################################################
#Pipeline ######################################################################################
ref_fa=             "0.ref/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
gtf=                "0.ref/genome_anno/Homo_sapiens.GRCh38.V102.withChr.gtf"
bed_DUX4=           "0.ref/DUX4patch/Homo_sapiens.GRCh38.V102.withChr.DUX4.bed"
ref_star=           "0.ref/index_star"
ref_fusioncatcher=  "0.ref/fusioncatcher/human_v102"
ref_cicero=         "/ref_genomes/CICERO/human/GRCh38/reference"
ref_RNApeg_fa=      ref_fa
ref_RNApeg_flat=    "0.ref/genome_anno/refFlat.txt"

set_chr=            ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y']]

cores_star=             8
cores_samtoolsSort=     8
cores_fusioncatcher=    8
cores_RNApeg=           8
cores_cicero=           8

################################################################################################
#Project #######################################################################################
dir_in=         "0.original"
dir_out=        "out_raw"

samplelist=     ['COH000456_D1','COH000893_D1']

################################################################################################
#Pipeline ######################################################################################
rule all:
    input:
        expand(dir_out + "/{sample}/TRANSCRIPTOME/log/Done.txt",sample=samplelist)

include: "rnaseq.smk"
