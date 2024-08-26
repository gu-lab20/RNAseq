################################################################################################
#Parameters #######################################################################################
dir_script=         "/Full_path_to/mdall/scripts/"
dir_ref=            "/Full_path_to/ref/"

ref_fa=             dir_ref + "genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
gtf=                dir_ref + "genome_anno/Homo_sapiens.GRCh38.V102.withChr.gtf"
bed_DUX4=           dir_ref + "DUX4patch/Homo_sapiens.GRCh38.V102.withChr.DUX4.bed"
ref_star=           dir_ref + "index_star"
ref_fusioncatcher=  dir_ref + "fusioncatcher/human_v102"
ref_RNApeg_fa=      dir_ref + "genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
ref_RNApeg_flat=    dir_ref + "genome_anno/refFlat.txt"

set_chr=            ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y']]

cores_star=             8
cores_samtoolsSort=     8
cores_fusioncatcher=    8
cores_RNApeg=           8
cores_cicero=           8

################################################################################################
#Project #######################################################################################

dir_in=config["dir_in"]
dir_out=config["dir_out"]
samplelist=config["sample"]

################################################################################################
#Pipeline ######################################################################################
rule all:
    input:
        expand(dir_out + "/{sample}/TRANSCRIPTOME/log/Done.txt",sample=samplelist)

include: dir_script+"rnaseq.smk"
