################################################################################################
#Parameters #######################################################################################
dir_script=         "/net/nfs-irwrsrchnas01/labs/zgu_grp/DataFolder/Seq/mdall/scripts/"

dir_ref=            "/net/nfs-irwrsrchnas01/labs/zgu_grp/DataFolder/Seq/mdall/ref/"
ref_fa=             dir_ref + "genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
gtf=                dir_ref + "genome_anno/Homo_sapiens.GRCh38.V102.withChr.gtf"
bed_DUX4=           dir_ref + "DUX4patch/Homo_sapiens.GRCh38.V102.withChr.DUX4.bed"
ref_star=           dir_ref + "index_star"
ref_fusioncatcher=  dir_ref + "fusioncatcher/human_v102"
ref_RNApeg_fa=      dir_ref + "genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
ref_RNApeg_flat=    dir_ref + "genome_anno/refFlat.txt"
ref_cicero=         dir_ref + "cicero"

set_chr=            ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y']]

cores_star=             8
cores_samtoolsSort=     8
cores_fusioncatcher=    8
cores_RNApeg=           8
cores_cicero=           8

################################################################################################
#Project #######################################################################################
dir_in=         "/net/nfs-irwrsrchnas01/labs/zgu_grp/DataFolder/Seq/mdall/test"
dir_out=        "/net/nfs-irwrsrchnas01/labs/zgu_grp/DataFolder/Seq/mdall/out_raw"

# samplelist=['COH000892_10M_D1']
# samplelist=['COH000892_5M_D1']
# samplelist=['COH000892_1M_D1']

# samplelist=['COH000892_1M_D1','COH000892_5M_D1']

# samplelist=['COH004898_D1.5M']
samplelist=['COH004922_D1.5M']

################################################################################################
#Pipeline ######################################################################################
rule all:
    input:
        expand(dir_out + "/{sample}/TRANSCRIPTOME/log/Done.txt",sample=samplelist)
        # expand(dir_out + "/{sample}/TRANSCRIPTOME/log/cicero.bmk",sample=samplelist)
# 
include: "scripts/rnaseq.smk"
include: "scripts/cicero.smk"

