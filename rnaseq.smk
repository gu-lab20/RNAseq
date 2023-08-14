import pandas as pd
dir_project = "/home/zgu_labs/pipeline_online/rnaseq/"

ref_fa=dir_project+"0.ref/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
gtf="0.ref/genome_anno/Homo_sapiens.GRCh38.V102.withChr.gtf"
bed_DUX4="0.ref/DUX4patch/Homo_sapiens.GRCh38.V102.withChr.DUX4.bed"
ref_star="0.ref/index_star"

ref_fusioncatcher="0.ref/fusioncatcher/human_v102"

# ref_fusioncatcher="/ref_genomes/fusioncatcher/human_v102_new/human_v102"

ref_cicero="/ref_genomes/CICERO/human/GRCh38/reference"

# ref_RNApeg_dir="/ref_genomes/CICERO/human/GRCh38/reference/Homo_sapiens/GRCh38_no_alt/mRNA/RefSeq,/ref_genomes/CICERO/human/GRCh38/reference/Homo_sapiens/GRCh38_no_alt/FASTA"
ref_RNApeg_dir=dir_project+"0.ref/genome,/ref_genomes/CICERO/human/GRCh38/reference/Homo_sapiens/GRCh38_no_alt/FASTA"

# ref_RNApeg_fa="/home/zgu_labs/pipeline_online/rnaseq/0.ref/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
ref_RNApeg_fa="/ref_genomes/CICERO/human/GRCh38/reference/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa"
ref_RNApeg_flat=dir_project+"0.ref/genome/refFlat.txt"



# set_chr=['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y']]

set_chr="chr1"

dir_out=    dir_project + "out_raw"; 
dir_in=     dir_project + "0.original"; 
id_info =   dir_project + "0.original/df_id.tsv"

samplelist=['COH000456_D1','COH000893_D1']

#output --##--##--##--##--
# key_file=dir_out + "/{sample}/TRANSCRIPTOME/log/DoneMapping.txt"
# key_file=dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.RNAseqCNV.bmk"
# key_file=dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.Cicero.bmk"
# key_file=dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.FusionCatcher.bmk"
key_file=dir_out + "/{sample}/TRANSCRIPTOME/log/Done.txt"

rule all:
    input:
        # key_file
        # expand(dir_code +"/{sample}.sh",sample=samplelist)
        # expand(key_file,sample=samplelist)
        expand(key_file,sample=samplelist[0]),
        # expand(key_file2,sample=samplelist)
        # expand(key_file,id="test1")

rule Star:
    input:
        fq1=        dir_in + '/{sample}.R1.fq.gz',
        fq2=        dir_in + '/{sample}.R2.fq.gz',
    output:
        bam=        temp(dir_out + "/{sample}/TRANSCRIPTOME/star/{sample}.bam"),
        junction=   dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.star.junction",
    benchmark:      dir_out + "/{sample}/TRANSCRIPTOME/log/star/{sample}.star.bmk",
    log:            dir_out + "/{sample}/TRANSCRIPTOME/star/{sample}.star.log",
    params:
        star_ref=       ref_star,
        rg=             "ID:{sample} PL:ILLUMINA PU:ILLUMINA LB:RNAseq SM:{sample}",
        dir_out=        dir_out + "/{sample}/TRANSCRIPTOME/star/",
        bam_tmp=        dir_out + "/{sample}/TRANSCRIPTOME/star/Aligned.out.bam",
        junction_tmp=   dir_out + "/{sample}/TRANSCRIPTOME/star/Chimeric.out.junction"
    threads: 8
    shell:
        '''
        STAR \
        --twopassMode Basic \
        --runThreadN {threads} \
        --readFilesCommand zcat \
        --limitOutSAMoneReadBytes 90000000 \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        --outSAMmapqUnique 60 \
        --outSAMmultNmax 1 \
        --outReadsUnmapped None \
        --outMultimapperOrder Random \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMattrRGline {params.rg} \
        --outFilterType BySJout \
        --outFilterMultimapNmax 100 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 4 \
        --alignInsertionFlush Right \
        --chimOutType Junctions \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 100 \
        --chimMultimapScoreRange 10 \
        --chimNonchimScoreDropMin 10 \
        --chimOutJunctionFormat 1 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --outFileNamePrefix {params.dir_out} \
        --genomeDir {params.star_ref} \
        --genomeLoad NoSharedMemory \
        --readFilesIn {input.fq1} {input.fq2} \
        --sjdbOverhang 100 &> {log}

        mv {params.bam_tmp} {output.bam}
        mv {params.junction_tmp} {output.junction}
        '''

rule samtools_sort:
    input:
        bam=        dir_out + "/{sample}/TRANSCRIPTOME/star/{sample}.bam",
        bmk=        dir_out + "/{sample}/TRANSCRIPTOME/log/star/{sample}.star.bmk",
    output:
        bam=        temp(dir_out + "/{sample}/TRANSCRIPTOME/samtools_sort/{sample}.bam"),
        bai=        temp(dir_out + "/{sample}/TRANSCRIPTOME/samtools_sort/{sample}.bam.bai"),
        dir=        temp(directory(dir_out + "/{sample}/TRANSCRIPTOME/samtools_sort/")),
    log:            dir_out + "/{sample}/TRANSCRIPTOME/log/samtools_sort.log",
    benchmark:      dir_out + "/{sample}/TRANSCRIPTOME/log/samtools_sort.bmk",
    params:
        dir_star=   dir_out + "/{sample}/TRANSCRIPTOME/star/",
    threads: 8
    shell:
        '''
        samtools sort -m 1G -@ {threads} -O bam -T {output.dir} -o {output.bam} {input.bam}
        samtools index {output.bam} {output.bai} &>> {log}
        if [ -d {params.dir_star} ]; then rm -r {params.dir_star}; fi
        '''

rule MarkDup:
    input:
        dir=        dir_out + "/{sample}/TRANSCRIPTOME/samtools_sort/",
        bam=        dir_out + "/{sample}/TRANSCRIPTOME/samtools_sort/{sample}.bam",
        bmk=        dir_out + "/{sample}/TRANSCRIPTOME/log/samtools_sort.bmk"
    output:
        bam=        dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        bai=        dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam.bai",
        metrics=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.MarkDup.metrics.txt",
    log:            dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.markDup.log"
    benchmark:      dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.markDup.bmk"
    params:
        bam_star=   dir_out + "/{sample}/TRANSCRIPTOME/star/{sample}.bam",
    shell:
        '''        
        gatk MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} &> {log}
        samtools index {output.bam} {output.bai} &>> {log}
        '''

rule HTSeq_count:
    input:
        bam=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        bmk=    dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.markDup.bmk",
    output:
        HTSeq_count=    dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.HTSeq",
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.HTSeq_count.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.HTSeq_count.bmk"
    params:
        gtf=    gtf
    shell:
        '''
        htseq-count -m intersection-strict -f bam -r pos -s no -a 10 {input.bam} {params.gtf} 1> {output.HTSeq_count} 2> {log}
        '''

rule HTSeqDUX4:
    input:
        bam=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        htseq=  dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.HTSeq",
        bmk=    dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.HTSeq_count.bmk"
    output:
        HTSeq=  dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.DUX4patched.HTSeq",
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.HTSeqDUX4.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.HTSeqDUX4.bmk"
    params:
        bed=    bed_DUX4
    shell:
        '''
        0.scripts/HTSeqDUX4patch.pl -f {params.bed} -b {input.bam} -h {input.htseq} -o {output.HTSeq}
        '''

rule FusionCatcher:
    input:
        fq= lambda wildcards: expand(dir_in + '/{sample}.R{i}.fq.gz',i=[1,2],sample=wildcards.sample)
    output:
        dir=    directory(dir_out + "/{sample}/TRANSCRIPTOME/FusionCatcher")
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.FusionCatcher.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.FusionCatcher.bmk"
    params:
        ref=        ref_fusioncatcher,
        dir_out=    dir_out + "/{sample}/TRANSCRIPTOME/FusionCatcher",
        fq1=        dir_in + '/{sample}.R1.fq.gz',
        fq2=        dir_in + '/{sample}.R2.fq.gz',
    threads: 8
    shell:
        '''        
        fusioncatcher.py -p {threads} -d {params.ref} -i {params.fq1},{params.fq2} -o {output.dir} &>{log}
        '''

rule RNApeg:
    input:
        bam = dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        bai = dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam.bai",
    output:
        RNApeg=     dir_out + "/{sample}/TRANSCRIPTOME/RNApeg/done",
        junctions=  dir_out + "/{sample}/TRANSCRIPTOME/RNApeg/{sample}.bam.junctions.tab.shifted.tab"
    log:            dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.RNApeg.log"
    benchmark:      dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.RNApeg.bmk"
    params:
        ref1=   ref_RNApeg_fa,
        ref2=   ref_RNApeg_flat,
        dir_out=dir_out + "/{sample}/TRANSCRIPTOME/RNApeg",
        dir_bam=dir_out + "/{sample}/TRANSCRIPTOME/bam",
        dir_ref=ref_RNApeg_dir,
    threads: 8
    shell:
        '''
        singularity run --containall \
        --bind {params.dir_ref},{params.dir_bam},{params.dir_out}:/results /packages/singularity-images/rnapeg.simg RNApeg.sh \
        -b {input.bam} -f {params.ref1} -r {params.ref2}  &>{log}
        touch {output.RNApeg}
        '''

rule Cicero:
    input:
        bam=        dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        bai=        dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam.bai",
        junctions=  dir_out + "/{sample}/TRANSCRIPTOME/RNApeg/{sample}.bam.junctions.tab.shifted.tab",
        bmk=        dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.RNApeg.bmk"
    output:
        done=   dir_out + "/{sample}/TRANSCRIPTOME/cicero/done.txt",
        dir=    directory(dir_out + "/{sample}/TRANSCRIPTOME/cicero")
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.Cicero.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.Cicero.bmk"
    params:
        ref=        ref_cicero,
        dir_bam=dir_out + "/{sample}/TRANSCRIPTOME/bam/",
        dir_junctions=dir_out + "/{sample}/TRANSCRIPTOME/RNApeg/",
        dir_cicero= dir_out + "/{sample}/TRANSCRIPTOME/cicero",
        dir_temp=   dir_out + "/{sample}/TRANSCRIPTOME/Cicero",
        dir_fusion= dir_out + "/{sample}/TRANSCRIPTOME/Cicero/CICERO_DATADIR/{sample}",
    threads: 8
    shell:
        '''
        if [ -d {params.dir_temp} ]; then rm -rf {params.dir_temp}; fi
        singularity exec --bind {params.ref},{params.dir_bam},{params.dir_junctions} /packages/singularity-images/cicero_0.3.0p2.sif Cicero.sh \
        -n {threads} -b {input.bam} -g GRCh38_no_alt -r {params.ref} -j {input.junctions} -s 2 -c 10 -o {params.dir_temp} &>{log}
        mv {params.dir_fusion}/*.txt {params.dir_cicero}
        touch {output.done}
        rm -rf {params.dir_temp}
        '''

rule SplitNCigarReads:
    input:
        bam=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        #bmk=   dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}_markDup.bmk",
    output:
        bam_split=  temp(dir_out + "/{sample}/TRANSCRIPTOME/SplitNCigarReads/{sample}.{chr}.split.bam"),
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/SplitNCigarReads/{sample}.{chr}.SplitNCigarReads.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/SplitNCigarReads/{sample}.{chr}.SplitNCigarReads.bmk"
    params:
        ref=ref_fa,
        chr="{chr}",
    shell:
        '''
        gatk SplitNCigarReads -R {params.ref} -I {input.bam} --intervals {params.chr} -O {output.bam_split} &> {log}
        '''

rule HaplotypeCaller:
    input:
        bam=    dir_out + "/{sample}/TRANSCRIPTOME/SplitNCigarReads/{sample}.{chr}.split.bam",
        bmk=    dir_out + "/{sample}/TRANSCRIPTOME/log/SplitNCigarReads/{sample}.{chr}.SplitNCigarReads.bmk"
    output:
        vcf=    dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/{sample}.{chr}.HaplotypeCaller.vcf"
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/HaplotypeCaller/{sample}.{chr}.HaplotypeCaller.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/HaplotypeCaller/{sample}.{chr}.HaplotypeCaller.bmk"
    params:
        ref=ref_fa,
        chr="{chr}",
    shell:
        '''
        gatk HaplotypeCaller --native-pair-hmm-threads {threads} \
        --dont-use-soft-clipped-bases \
        -R {params.ref} --intervals {params.chr} -I {input.bam} \
        -O {output.vcf} &>> {log}
        '''

rule HaplotypeCaller_merge:
    input:
        vcf=lambda wildcards: expand(
            dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/{sample}.{chr}.HaplotypeCaller.vcf",
            chr=set_chr,sample=wildcards.sample),
        bmk=lambda wildcards: expand(
            dir_out + "/{sample}/TRANSCRIPTOME/log/HaplotypeCaller/{sample}.{chr}.HaplotypeCaller.bmk",
            chr=set_chr,sample=wildcards.sample)
    output:
        vcf=    dir_out + "/{sample}/TRANSCRIPTOME/Mutation/{sample}.HaplotypeCaller.vcf",
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.HaplotypeCaller_merge.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.HaplotypeCaller_merge.bmk"
    params:
        filename_prefix=    dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/{sample}.*.HaplotypeCaller.vcf",
        idx_prefix=         dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/{sample}.*.HaplotypeCaller.vcf.idx",
        filename_list=      dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/vcf.list",
        dir_HP=             dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller",
    shell:
        '''
        ls {params.filename_prefix} > {params.filename_list}
        gatk MergeVcfs -I {params.filename_list} -O {output.vcf} &> {log}
        rm -rf {params.dir_HP}
        '''

rule RNAseqCNV:
    input:
        vcf=        dir_out + "/{sample}/TRANSCRIPTOME/Mutation/{sample}.HaplotypeCaller.vcf",
        bmk_vcf=    dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.HaplotypeCaller_merge.bmk",
        counts=     dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.DUX4patched.HTSeq",
        bmk_counts= dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.HTSeqDUX4.bmk",
    output:
        cnv=    directory(dir_out + "/{sample}/TRANSCRIPTOME/RNAseqCNV/{sample}"),
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.RNAseqCNV.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.RNAseqCNV.bmk"
    params:
        config_file=    dir_out + "/{sample}/TRANSCRIPTOME/RNAseqCNV/config",
        outdir=         dir_out + "/{sample}/TRANSCRIPTOME/RNAseqCNV",
        count_dir=      dir_out + "/{sample}/TRANSCRIPTOME/HTSeq",
        snv_dir=        dir_out + "/{sample}/TRANSCRIPTOME/Mutation",
        meta_file=      dir_out + "/{sample}/TRANSCRIPTOME/RNAseqCNV/metadata",
        sample=         "{sample}",
        htseq_file=     "{sample}.DUX4patched.HTSeq",
        vcf_file=       "{sample}.HaplotypeCaller.vcf"
    shell:
        '''
        cat <(echo 'out_dir = "{params.outdir}"') <(echo 'count_dir = "{params.count_dir}"') <(echo 'snv_dir = "{params.snv_dir}"')  > {params.config_file}
        cat <(echo '{params.sample} {params.htseq_file} {params.vcf_file}') > {params.meta_file}
        Rscript_GEP /home/zgu_labs/bin/R/RNAseq/RNAseqCNV.R {params.config_file} {params.meta_file} &> {log}
        '''

rule Done:
    input:
        HTSeq_count=    dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.DUX4patched.HTSeq",
        vcf=            dir_out + "/{sample}/TRANSCRIPTOME/Mutation/{sample}.HaplotypeCaller.vcf",
        cicero=         dir_out + "/{sample}/TRANSCRIPTOME/cicero/done.txt",
        fusioncatcher=  dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.FusionCatcher.bmk",
        RNAseqCNV=      dir_out + "/{sample}/TRANSCRIPTOME/log/{sample}.RNAseqCNV.bmk",
    output:
        dir_out + "/{sample}/TRANSCRIPTOME/log/Done.txt"
    params:
        dir_temp=dir_out + "/{sample}/TRANSCRIPTOME/temp/",
        mapping_file_prefix=dir_out + "/temp/TRANSCRIPTOME/log/*/{sample}",
        log_dir=dir_out + "/{sample}/TRANSCRIPTOME/log/"
    shell:
        '''
        touch {output}
        rm -rf {params.dir_temp}
        #mv {params.mapping_file_prefix}* {params.log_dir}
        '''



























