#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )
ids=params.sraid


    
process sratool {
    input:

    output:

    script:

    """
    """

    }

process staralign {
    
    input:

    output:

    
    script:
    """
    STAR --genomeDir /fdb/STAR_indices/2.7.10b/UCSC/hg38/genes-100/ \
    --readFilesIn %s --readFilesCommand zcat --soloType CB_UMI_Simple \
    --soloCBwhitelist 3M-february-2018.txt --runThreadN 16 \
    --soloUMIlen 12 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
    --outSAMattributes NH HI nM AS CR CY CB UR UY UB RG sS sQ sM \
    --outFileNamePrefix %s --outSAMattrRGline \"ID:%s\" \
    --outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 210 \
    --limitBAMsortRAM 42186424882"

    sambamba index %s -@ 18

    """

}


process cellsnplite {
    input:

    output:


    script:

    """
    """
}

workflow {
    Channel
    .fromSRA(ids)
    .view()
    


    

}

/*

ids = ['ERR908507', 'ERR908506', 'ERR908505']
Channel
    .fromSRA(ids)
    .view()

module load sratoolkit
awk -F "," '{if (NR!=1) print "prefetch " $1'} SraRunTable.txt  >dl.sh
sh dl.sh

##error in SRR14037761
grep "SRR14037779\|SRR14037780\|SRR14037793\|SRR14037794" dl.sh >redl.sh

####Create loop to dump fastq.gz
#library(tidyverse)
#sra=read_csv("SraRunTable.txt")
outc=sprintf("/data/yesudhasd2/PRJNA716349-scRNA-bladder-cancer/darryl/parallel-fastq-dump/parallel-fastq-dump --sra-id %s/%s.sra --threads 8 --outdir %s/ --split-files --gzip",sra$Run,sra$Run,sra$Run)
write.table(outc,"outc.swarm",row.names=F,col.names=F,quote=F)
*/