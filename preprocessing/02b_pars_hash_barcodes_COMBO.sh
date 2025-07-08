#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=1GB
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 02b_parse_hash_barcodes
#PBS -q batch
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

#-------------------------------------------------------------------------------
# Parse Hash Barcodes
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts
HASH_BARCODES_FILE=$SCRIPTS_DIR/hashSampleSheet.txt

SAMPLE_NAME="sciPlex"

cd ${WORKING_DIR}

mkdir ${WORKING_DIR}/output/combo/hashRDS

module load datamash/1.7
module load R/3.6.1

zcat ${WORKING_DIR}/output/combo/hash-fastq/*.gz \
    | datamash -g 2,3,4 count 3 \
    | datamash -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > ${WORKING_DIR}/output/combo/hashRDS/hashReads.per.cell


zcat ${WORKING_DIR}/output/combo/hash-fastq/*.gz \
    | uniq \
    | datamash -g 2,3,4 count 3 \
    | datamash -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > ${WORKING_DIR}/output/combo/hashRDS/hashUMIs.per.cell


Rscript $SCRIPTS_DIR/knee-plot.R            \
    ${WORKING_DIR}/output/combo/hashRDS/hashUMIs.per.cell          \
    ${WORKING_DIR}/output/combo/hashRDS


zcat ${WORKING_DIR}/output/combo/hash-fastq/*.gz \
    | uniq \
    | datamash -g 1,2,4,5 count 3  \
    > ${WORKING_DIR}/output/combo/hashRDS/hashTable.out 


paste ${WORKING_DIR}/output/combo/hashRDS/hashUMIs.per.cell  ${WORKING_DIR}/output/combo/hashRDS/hashReads.per.cell \
     | cut -f 1,2,6,3 \
     | awk 'BEGIN {OFS="\t";} {dup = 1-($3/$4); print $1,$2,$3,$4,dup;}' \
     > ${WORKING_DIR}/output/combo/hashRDS/hashDupRate.txt