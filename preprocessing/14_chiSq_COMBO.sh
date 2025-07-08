#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=15GB
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 14_chiSq
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

#-------------------------------------------------------------------------------
# MAKE UMI COUNT MATRIX
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Using the knee plots, define cells that are background cells and 
# those cells that are intermediate cells
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts
HASH_BARCODES_FILE=$SCRIPTS_DIR/hashSampleSheet.txt

cd ${WORKING_DIR}
BATCH_SIZE=4

module load R/3.6.1
 
UMI_PER_CELL_CUTOFF=750
UMI_PER_CELL_LOWER=300


 awk '$3 >= CUTOFF { print $2 }' CUTOFF=$UMI_PER_CELL_CUTOFF $WORKING_DIR/output/combo/UMIs.per.cell.barcode \
    > $WORKING_DIR/output/combo/hashRDS/real.cells.csv

awk '$3 <= CUTOFF_LOW {print $2}' CUTOFF_LOW=$UMI_PER_CELL_LOWER $WORKING_DIR/output/combo/UMIs.per.cell.barcode \
    > $WORKING_DIR/output/combo/hashRDS/background.cells.csv


awk '($3 > CUTOFF_LOW) && ($3 < CUTOFF_HIGH) { print $2}
    ' CUTOFF_LOW=$UMI_PER_CELL_LOWER CUTOFF_HIGH=$UMI_PER_CELL_CUTOFF $WORKING_DIR/output/combo/UMIs.per.cell.barcode \
    > $WORKING_DIR/output/combo/hashRDS/intermediate.cells.csv


Rscript $SCRIPTS_DIR/chiSq_2lvl.R                           \
        $WORKING_DIR/output/combo/hashRDS/background.cells.csv           \
        $WORKING_DIR/output/combo/hashRDS/real.cells.csv                 \
        $WORKING_DIR/output/combo/hash-fastq/                          \
        $HASH_BARCODES_FILE                                 \
        $WORKING_DIR/output/combo/cds.RDS                                \
        $WORKING_DIR/output/combo/

