#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=3GB
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 13c_Make_final_UMI_count_matrix
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

#-------------------------------------------------------------------------------
# MAKE UMI COUNT MATRIX
#-------------------------------------------------------------------------------

WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts
GENE_MODEL_DIR=/group/card2/Neda/MCRI_LAB/must-do-projects/GENERAL_CODES

cd ${WORKING_DIR}

module load R/3.6.1
source $SCRIPTS_DIR/bin/activate

Rscript $SCRIPTS_DIR/makeCDS.R \
        $WORKING_DIR/output/combo/UMI.count.matrix \
        $GENE_MODEL_DIR/gene.annotations \
        $WORKING_DIR/output/combo/cell.annotations \
        $WORKING_DIR/output/combo
