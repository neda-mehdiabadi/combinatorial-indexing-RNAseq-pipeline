#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=3GB
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 13b_Make_final_UMI_count_matrix
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
BATCH_SIZE=4

cat ${WORKING_DIR}/output/combo/UMI-count-rollup/* | gzip > ${WORKING_DIR}/output/combo/prelim.UMI.count.rollup.gz

#
# Make samples.to.exclude file.
# Each line lists a sample name to exclude.
# You can leave it empty.
#

touch ${WORKING_DIR}/output/combo/samples.to.exclude

#
# This uses the UMI_PER_CELL_CUTOFF variable
# defined in the knee plot section.
#

UMI_PER_CELL_CUTOFF=0
echo "UMI_PER_CELL_CUTOFF = $UMI_PER_CELL_CUTOFF"
module load datamash/1.7

gunzip < ${WORKING_DIR}/output/combo/prelim.UMI.count.rollup.gz \
| datamash -g 1 sum 3 \
| tr '|' '\t' \
| awk -v CUTOFF=$UMI_PER_CELL_CUTOFF '
    ARGIND == 1 {
        exclude[$1] = 1;
    } $3 >= CUTOFF && !($1 in exclude) {
        print $2 "\t" $1; 
    }' ${WORKING_DIR}/output/combo/samples.to.exclude - \
| sort -k1,1 -S 4G \
> ${WORKING_DIR}/output/combo/cell.annotations

gunzip < ${WORKING_DIR}/output/combo/prelim.UMI.count.rollup.gz \
| tr '|' '\t' \
| awk '{
    if (ARGIND == 1) {
        gene_idx[$1] = FNR;
    } else if (ARGIND == 2) {
        cell_idx[$1] = FNR;
    } else if ($2 in cell_idx) {
        printf "%d\t%d\t%d\n",
            gene_idx[$3], cell_idx[$2], $4;
    }
}'  $GENE_MODEL_DIR/gene.annotations  ${WORKING_DIR}/output/combo/cell.annotations - \
> ${WORKING_DIR}/output/combo/UMI.count.matrix