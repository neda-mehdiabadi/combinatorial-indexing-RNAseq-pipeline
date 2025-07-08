#-------------------------------------------------------------------------------
# MAKE UMI COUNT MATRIX
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts

cd ${WORKING_DIR}
BATCH_SIZE=4

mkdir ${WORKING_DIR}/output/combo/file-lists-for-UMI-count-rollup
mkdir ${WORKING_DIR}/output/combo/UMI-count-rollup

ls ${WORKING_DIR}/output/combo/unique-read-to-gene-assignments/ | split -l $BATCH_SIZE -d - ${WORKING_DIR}/output/combo/file-lists-for-UMI-count-rollup/

ls ${WORKING_DIR}/output/combo/file-lists-for-UMI-count-rollup | while read BATCH; do
   qsub $SCRIPTS_DIR/UMI-count-rollup.sh        \
   -v INPUT_DIR=${WORKING_DIR}/output/combo/unique-read-to-gene-assignments,FILE_LIST=${WORKING_DIR}/output/combo/file-lists-for-UMI-count-rollup/$BATCH,OUTPUT_DIR=${WORKING_DIR}/output/combo/UMI-count-rollup
done
