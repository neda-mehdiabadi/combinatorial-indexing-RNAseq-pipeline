#-------------------------------------------------------------------------------
# Trim poly-A tails for large (RNA libraries) fragment
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts

cd ${WORKING_DIR}
mkdir ${WORKING_DIR}/output/combo/trimmed-fastq
mkdir ${WORKING_DIR}/output/combo/file-lists-for-trimming
BATCH_SIZE=4

ls ${WORKING_DIR}/output/combo/combined-fastq/ | split -l $BATCH_SIZE -d - ${WORKING_DIR}/output/combo/file-lists-for-trimming/

ls ${WORKING_DIR}/output/combo/file-lists-for-trimming | while read BATCH; do
    qsub $SCRIPTS_DIR/run-trim-galore.sh \
    -v INPUT_DIR=$WORKING_DIR/output/combo/combined-fastq,FILE_LIST=$WORKING_DIR/output/combo/file-lists-for-trimming/$BATCH,OUTPUT_DIR=$WORKING_DIR/output/combo/trimmed-fastq
    done
