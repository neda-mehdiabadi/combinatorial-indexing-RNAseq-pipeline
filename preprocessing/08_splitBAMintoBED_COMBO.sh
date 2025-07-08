#-------------------------------------------------------------------------------
# Split reads in BAM files into BED intervals
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts

BATCH_SIZE=4
cd ${WORKING_DIR}
mkdir ${WORKING_DIR}/output/combo/file-lists-for-rmdup
mkdir ${WORKING_DIR}/output/combo/aligned-reads-rmdup-split-bed

ls ${WORKING_DIR}/output/combo/aligned-reads-filtered-sorted/ | grep "[.]bam$" | split -l $BATCH_SIZE -d - ${WORKING_DIR}/output/combo/file-lists-for-rmdup/

ls ${WORKING_DIR}/output/combo/file-lists-for-rmdup | while read BATCH; do
    qsub $SCRIPTS_DIR/rmdup-and-make-split-bed.sh    \
    -v INPUT_DIR=${WORKING_DIR}/output/combo/aligned-reads-filtered-sorted,FILE_LIST=${WORKING_DIR}/output/combo/file-lists-for-rmdup/$BATCH,SCRIPTS_DIR=$SCRIPTS_DIR/,OUTPUT_DIR=${WORKING_DIR}/output/combo/aligned-reads-rmdup-split-bed
done