#-------------------------------------------------------------------------------
# Filter ambiguously-mapped reads and sort BAM files
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts

cd ${WORKING_DIR}
BATCH_SIZE=4

mkdir ${WORKING_DIR}/output/combo/file-lists-for-samtools-sort
mkdir ${WORKING_DIR}/output/combo/aligned-reads-filtered-sorted

ls ${WORKING_DIR}/output/combo/aligned-reads/ | grep "[.]Aligned[.]out[.]bam$" | split -l $BATCH_SIZE -d - ${WORKING_DIR}/output/combo/file-lists-for-samtools-sort/

ls ${WORKING_DIR}/output/combo/file-lists-for-samtools-sort | while read BATCH; do
    qsub $SCRIPTS_DIR/samtools-filter-sort.sh    \
    -v INPUT_DIR=${WORKING_DIR}/output/combo/aligned-reads,FILE_LIST=${WORKING_DIR}/output/combo/file-lists-for-samtools-sort/$BATCH,OUTPUT_DIR=${WORKING_DIR}/output/combo/aligned-reads-filtered-sorted
done


