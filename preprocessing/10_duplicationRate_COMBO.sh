#-------------------------------------------------------------------------------
# Compute the duplication rate and proportion of reads that are from rRNA
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts

cd ${WORKING_DIR}
BATCH_SIZE=4

module load datamash/1.7

mkdir ${WORKING_DIR}/output/combo/UMI-counts-by-sample
mkdir ${WORKING_DIR}/output/combo/file-lists-for-UMI-counting

ls ${WORKING_DIR}/output/combo/aligned-reads-rmdup-split-bed/ | while read FILE; do
    PCR_WELL=`basename $FILE .bed`
    echo "$PCR_WELL"
done \
| split -l $BATCH_SIZE -d - ${WORKING_DIR}/output/combo/file-lists-for-UMI-counting/

ls ${WORKING_DIR}/output/combo/file-lists-for-UMI-counting | while read BATCH; do
    qsub $SCRIPTS_DIR/count-UMI-per-sample.sh    \
    -v FILTERED_READS_DIR=${WORKING_DIR}/output/combo/aligned-reads-filtered-sorted,RMDUP_SPLIT_BED_DIR=${WORKING_DIR}/output/combo/aligned-reads-rmdup-split-bed,BATCH=${WORKING_DIR}/output/combo/file-lists-for-UMI-counting/$BATCH,OUTPUT_DIR=${WORKING_DIR}/output/combo/UMI-counts-by-sample
done

