
#-------------------------------------------------------------------------------
# Parse Hash Barcodes
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts
HASH_BARCODES_FILE=$SCRIPTS_DIR/hashSampleSheet.txt

cd ${WORKING_DIR}
mkdir ${WORKING_DIR}/output/combo/hash-file-lists-for-trimming
mkdir ${WORKING_DIR}/output/combo/hash-fastq
mkdir ${WORKING_DIR}/output/combo/hash-logs
BATCH_SIZE=1
ls $WORKING_DIR/output/combo/combined-fastq/ | split -l $BATCH_SIZE -d - ${WORKING_DIR}/output/combo/hash-file-lists-for-trimming/
ls $WORKING_DIR/output/combo/hash-file-lists-for-trimming | while read BATCH; do
    qsub $SCRIPTS_DIR/parse_hash.sh \
    -v INPUT_DIR=$WORKING_DIR/output/combo/combined-fastq,R1_FILE_LIST=$WORKING_DIR/output/combo/hash-file-lists-for-trimming/$BATCH,SCRIPTS_DIR=$SCRIPTS_DIR/,RT_OLIGO_LIST=$HASH_BARCODES_FILE,INDEXING_KEY=$WORKING_DIR/output/combinatorial.indexing.key,OUTPUT_DIR=$WORKING_DIR/output/combo/hash-fastq,LOGS_DIR=$WORKING_DIR/output/combo/hash-logs
    done


