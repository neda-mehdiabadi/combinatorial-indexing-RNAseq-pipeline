#-------------------------------------------------------------------------------
# Assign reads to genes, using the BED files as input
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts
GENE_MODEL_DIR=/group/card2/Neda/MCRI_LAB/must-do-projects/GENERAL_CODES

cd ${WORKING_DIR}
BATCH_SIZE=4

mkdir ${WORKING_DIR}/output/combo/file-lists-for-assign-reads-to-genes
mkdir ${WORKING_DIR}/output/combo/unique-read-to-gene-assignments

ls ${WORKING_DIR}/output/combo/aligned-reads-rmdup-split-bed/ | grep "[.]bed$" | split -l $BATCH_SIZE -d - ${WORKING_DIR}/output/combo/file-lists-for-assign-reads-to-genes/

ls ${WORKING_DIR}/output/combo/file-lists-for-assign-reads-to-genes | while read BATCH; do
    qsub $SCRIPTS_DIR/assign-reads-to-genes.sh       \
    -v INPUT_DIR=${WORKING_DIR}/output/combo/aligned-reads-rmdup-split-bed,FILE_LIST=${WORKING_DIR}/output/combo/file-lists-for-assign-reads-to-genes/$BATCH,EXON_BED=$GENE_MODEL_DIR/latest.exons.bed,GENE_BED=$GENE_MODEL_DIR/latest.genes.bed,SCRIPTS_DIR=$SCRIPTS_DIR/,OUTPUT_DIR=${WORKING_DIR}/output/combo/unique-read-to-gene-assignments
done