#-------------------------------------------------------------------------------
# Make knee plots -- Sanjay
#-------------------------------------------------------------------------------
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts

cd ${WORKING_DIR}


module load python/3.7.3
module load R/3.3.2
source $SCRIPTS_DIR/bin/activate

START_TIME=$SECONDS
cat ${WORKING_DIR}/output/combo/unique-read-to-gene-assignments/* \
    | python  $SCRIPTS_DIR/count_UMIs.py --cell - \
    > ${WORKING_DIR}/output/combo/UMIs.per.cell.barcode
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Processed $FILE in $ELAPSED_TIME seconds"


mkdir ${WORKING_DIR}/output/combo/knee-plots

Rscript $SCRIPTS_DIR/knee-plot.R            \
    ${WORKING_DIR}/output/combo/UMIs.per.cell.barcode                   \
    ${WORKING_DIR}/output/combo/knee-plots
