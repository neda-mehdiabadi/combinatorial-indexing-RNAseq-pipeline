#-------------------------------------------------------------------------------
# Align reads using STAR
# Need to download STAR index
# Not included in github repository  HINT: 3 Gb ram for 20 minutes using 10 cores taken for 2 runs (now we have 192 runs!!!)
#-------------------------------------------------------------------------------
#####################################################################
WORKING_DIR=/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-1st-30.04.21/output
PATH_TO_GIT_REPO=/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-1st-30.04.21/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts
STAR_INDEX=/group/card2/Neda/MCRI_LAB/must-do-projects/GENERAL_CODES/star/hg38_GRCh38.100.premRNA/star/

cd ${WORKING_DIR}

mkdir unmapped-aligned-reads

qsub $SCRIPTS_DIR/STAR-alignReads-unmapped.sh \
     -v INPUT=$WORKING_DIR/trimmed-fastq,INDEX=$STAR_INDEX,OUTPUT=$WORKING_DIR/unmapped-aligned-reads

