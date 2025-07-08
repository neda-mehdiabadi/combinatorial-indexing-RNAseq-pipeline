#-------------------------------------------------------------------------------
# Align reads using STAR
# Need to download STAR index
# Not included in github repository  HINT: 3 Gb ram for 20 minutes using 10 cores taken for 2 runs (now we have 192 runs!!!)
#-------------------------------------------------------------------------------
#####################################################################
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts
STAR_INDEX=/group/card2/Neda/MCRI_LAB/must-do-projects/GENERAL_CODES/refgenome/hg38_GRCh38.105.100_and_mm39_GRCh39.105.100.premRNA/star/

cd ${WORKING_DIR}

mkdir ${WORKING_DIR}/output/combo/aligned-reads

qsub $SCRIPTS_DIR/STAR-alignReads.sh \
     -v INPUT=$WORKING_DIR/output/combo/trimmed-fastq,INDEX=$STAR_INDEX,OUTPUT=$WORKING_DIR/output/combo/aligned-reads

