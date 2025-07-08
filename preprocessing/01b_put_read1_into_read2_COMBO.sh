#!/bin/bash

##########################
#                        #
#   The PBS directives   #
#                        #
##########################

# Define the shell in which your jobs should run. Shouldn't really be changed
# unless you have a very specific reason for doing so
#PBS -S /bin/bash


# Define the name for the job
#PBS -N 01a_put_read1_into_read2_RENAME_COMBO fastq



# Defining the wall time for the job
#PBS -l walltime=2:00:00




# Selecting which queue to send the job to
#PBS -q batch





# Defining the amount of memory you require
#PBS -l mem=1GB



# Defining email notifications
## a = notify when job aborts (default)
## b = notify when job begins
## e = notify when job terminates
## n = no mail at all (If specified takes precedence)
#PBS -m aeb

# Define the email address to be used in correspondence
#PBS -M neda90rahmani@gmail.com


# Define the number of nodes and cores you require
#PBS -l nodes=1:ppn=1


# Used to define which project job is associated with
#PBS -A card2

### END OF PBS OPTIONS

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#echo ------------------------------------------------------
#echo -n 'Job is running on node '; cat $PBS_NODEFILE
#echo ------------------------------------------------------
#echo PBS: qsub was run on $PBS_O_HOST
#echo PBS: originating queue is $PBS_O_QUEUE
#echo PBS: executing queue is $PBS_QUEUE
#echo PBS: working directory is $PBS_O_WORKDIR
#echo PBS: execution mode is $PBS_ENVIRONMENT
#echo PBS: job identifier is $PBS_JOBID
#echo PBS: job name is $PBS_JOBNAME
#echo PBS: node file is $PBS_NODEFILE
#echo PBS: current home directory is $PBS_O_HOME
#echo PBS: temporary directory on node is $TMPDIR
#echo PBS: PATH = $PBS_O_PATH
#echo ------------------------------------------------------

#-------------------------------------------------------------------------------
# Put read 1 info (RT well, UMI) into read 2 read name
# 3 Level sci-RNA Seq 
#-------------------------------------------------------------------------------
#####################################################################
WORKING_DIR=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv
PATH_TO_GIT_REPO=/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master
SCRIPTS_DIR=$PATH_TO_GIT_REPO/process_from_raw/sci-RNA-seq-pipeline-scripts
RT_BARCODES_FILE=$SCRIPTS_DIR/RT_384_3lvl
LIG_BARCODES_FILE=$SCRIPTS_DIR/ligation_384_3lvl

SAMPLE_NAME="sciPlex"
echo "$SAMPLE_NAME" > ${WORKING_DIR}/output/combinatorial.indexing.key
BATCH_SIZE=4
cd ${WORKING_DIR}

mkdir ${WORKING_DIR}/output/combo/combined-fastq
mkdir ${WORKING_DIR}/output/combo/file-lists-for-r1-info-munging
mkdir ${WORKING_DIR}/output/combo/put-r1-info-in-r2-logs

ls ${WORKING_DIR}/data/combo/ | grep _R1 | grep -v Undetermined | split -l $BATCH_SIZE -d - ${WORKING_DIR}/output/combo/file-lists-for-r1-info-munging/
   
ls ${WORKING_DIR}/output/combo/file-lists-for-r1-info-munging | while read BATCH; do
   qsub $SCRIPTS_DIR/put-read1-info-in-read2_3Level.sh \
   -v INPUT_DIR=$WORKING_DIR/data/combo,R1_FILE_LIST=$WORKING_DIR/output/combo/file-lists-for-r1-info-munging/$BATCH,SCRIPTS_DIR=$SCRIPTS_DIR/,RT_OLIGO_LIST=$RT_BARCODES_FILE,LIG_OLIGO_LIST=$LIG_BARCODES_FILE,INDEXING_KEY=$WORKING_DIR/output/combinatorial.indexing.key,OUTPUT_DIR=$WORKING_DIR/output/combo/combined-fastq,LOGS_DIR=$WORKING_DIR/output/combo/put-r1-info-in-r2-logs
   done





