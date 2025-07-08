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
#PBS -N 00_fastqc



# Defining the wall time for the job
#PBS -l walltime=24:00:00




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


#>>>>>>>>>> Steps to excecute your job goes here <<<<<<<<<<<<<<<<<<<<<

#-------------------------------------------------------------------------------
# fastqc *.fastq.gz for both short(hashing molecules) and large (RNA libraries) fragments
#-------------------------------------------------------------------------------
#####################################################################
WORKING_DIR=/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-5th-DSP-TPM

cd ${WORKING_DIR}
module load fastqc/0.11.8

mkdir ${WORKING_DIR}/output/large_fragment_fastqc_results-firstBatch
#mkdir ${WORKING_DIR}/output/large_fragment_fastqc_R1_fastq_logs-firstBatch

ls ${WORKING_DIR}/data/230329_A01221_0173_AHKG2GDMXY-large-fragment-firstBatch/ | grep '\_R1.fastq.gz$' > ${WORKING_DIR}/output/R1_files.txt

cat $WORKING_DIR/output/R1_files.txt | while read R1_FILE; do
    R2_FILE=`echo "$R1_FILE" | sed 's/_R1/_R2/'`
    fastqc $WORKING_DIR/data/230329_A01221_0173_AHKG2GDMXY-large-fragment-firstBatch/$R1_FILE -o $WORKING_DIR/output/large_fragment_fastqc_results-firstBatch/
    fastqc $WORKING_DIR/data/230329_A01221_0173_AHKG2GDMXY-large-fragment-firstBatch/$R2_FILE -o $WORKING_DIR/output/large_fragment_fastqc_results-firstBatch/
done 





