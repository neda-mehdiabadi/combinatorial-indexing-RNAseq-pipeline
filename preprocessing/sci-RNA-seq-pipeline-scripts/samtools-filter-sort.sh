#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=15GB
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=10
#PBS -N 06_filterreads_sortBAM
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

#The following 2 lines are needed for temp memories to be stored in a space with enough space, otherwise they will be stored in your personal drive with limited space and will send an error.
WORKING_DIR=/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-6th-TPM1-DSPKI
cd ${WORKING_DIR}

module load samtools/1.9

cat $FILE_LIST | while read FILE; do
    SAMPLE=`basename "$FILE" .Aligned.out.bam`

    samtools view -bh -q 30 -F 4 $INPUT_DIR/$FILE \
    | samtools sort -@ 10 - \
    >$OUTPUT_DIR/$SAMPLE.bam

done

