#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=3GB
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=16
#PBS -N 04_align_with_STAR
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

WORKING_DIR=/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-6th-TPM1-DSPKI
cd ${WORKING_DIR}

module load star/2.5.3a

STAR --genomeDir $INDEX --genomeLoad LoadAndExit

ls $INPUT | grep "[.]fq[.]gz$" | while read FILE; do
    SAMPLE=`basename "$FILE" _trimmed.fq.gz`
    STAR \
        --runThreadN 16 \
        --genomeDir $INDEX \
        --genomeLoad LoadAndKeep \
        --readFilesIn $INPUT/$FILE \
        --readFilesCommand zcat \
        --outFileNamePrefix $OUTPUT/$SAMPLE. \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif 
done

STAR --genomeDir $INDEX --genomeLoad Remove

