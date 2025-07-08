#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=3GB
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=10
#PBS -N STAR_alignment_unmapped
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

module load star/2.5.3a

STAR --genomeDir $INDEX --genomeLoad LoadAndExit

ls $INPUT | grep "[.]fq[.]gz$" | while read FILE; do
    SAMPLE=`basename "$FILE" _trimmed.fq.gz`
    STAR \
        --runThreadN 10 \
        --genomeDir $INDEX \
        --genomeLoad NoSharedMemory \
        --readFilesIn $INPUT/$FILE \
        --readFilesCommand zcat \
        --outFileNamePrefix $OUTPUT/$SAMPLE. \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outReadsUnmapped Fastx
done

STAR --genomeDir $INDEX --genomeLoad Remove

