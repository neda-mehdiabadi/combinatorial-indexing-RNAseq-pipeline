#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=1GB
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 08_splitBAMintoBED
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

module load samtools/1.9
module load bedtools/2.29.0
#module purge
#module load modules modules-init modules-gs
#module load pypy/3.5.6.0
#source /net/trapnell/vol1/home/sanjays/bin/sciRNAseq/pipeline_scripts/pypy_env/bin/activate

# this makes "sort" case sensitive
export LC_ALL=C
START_TIME=$SECONDS
cat $FILE_LIST | while read FILE; do
    PCR_WELL=`basename "$FILE" .bam`

    samtools view -h $INPUT_DIR/$FILE \
    | awk -f $SCRIPTS_DIR/rmdup.awk \
    | samtools view -bh \
    | bedtools bamtobed -i - -split \
    | sort -k1,1 -k2,2n -k3,3n -S 3G \
    >$OUTPUT_DIR/$PCR_WELL.bed

    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    echo "Processed $FILE in $ELAPSED_TIME seconds"
done

