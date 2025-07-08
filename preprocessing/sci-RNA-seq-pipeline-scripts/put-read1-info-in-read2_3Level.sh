#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=1GB
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 01b_read1_info_in_read2_3Level
#PBS -q batch
#PBS -A card2

# notifications for events: b for job begin, e for end, or a when job is aborted
#PBS -m abe
#PBS -M neda90rahmani@gmail.com


BATCH_ID=`basename "$R1_FILE_LIST"`

cat $R1_FILE_LIST | while read R1_FILE; do
    R2_FILE=`echo "$R1_FILE" | sed 's/_R1_/_R2_/'`
    PCR_COMBO=`echo "$R1_FILE" | cut -d '_' -f 1`
    paste \
        <(gunzip <$INPUT_DIR/$R1_FILE) \
        <(gunzip <$INPUT_DIR/$R2_FILE) \
    | awk -f $SCRIPTS_DIR/put-read1-info-in-read2_3Level.awk -v PCR_COMBO="$PCR_COMBO" \
        $RT_OLIGO_LIST $LIG_OLIGO_LIST $INDEXING_KEY - \
    | gzip >$OUTPUT_DIR/$PCR_COMBO.fastq.gz

    echo "Processed $PCR_COMBO">$OUTPUT_DIR/processedPCRcombo.txt
done 2>$LOGS_DIR/$BATCH_ID