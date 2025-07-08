#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=1GB
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 02a_parse_hash_barcodes
#PBS -q batch
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

BATCH_ID=`basename "$R1_FILE_LIST"`
cat $R1_FILE_LIST | while read R1_FILE; do
    PCR_COMBO=`echo "$R1_FILE" | cut -d '.' -f 1`

    zcat $INPUT_DIR/$R1_FILE \
    | awk -f $SCRIPTS_DIR/parseHash.awk -v PCR_COMBO="$PCR_COMBO" \
        $RT_OLIGO_LIST $INDEXING_KEY - \
    | sed -e 's/|/,/g' \
    | awk 'BEGIN {FS=","; OFS="\t";} {print $2,$3"_"$4"_"$5,$6,$7,$8}' \
    | sort -S 16G -k1,1 -k2,2 -k4,4 -k3,3 \
    | gzip > $OUTPUT_DIR/$PCR_COMBO.hash.gz

    echo "Processed $PCR_COMBO">${INPUT_DIR}/processed.txt
done 2> $LOGS_DIR/$BATCH_ID
#To redirect stderr, we do: 2>
