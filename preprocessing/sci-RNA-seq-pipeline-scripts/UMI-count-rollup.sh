#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=1GB
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 13a_Make_final_UMI_count_matrix
#PBS -q batch
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

module load datamash/1.7

cat $FILE_LIST | while read FILE; do
    awk '$3 == "exonic" || $3 == "intronic" {
        split($1, arr, "|");
        printf "%s|%s_%s_%s\t%s\n",
            arr[2], arr[3], arr[4], arr[5], $2;
    }' $INPUT_DIR/$FILE \
    | sort -k1,1 -k2,2 -S 2G \
    | datamash -g 1,2 count 2 \
    >$OUTPUT_DIR/$FILE
done

