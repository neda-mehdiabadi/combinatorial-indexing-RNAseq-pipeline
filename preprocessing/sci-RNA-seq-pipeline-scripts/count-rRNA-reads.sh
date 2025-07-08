#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=15GB
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=10
#PBS -N 07_countrRNA
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

module load bedtools/2.29.0

cat $FILE_LIST | while read FILE; do
    PCR_WELL=`basename "$FILE" .Aligned.out.bam`

    bedtools intersect -a $INPUT_DIR/$FILE -b $RRNA_BED -c -nonamecheck -bed \
    | awk '{
        split($4, arr, "|");
        if (!seen[arr[1]]) {
            if ($NF > 0)
                rrna_count[arr[2]]++;
            total_count[arr[2]]++;
            seen[arr[1]] = 1;
        }
    } END {
        for (sample in total_count)
            printf "%s\t%d\t%d\n",
                sample, rrna_count[sample], total_count[sample];
    }' \
    >$OUTPUT_DIR/$PCR_WELL

    echo "Processed $FILE">${INPUT_DIR}/processed.txt
done

