#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=1GB
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 10_duplicationRate
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

module load samtools/1.9
module load datamash/1.7

cat $BATCH | while read PCR_WELL; do
    awk '{
        split($4, arr, "|");
        if (!seen[arr[1]]) {
            seen[arr[1]] = 1;
            count[arr[2]]++;
        }
    } END {
        for (sample in count)
            print sample "\t" count[sample];
    }' $RMDUP_SPLIT_BED_DIR/$PCR_WELL.bed \
    | sort -k1,1 \
    >$OUTPUT_DIR/$PCR_WELL.UMI.count

    samtools view $FILTERED_READS_DIR/$PCR_WELL.bam \
    | cut -d '|' -f 2 \
    | datamash -g 1 count 1 \
    | sort -k1,1 -S 2G \
    | datamash -g 1 sum 2 \
    >$OUTPUT_DIR/$PCR_WELL.read.count

    echo "Processed $PCR_WELL"
done

