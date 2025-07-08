#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=1GB
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 09_assign_readsTo_genes_usingBED
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com


module load bedtools/2.29.0
module load datamash/1.7
module load python/3.7.3

cat $FILE_LIST | while read FILE; do
    PCR_WELL=`basename "$FILE" .bed`

    bedtools map \
        -a $INPUT_DIR/$FILE \
        -b $EXON_BED \
        -nonamecheck -s -f 0.95 -c 5 -o distinct -delim '|' \
    | bedtools map \
        -a - -b $GENE_BED \
        -nonamecheck -s -f 0.95 -c 5 -o distinct -delim '|' \
    | sort -k4,4 -k2,2n -k3,3n -S 3G \
    | datamash -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8 \
    | python $SCRIPTS_DIR/assign-reads-to-genes.py $GENE_BED \
    >$OUTPUT_DIR/$PCR_WELL

done

