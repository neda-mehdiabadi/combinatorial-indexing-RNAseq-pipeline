#!/bin/bash

# first line tells the shell program what prgram to use, in this case bash
# Lines starting with # are usually comments, excepting the above and #PBS lines
# lines with #PBS are metadata for queue regading resources etc. required

#PBS -l mem=1GB
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=1
#PBS -N 03_Trim_poly_A_tails_COMBO
#PBS -q batch
#PBS -A card2
# notifications for events: b for job begin, e for end, or a when job is aborted

#PBS -m abe
#PBS -M neda90rahmani@gmail.com

module load perl/5.30.0

module load python/3.7.3
module load cutadapt

module load trimgalore/0.6.5

cat $FILE_LIST | while read FILE; do
    trim_galore $INPUT_DIR/$FILE -a AAAAAAAA --three_prime_clip_R1 1 --gzip -o $OUTPUT_DIR
done

