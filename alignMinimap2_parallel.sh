#!/bin/bash -eu

# Matthew Bashton 2017
# Takes a directory of input .fastq.gz and runs the Minimap2 aligner on them
# using GNU parallel.  Output is then converted to sroted index .bam

hostname
date

tput bold
echo -ne "\nFinding all fastq.gz and generating jobs list\n"
tput sgr0

# CPU cores for GNU parallel.
CPU=32

# Threads for Minmap2
THREADS=1

# Index file for Minmap2
IDX="~/Ens91/GRCh38.PA.Ens91.mmi"

# Job_list_file
JL_FILE="Minimap2_jobs.txt"

# Remove existing jobs file if exists
if [ -e $JL_FILE ]
then
    rm $JL_FILE
fi

# Create BAM oputput dir
BAMOUT_DIR="../BAM_Minmap2_B"
if [[ ! -e $BAMOUT_DIR ]]
then
    mkdir $BAMOUT_DIR
fi

# Make array of all fastq.gz
FASTQ_LIST=( $(ls -1 *.fastq.gz) )
NO_FASTQ=${#FASTQ_LIST[@]}
COUNT=1

# Generate jobs list
for FASTQ in ${FASTQ_LIST[@]}
do
    B_NAME=$(basename $FASTQ .fastq.gz)
    echo -ne "Found $COUNT of $NO_FASTQ, $FASTQ: creating Minimap2 and samtools jobs"
    echo -ne "\n"
    echo "minimap2 -ax map-ont -L $IDX $FASTQ > $B_NAME.sam; samtools sort $B_NAME.sam > $BAMOUT_DIR/$B_NAME.bam; samtools index $BAMOUT_DIR/$B_NAME.bam; rm $B_NAME.sam" >> $JL_FILE
    ((COUNT++))
done

# Run jobs via GNU parallel
tput bold
echo -ne "\nRunning jobs via GNU parallel over $CPU cores\n"
tput sgr0
LOG_F=$(basename $JL_FILE .txt)
parallel --progress --jobs $CPU --joblog $LOG_F.log < $JL_FILE
echo -ne "\nDone! - Run logged to $LOG_F.log\n"
