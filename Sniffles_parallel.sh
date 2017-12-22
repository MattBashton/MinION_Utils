#!/bin/bash -eu

# Matthew Bashton 2017
# Takes a directory of input .fastq.gz and runs the Sniffles SV detector on them
# using GNU parallel.  Output is VCF.

hostname
date

tput bold
echo -ne "\nFinding all .bam and generating jobs list\n"
tput sgr0

# CPU cores for GNU parallel.
CPU=10

# Threads for Sniffles
THREADS=1

# Job_list_file
JL_FILE="Sniffles_jobs.txt"

# Remove existing jobs file if exists
if [ -e $JL_FILE ]
then
    rm $JL_FILE
fi

# Create VCF oputput dir
VCFOUT_DIR="../VCF_sniffles"
if [[ ! -e $VCFOUT_DIR ]]
then
    mkdir $VCFOUT_DIR
fi

# Make array of all bam
BAM_LIST=( $(ls -1 *.bam) )
NO_BAM=${#BAM_LIST[@]}
COUNT=1

# Generate jobs list
for BAM in ${BAM_LIST[@]}
do
    B_NAME=$(basename $BAM .bam)

    echo -ne "Found $COUNT of $NO_BAM, $BAM: creating Sniffles jobs"
    echo -ne "\n"
    # Note change parameters to suite your own experimental design / SV type
    echo "sniffles -t $THREADS -s 3 -q 5 -l 30 -m $B_NAME.bam -v $VCFOUT_DIR/$B_NAME.vcf" >> $JL_FILE
    ((COUNT++))
done

# Run jobs via GNU parallel
tput bold
echo -ne "\nRunning jobs via GNU parallel over $CPU cores\n"
tput sgr0
LOG_F=$(basename $JL_FILE .txt)
parallel --progress --jobs $CPU --joblog $LOG_F.log < $JL_FILE
echo -ne "\nDone! - Run logged to $LOG_F.log\n"
