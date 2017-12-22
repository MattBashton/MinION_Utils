#!/bin/bash -eu

# Matthew Bashton 2017
# Takes a directory of input .fastq.gz and runs the NanoSV on them
# using GNU parallel.  Also requires sambamba.  Output is VCF.

hostname
date

tput bold
echo -ne "\nFinding all .bam and generating jobs list\n"
tput sgr0

# CPU cores for GNU parallel.
CPU=10

# Job_list_file
JL_FILE="nanosv_jobs.txt"

# Remove existing jobs file if exists
if [ -e $JL_FILE ]
then
    rm $JL_FILE
fi

# Create VCF oputput dir
VCFOUT_DIR="../VCF_nanosv"
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
    
    echo -ne "Found $COUNT of $NO_BAM, $BAM: creating NanoSV jobs"
    echo -ne "\n"
    # Note change parameters to suite your own experimental design / SV type 
    echo "nanosv.pl -sambamba /usr/local/bin/sambamba -pid 0.8 -count 4 -distance 100 -pidf 0.8 -cluster 10 barcode01.bam -gap 100 > $VCFOUT_DIR/$B_NAME.vcf" >> $JL_FILE
    ((COUNT++))
done

# Run jobs via GNU parallel
tput bold
echo -ne "\nRunning jobs via GNU parallel over $CPU cores\n"
tput sgr0
LOG_F=$(basename $JL_FILE .txt)
parallel --progress --jobs $CPU --joblog $LOG_F.log < $JL_FILE
echo -ne "\nDone! - Run logged to $LOG_F.log\n"
