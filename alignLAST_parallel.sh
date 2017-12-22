#!/bin/bash -eu

# Matthew Bashton 2017
# Takes a directory of input .fastq.gz and runs the LAST aligner on them
# using GNU parallel.  Output is then converted to sroted index .bam

hostname
date

tput bold
echo -ne "\nFinding all fastq.gz and generating jobs list\n"
tput sgr0

# CPU cores for GNU parallel.
CPU=10

# Threads for last-train
LTHREADS=10

# Threads for LAST and sambamba
THREADS=1

# Index file for LAST
IDX="/home/bashton/Ens91/GRCh38.PA.Ens91.LAST"

# Job_list_file
JL_FILE="LAST_jobs.txt"

# Remove existing jobs file if exists
if [ -e $JL_FILE ]
then
    rm $JL_FILE
fi

# Create BAM oputput dir
BAMOUT_DIR="../BAM_LAST"
if [[ ! -e $BAMOUT_DIR ]]
then
    mkdir $BAMOUT_DIR
fi

# Make array of all fastq.gz
FASTQ_LIST=( $(ls -1 *.fastq.gz) )
NO_FASTQ=${#FASTQ_LIST[@]}
COUNT=1

# If no LAST_PARAMS file then generate one using the first file in the list
if [[ ! -e LAST_params ]]
then
    echo -ne "\n - Running last-train on $LTHREADS threads, with index $IDX and fastq file: ${FASTQ_LIST[0]}\n"
    last-train -Q1 -P $LTHREADS $IDX ${FASTQ_LIST[0]} > LAST_params
    echo -ne " - Output saved to LAST_params\n\n"
fi


# Generate jobs list
for FASTQ in ${FASTQ_LIST[@]}
do
    B_NAME=$(basename $FASTQ .fastq.gz)
    echo -ne "Found $COUNT of $NO_FASTQ, $FASTQ: creating LAST and sambamba jobs"
    echo -ne "\n"
    echo "lastal -Q1 -P $THREADS -p LAST_params $IDX $FASTQ | last-split > $B_NAME.maf; maf-convert -d sam $B_NAME.maf | sambamba view -S -f=bam /dev/stdin | sambamba sort /dev/stdin -o $BAMOUT_DIR/$B_NAME.bam; rm $B_NAME.maf; sambamba index $BAMOUT_DIR/$B_NAME.bam; sambamba index -t $THREADS $BAMOUT_DIR/$B_NAME.bam " >> $JL_FILE
    ((COUNT++))
done

# Run jobs via GNU parallel
tput bold
echo -ne "\nRunning jobs via GNU parallel over $CPU cores\n"
tput sgr0
LOG_F=$(basename $JL_FILE .txt)
parallel --progress --jobs $CPU --joblog $LOG_F.log < $JL_FILE
echo -ne "\nDone! - Run logged to $LOG_F.log\n"
