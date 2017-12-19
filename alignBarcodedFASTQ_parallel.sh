#!/bin/bash -eu

# Matthew Bashton 2017
# Takes barcoded Albacore output and give you unified FASTQ.gz and
# sorted index BAM files per barcode using minimap2 and samtools.

# Needs to be run in the pass/ directory of your workspace/ dir of Albacore
# output.  First stage creates a jobs list which is supplied to GNU parallel in
# the second stage.

# Tested with minimap2.6 and samtools 1.6 also requires pigz parallel Gzip
# gzip can be used in place if need be.

hostname
date

tput bold
echo -ne "\nFinding all FASTQ and generating jobs list\n"
tput sgr0

# CPU cores
CPU=10

# gzip
GZIP="pigz"

# Index file for minimap2
IDX="~/Ens91/GRCh38.PA.Ens91.mmi"

# Save PWD
ORIG_DIR=$PWD

# Job_list_file
JL_FILE="minimap2_jobs.txt"

# Remove existing jobs file if exists
if [ -e $JL_FILE ]
then
    rm $JL_FILE
fi

# Create BAM oputput dir
BAMOUT_DIR="../BAM"
if [[ ! -e $BAMOUT_DIR ]]
then
    mkdir $BAMOUT_DIR
fi

# Create FASTQ output dir
FQOUT_DIR="../FASTQ"
if [[ ! -e $FQOUT_DIR ]]
then
    mkdir $FQOUT_DIR
fi

# Make array of all FASTQ locations
DIR_LIST=( $(ls -d */ | sed 's/\///') )
NO_DIR=${#DIR_LIST[@]}
COUNT=1

# Generate jobs list
for DIR in ${DIR_LIST[@]}
do
    echo -ne "Working on dir $COUNT of $NO_DIR, $DIR: "
    cd $DIR
    FQ_LIST=( $(ls -1 *.fastq 2>/dev/null) )
    NO_FQ=${#FQ_LIST[@]}
    echo -ne "found $NO_FQ FASTQ file(s)"
    if [ $NO_FQ -eq 0 ]
    then
	echo -ne " - no FASTQ for this barcode, skiping"
    elif
	[ $NO_FQ -eq 1 ]
    then
	echo -ne " - only one FASTQ file: "
	B_NAME=$(basename $FQ_LIST .fastq)
	echo -ne "$FQ_LIST"
	echo "cd $PWD; cat $B_NAME.fastq | $GZIP -fc > $DIR.fastq.gz; minimap2 -ax map-ont -L $IDX $DIR.fastq.gz > $DIR.sam; samtools sort $DIR.sam > ../$BAMOUT_DIR/$DIR.bam; samtools index ../$BAMOUT_DIR/$DIR.bam; mv $DIR.fastq.gz ../$FQOUT_DIR/; rm *.sam; rm *.fastq; cd $ORIG_DIR" >> $ORIG_DIR/$JL_FILE

    elif
	[ $NO_FQ -gt 1 ]
    then
	echo -ne " - more than one FASTQ file:"
	for FQ in ${FQ_LIST[@]}
	do
	    echo -ne " $FQ"
	done
	echo "cd $PWD; cat *.fastq | $GZIP -fc > $DIR.fastq.gz; minimap2 -ax map-ont -L $IDX $DIR.fastq.gz > $DIR.sam; samtools sort $DIR.sam > ../$BAMOUT_DIR/$DIR.bam; samtools index ../$BAMOUT_DIR/$DIR.bam; mv $DIR.fastq.gz ../$FQOUT_DIR/; rm *.sam; rm *.fastq; cd $ORIG_DIR" >> $ORIG_DIR/$JL_FILE
    fi
    cd $ORIG_DIR
    echo -ne "\n"
    ((COUNT++))
done

# Run jobs via GNU parallel
tput bold
echo -ne "\nRunning jobs via GNU parallel over $CPU cores\n"
tput sgr0
LOG_F=$(basename $JL_FILE .txt)
parallel --progress --jobs $CPU --joblog $LOG_F.log < $JL_FILE
tput bold
echo -ne "\nDone! - Run logged to $LOG_F.log\n"
tput sgr0
