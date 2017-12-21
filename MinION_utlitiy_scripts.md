## MinION utlitiy scripts

A collection of various scripts for common minION related tasks.

* `alignBarcodedFATQ_parallel.sh` - takes basecalled output of Albacore and returns a sorted single index BAM and gzipped FASTQ for each barcode. Uses: minimap2 2.6, samtools 1.6 and pigz, NOTE: requires you to use the `--barcode` and `-o fastq` with `read_fast5_basecaller.py` also this will delete your original uncompressed .fastq files after use.
* `alignNGM-LR_parallel.sh` aligns FATQ produced by above using the NGM-LR aligner.  NGM-LR bam useful for running Sniffles uses samtools 1.6 to generate sorted index BAM.