## MinION utility scripts

A collection of various scripts for common minION related tasks.

### For processing barcoded reads using GNU parallel:

* `alignBarcodedFATQ_parallel.sh` - takes basecalled output of Albacore and returns a sorted single index BAM and gzipped FASTQ for each barcode. Uses: minimap2, samtools and pigz, NOTE: requires you to use the `--barcode` and `-o fastq` with `read_fast5_basecaller.py` also this will delete your original uncompressed `.fastq` files after use.

### SV detection with multiple input BAM and GNU parallel

* `alignNGM-LR_parallel.sh` aligns all `.fastq.gz` in a directory produced by above using the [NGM-LR](https://github.com/philres/ngmlr) aligner.  NGM-LR bam useful for running Sniffles uses samtools 1.7 to generate sorted index BAM.

* `Sniffles_parallel.sh` runs the [Sniffles](https://github.com/fritzsedlazeck/Sniffles) SV detector program on all BAM files generated using the above script, output will be saved to VCF_sniffles directory in the preceding parent directory.  _Note change parameters here to suite your own experimental design / SV detection issues_.

* `alignLAST_parallel.sh` runs the [LAST](http://last.cbrc.jp/) aligner on all `.fastq.gz` files in a directory requires [sambamba](http://lomereiter.github.io/sambamba/).  _Note_ `last-train` _is run to optimise alignment params on the first BAM file in the directory_.

* `NanoSV_parallel.sh` runs [NanoSV](https://github.com/mroosmalen/nanosv) SV detector on all BAM files generated using the above script, output will be saved to a NanoSV_VCF directory in the preceding parent directory.  _Note change parameters here to suite your own experimental design / SV detection issues_.

* `alignMinimap2_parallel.sh` runs [Minimap2](https://github.com/lh3/minimap2) only, assuming you've already run the `alignBarcodedFATQ_parallel.sh` and perhaps want to try different alignment params or a different genome, input is a directory of `.fastq.gz` as with the `alignNGM-LR_parallel.sh` and `alignLAST_parallel.sh` script.
