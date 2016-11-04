# Predict Non-methylated Genomic Regions using a Spectrum Kernel SVM

These scripts can be used to predict non-methylated regions of the genome using
a spectrum kernel support vector machine. The SVM must first be trained on a
set of known non-methylated regions.

## Software Dependencies

- [Shogun machine learning toolbox](http://shogun-toolbox.org/) (version 3.2.1) and its R library
- R and the following R packages:
    - optparse
    - digest
    - Rsamtools
    - GenomicRanges
    - rtracklayer
    - seqinr

## Usage

You can run the train and predict script with the `--help` argument in order to
get a list of all available options.

```bash
$ ./train-and-predict.R --help


Usage: ./train-and-predict.R [options]


Options:
        -t TRAINSEQS, --trainseqs=TRAINSEQS
                Full path to a fasta format genome

        -f FOREGROUND, --foreground=FOREGROUND
                A bed file with the locations of NMIs

        -b BACKGROUND, --background=BACKGROUND
                A bed file with the locations of background sequences

        -c C, --c=C
                SVM soft margin penalty

        -k K, --k=K
                K-mer length

        -w WINDOW, --window=WINDOW
                Split the sequences into windows of this many bases

        -p PREDICTSEQS, --predictseqs=PREDICTSEQS
                Fasta file to do predictions on

        -j CORES, --cores=CORES
                Number of processes to use

        -d SEED, --seed=SEED
                Set the random seed at the beginning of the script

        -h, --help
                Show this help message and exit

```

Example invocation:

```bash

$ ./train-and-predict.R -t ~/example/hg19.fa \
                        -f ~/example/nmis-chr1-small.bed \
			-b ~/example/notnmis-chr1-small.bed \
			-c 1.0 -k 3 -w 750 \
			-p ~/example/hg19-chr10_11.fa \
			-j 20
```
