#!/usr/bin/env Rscript

library(optparse)
library(digest)

option_list <- list(
    make_option(c("-t", "--trainseqs"), default="mm9.fa",
                help="Full path to a fasta format genome"),
    make_option(c("-f", "--foreground"), default="GSM1064678_mm_testes_nmi.bed",
                help="A bed file with the locations of NMIs"),
    make_option(c("-b", "--background"), default="not-nmis.bed",
                help="A bed file with the locations of background sequences"),
    make_option(c("-c", "--c"), default=1.0,
                help="SVM soft margin penalty"),
    make_option(c("-k", "--k"), default=2,
                help="K-mer length"),
    make_option(c("-w", "--window"), default=750,
                help="Split the sequences into windows of this many bases"),
    make_option(c("-p", "--predictseqs"), default="hg19.fa",
                help="Fasta file to do predictions on"),
    make_option(c("-j", "--cores"), default=1,
                help="Number of processes to use"),
    make_option(c("-d", "--seed"), default=531,
                help="Set the random seed at the beginning of the script")
    )
conf <- parse_args(OptionParser(option_list = option_list))
set.seed(conf$seed)

options(cores = conf$cores)

library(Rsamtools)
library(GenomicRanges)
library(shogun)
library(rtracklayer)
library(seqinr)

# --------------------------------------------------------------------------
# TRAIN
# --------------------------------------------------------------------------

source("../fit-spectrum-svm.R")

# Read in training sequences
genome <- FaFile(conf$trainseqs)
nmis <- import.bed(conf$foreground)
bg <- import.bed(conf$background)

# Trim off the ends that don't fit into a full window
width(nmis) <- width(nmis) - (width(nmis) %% conf$window)
nmis <- nmis[width(nmis) > 0]
width(bg) <- width(bg) - (width(bg) %% conf$window)
bg <- bg[width(bg) > 0]
# Bin the regions into windows
nmis_binned <- unlist(tile(nmis, width=conf$window))
bg_binned <- unlist(tile(bg, width=conf$window))

x <- toupper(c(Rsamtools::scanFa(genome, nmis_binned),
               Rsamtools::scanFa(genome, bg_binned)))

y <- factor(c(rep("NMI", length(nmis_binned)),
              rep("Background", length(bg_binned))), levels=c("NMI", "Background"))

stopifnot(length(x) == length(y))

# Fit model
params <- list()
params$C <- conf$c
params$k <- conf$k
# Extra parameters that we will hard code for now
params$alphabet <- "SNP"
params$usesign <- FALSE
params$eps <- 1e-4

cat("Fitting model (this may take a while)... ")
fit <- fit_spectrum_svm(matrix(x, ncol=1), y, params)
cat("done!\n")

# Save model for use with predict.R
saveRDS(fit, "fit.rds")

# Clean up a bit
rm(x, y)

# --------------------------------------------------------------------------
# PREDICT
# --------------------------------------------------------------------------

source("../predict-spectrum-svm.R")

cat("Loading sequences to predict on (can be slow when using an entire genome)... ")
pred_seqs <- seqinr::read.fasta(conf$predictseqs, as.string=TRUE)
cat("done!\n")

# Predict new windows
cat("Predicting NMIs... ")
preds <- predict_spectrum_svm(fit, unlist(pred_seqs), conf$window, param=params)
cat("done!\n")

saveRDS(preds, "preds.rds")

# TODO: generate a csv/bigWig file with the predictions, and a .bed file with
# windows > 0 score

