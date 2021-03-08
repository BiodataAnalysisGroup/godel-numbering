# clear
cat("\014")
rm(list = ls())
dev.off(dev.list()["RStudioGD"])

# Libraries
library(seqinr)
library(MASS)
library(ggplot2)
library(ggpubr)

load('data/biom_data_150bp.rda')

sequences <- biom_data_150bp[1]$sequence

# Creates random sequences
seq_list_to_fasta<- function(seq_list, fileNameRandSeqs) {
  
  count <- length(seq_list)
  
  sink(fileNameRandSeqs)
  
  for (i in 1:count) {
    cat(">Seq", i, "\n", sep = "")
    seqX <- seq_list[i]
    cat(paste(seqX,collapse=""), "\n", sep = "")
  }
  sink()
}

# Save as fasta_file
seq_list_to_fasta(sequences, 'data/realSeqs_biom_150bp.fasta')