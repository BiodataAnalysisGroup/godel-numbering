# clear
cat("\014")
rm(list = ls())
dev.off(dev.list()["RStudioGD"])

# libraries
library(microseq)
library(stringr)
library(seqinr)
library(MASS)
library(ggplot2)
library(ggpubr)
library(combinat)

#-----------------   SPECIFYING FILE -------------------------------------------
# folder where ERR2681749_1.fastq.gz and ERR2681763_1.fastq.gz are stored
emtab_data_folder <- paste(dirname(getwd()), 'E-MTAB-6962', sep = '/')

# files
file1 <- 'ERR2681749_1.fastq.gz' # S10_L001 - subset 6
file2 <- 'ERR2681763_1.fastq.gz' # S1_L001 - subset 8

# choose file
# file1: ERR2681749_1.fastq.gz - S10_L001
# file2: ERR2681763_1.fastq.gz - S1_L001
file_chosen <- file2


if (file_chosen == 'ERR2681749_1.fastq.gz'){
  file_code = 'S10_L001'
} else if (file_chosen == 'ERR2681763_1.fastq.gz'){
  file_code = 'S1_L001'
} else {
  print("error, wrong file")
}
filepath <- paste(emtab_data_folder, file_chosen, sep = '/')


#---------------- FUNCTIONS ----------------------------------------------------

# Returns prime numbers up to integer n
sieve <- function(n) {
  n <- as.integer(n)
  if(n > 1e6) stop("n too large")
  primes <- rep(TRUE, n)
  primes[1] <- FALSE
  last.prime <- 2L
  for(i in last.prime:floor(sqrt(n)))
  {
    primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
    last.prime <- last.prime + min(which(primes[(last.prime+1):n]))
  }
  which(primes)
}

# Reads sequence from input file
seqList <- function(x, type) {
  if (type == "DNA") {
    as.list(read.fasta(x,
                       as.string = TRUE,
                       seqtype="DNA",
                       seqonly = TRUE,
                       strip.desc = TRUE))
  }
  else {
    if (type == "AA") {
      as.list(read.fasta(x,
                         as.string = FALSE,
                         seqtype="AA",
                         seqonly = FALSE,
                         strip.desc = TRUE))
    }
    else {
      NULL
    }
  }
}

# Generates random DNA sequences
createRandomSequencesBasedOnDistr <- function(count, length, prob=c(0.25,0.25,0.25,0.25), fileNameRandSeqs) {
  
  sink(fileNameRandSeqs)
  for (i in 1:count) {
    cat(">Seq", i, "\n", sep = "")
    seqX <- sample(c("A","C","G","T"), length, rep=TRUE, prob)
    cat(paste(seqX,collapse=""), "\n", sep = "")
  }
  sink()
}

# Generating random permutations of godel number assignments
createRandomSequenceValues <- function(seedList, type) {
  
  if ( type == "DNA" ) {
    results <- matrix(, nrow = length(seedList), ncol = 4)
  }
  else if ( type == "AA" ) {
    results <- matrix(, nrow = length(seedList), ncol = 23)
  }
  
  for (i in 1:length(seedList)) {
    if ( type == "DNA" ) {
      results[i,] <- sample(seq(from = 1, to = 4, by = 1), size = 4, replace = FALSE)
    }
    else if ( type == "AA" ) {
      results[i,] <- sample(seq(from = 1, to = 23, by = 1), size = 23, replace = FALSE)
    }
  }
  return(results)
}

assignSets <- function(randSequenceValues, type) {
  
  stringValues <- ""
  
  if (type == "DNA") {
    #DNA Set
    cat("Sequence Values: (")
    
    for (i in 1:4) {
      stringValues <- paste(stringValues, randSequenceValues[i], ", ")
    }
    cat(stringValues)
    cat(")\n")
    assign("A", randSequenceValues[1], envir = .GlobalEnv)
    assign("C", randSequenceValues[2], envir = .GlobalEnv)
    assign("G", randSequenceValues[3], envir = .GlobalEnv)
    assign("T", randSequenceValues[4], envir = .GlobalEnv)
  }
  else if ( type == "AA") {
    # AA Set 1
    cat("Sequence Values: (")
    for (i in 1:23) {
      stringValues <- paste(stringValues, randSequenceValues[i], ", ")
    }
    cat(stringValues)
    cat(")\n")
    assign("A", randSequenceValues[1], envir = .GlobalEnv) # alanine, ala
    assign("R", randSequenceValues[2], envir = .GlobalEnv) # arginine, arg
    assign("N", randSequenceValues[3], envir = .GlobalEnv) # asparagine, asn
    assign("D", randSequenceValues[4], envir = .GlobalEnv) # aspartic acid, asp
    assign("B", randSequenceValues[5], envir = .GlobalEnv) # sparagine or aspartic acid, asx
    assign("C", randSequenceValues[6], envir = .GlobalEnv) # cysteine, cys
    assign("E", randSequenceValues[7], envir = .GlobalEnv) # glutamic acid, glu
    assign("Q", randSequenceValues[8], envir = .GlobalEnv) # glutamine, gln
    assign("Z", randSequenceValues[9], envir = .GlobalEnv) # glutamine or glutamic acid, glx
    assign("G", randSequenceValues[10], envir = .GlobalEnv) # glycine, gly
    assign("H", randSequenceValues[11], envir = .GlobalEnv) # histidine, his
    assign("I", randSequenceValues[12], envir = .GlobalEnv) # isoleucine, ile
    assign("L", randSequenceValues[13], envir = .GlobalEnv) # leucine, leu
    assign("K", randSequenceValues[14], envir = .GlobalEnv) # lysine, lys
    assign("M", randSequenceValues[15], envir = .GlobalEnv) # methionine, met
    assign("F", randSequenceValues[16], envir = .GlobalEnv) # phenylalanine, phe
    assign("P", randSequenceValues[17], envir = .GlobalEnv) # proline, pro
    assign("S", randSequenceValues[18], envir = .GlobalEnv) # serine, ser
    assign("T", randSequenceValues[19], envir = .GlobalEnv) # threonine, thr
    assign("W", randSequenceValues[20], envir = .GlobalEnv) # tryptophan, trp
    assign("Y", randSequenceValues[21], envir = .GlobalEnv) # tyrosine, tyr
    assign("V", randSequenceValues[22], envir = .GlobalEnv) # valine, val
    assign("X", randSequenceValues[23], envir = .GlobalEnv) # undetermined
  }
  
  stringValues
}

# Statistics function
godelStatistics <- function(x) {
  
  if (logOutput) {
    cat(paste("  Min  : ", summary(x)[1], "\n"))
    cat(paste("1st Qu.: ", summary(x)[2], "\n"))
    cat(paste("Median : ", summary(x)[3], "\n"))
    cat(paste(" Mean  : ", summary(x)[4], "\n"))
    cat(paste("3rd Qu.: ", summary(x)[5], "\n"))
    cat(paste("  Max  : ", summary(x)[6], "\n"))
    cat(paste("St. Dev: ", sd(x), "\n"))
  }
  
  results <- list(minG = summary(x)[1], firstQG = summary(x)[2],
                  medianG = summary(x)[3], meanG = summary(x)[4], 
                  thirdQG = summary(x)[5], maxG = summary(x)[6],
                  stdG = sd(x))
  
  return(results)
}

calculateGodelNumbers <- function(sequences, primes, encoding){
  
  
  chars = str_split(sequences, "", simplify = TRUE)
  
  chars.enc = matrix(data = 0, nrow = nrow(chars), ncol = ncol(chars))
  
  chars.enc = chars.enc + (chars == "A") * encoding[1]
  chars.enc = chars.enc + (chars == "C") * encoding[2]
  chars.enc = chars.enc + (chars == "G") * encoding[3]
  chars.enc = chars.enc + (chars == "T") * encoding[4]
  
  primes = log(primes)
  
  chars.enc = chars.enc %*% diag(primes)
  godel_nums <- rowSums(chars.enc) # gn_new
  return(godel_nums)
}

#----------------------------- Main code ---------------------------------------
# readfastq - checking
sequences <- readFastq(filepath)

# Lengths of sequences
lens <- apply(X = sequences[,2], FUN = str_length, MARGIN = 1)
unique(lens)
#table(lens) # sequences with 151 bases are the most frequent

# Keep only sequences of length 151
sequences <- sequences[which(lens == 151 ),]

#----------------------------- Analysis ----------------------------------------

primes <- sieve(20000)            # length of primes should be >= max sequence length
numberOfEncodings <- 24;          # number of different assignments of letters for Godel numbers (we have DNA data)
type <- "DNA"                     # type of data
logOutput <- FALSE                # debug output
encodings <- permn(c(1,2,3,4))
limit  <- 100000
seqLengthLimit <- 151
sequences <- sequences[1:limit,]

encodingValues <- matrix(0, nrow = length(encodings), ncol = 4)
for (i in 1:nrow(encodingValues)){
  encodingValues[i,] <- encodings[[i]]
}

seqnames <- sequences$Header
sequences <- sequences$Sequence

sizeExp <- length(sequences)
selectedSequences <- c(1:sizeExp)

godelValuePoints <- data.frame(matrix(0, ncol = numberOfEncodings, nrow = sizeExp))
namesList <- list()
for (i in 1:numberOfEncodings) {
  namesList[i] <- paste('godel_log_pos', i, sep = '')
}
godelValuePoints <- setNames(godelValuePoints, namesList)
godelValuePoints$seqNames <- seqnames

stringAssignValues <- vector(mode="character", length=numberOfEncodings)

len_of_sequences <- str_length(sequences[1])
primes <- primes[1:len_of_sequences]
primes <- as.numeric(primes)

# Main loop - godel number calculation
for (indexPos in 1:numberOfEncodings) {
  
  stringAssignValues[indexPos] <- assignSets(encodingValues[indexPos,], type)
  encoding <- encodingValues[indexPos,]
  godelValuePoints[, paste('godel_log_pos', indexPos, sep = '')] <- calculateGodelNumbers(sequences, primes, encoding)
  
}

statsPos <- data.frame(matrix(0, ncol = 7, nrow = numberOfEncodings))
statsPos <- setNames(statsPos, c("minG", "firstQG", "medianG", "meanG", "thirdQG", "maxG", "stdG"))
for (indexPos in 1:numberOfEncodings) {
  statsPos[indexPos, ] <- godelStatistics(godelValuePoints[, paste('godel_log_pos', indexPos, sep = '')])
}

statsPos_rownames <- c()
for (i in 1:numberOfEncodings){
  this_enc <- paste(encodingValues[i,], collapse = '')
  statsPos_rownames <- c(statsPos_rownames, paste("enc_", this_enc, sep = ''))
}

rownames(statsPos) <- statsPos_rownames


# Theoretical distribution parameters
P1 <- sum(log(primes[1:seqLengthLimit]))
P2 <- sum((log(primes[1:seqLengthLimit]))^2)
theoreticalMeanEqual <- P1*2.5
theoreticalStdEqual <- sqrt(P2*1.25)
theoretical_dist <- rnorm(limit, theoreticalMeanEqual, theoreticalStdEqual)

# comparisons
para <- matrix(nrow=numberOfEncodings, ncol = 2)
for (indexPos in 1:numberOfEncodings) {
  fit <- fitdistr(godelValuePoints[, indexPos], "normal")
  para[indexPos, 1] <- fit$estimate["mean"]
  para[indexPos, 2] <- fit$estimate["sd"]
}


statsPos$stringAssignValues <- as.factor(stringAssignValues)
statsPos$fit_estimate <- ((para[,1]-theoreticalMeanEqual)/para[,1])*100
statsPos$fit_sd <- ((para[,2]-theoreticalStdEqual)/para[,2])*100

statsPos$theoretical_mean <- theoreticalMeanEqual
statsPos$theoretical_std <-theoreticalStdEqual

# t-test
statsPos$p.value <- NaN
for (i in 1:numberOfEncodings){
  t_test_results <- t.test(godelValuePoints[,i], theoretical_dist)
  statsPos[i,]$p.value <- t_test_results$p.value
}


#-------------------- Main histogram plotting ----------------------------------
dir.create('plots')
indexPos = 1 #7,9,15,16,17,18,19,21
binwidthPlot = 5
this_enc <- paste(encodingValues[indexPos,], collapse = '')

# comparisons file save
filepath <- paste('plots/comparisons_human_',file_code, '.csv', sep = '' )
write.csv(as.data.frame(statsPos), file = filepath)

# Plots
filepath <- paste('plots/histogram_human_', file_code,'_encoding_', this_enc,'.png', sep = '' )
png(file = filepath, width=800, height=600)
my_plot <- ggplot(godelValuePoints, aes(x=godelValuePoints[, indexPos]), environment = environment()) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=binwidthPlot)+
  geom_density(alpha=.2, fill="#FF6666") +
  stat_function(fun = dnorm, args = list(mean = para[indexPos, 1], sd = para[indexPos, 2]), color = "darkred", size = 2, linetype = "dotdash") +
  labs(title= paste("Histogram and Density plot of Godel numbers - encoding ", this_enc, sep = ''),x="Godel numbers", y = "Density")+
  scale_color_brewer(palette="Accent") + 
  theme_minimal()
print(my_plot)
dev.off()



