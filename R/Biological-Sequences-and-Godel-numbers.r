# clear
cat("\014")
rm(list = ls())
dev.off(dev.list()["RStudioGD"])

# Libraries
library(seqinr)
library(MASS)
library(ggplot2)
library(ggpubr)
library(stringr)
library(combinat)

#--------------------- SELECT DATA ---------------------------------------------
# 'artificial' and 'artificial-non-uni' for artificial dataset, 'real' for real-world data (biom_data_150bp.rda)
dataset_type <- "real"

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
  
  primes = log2(primes)
  
  chars.enc = chars.enc %*% diag(primes)
  godel_nums <- rowSums(chars.enc) # gn_new
  return(godel_nums)
  
  
}


#----------------------------- Main code ---------------------------------------
primes <- sieve(20000)            # length of primes should be >= max sequence length
logOutput <- FALSE                # debug output
replicate <- TRUE                 # if we need to set specific seed numbers
type <- "DNA"                     # type of data

# encodings
encodings <- permn(c(1,2,3,4))
numberOfEncodings <- 24;              # number of different assignments of letters for Godel numbers (we have DNA data)

# The following values should execute within 6 minutes with a good enough resolution
seqLengthLimit <- 150
numberOfArtificialSeqs <- 90000     # "A","C","G","T"

# Creating artificial data set
# createRandomSequencesBasedOnDistr(numberOfArtificialSeqs, seqLengthLimit, c(0.3,0.3,0.2,0.2), "data/artificialSeq-non-unis.fasta")

# Loading a fasta file - artificial dataset
if (dataset_type == "artificial"){
  ompGene.list <- seqList("data/artificialSeqs.fasta", type)
  ompGene.list.raw <- ompGene.list
} else if (dataset_type == "real"){
  ompGene.list <- seqList("data/realSeqs_biom_150bp.fasta", type)
} else if (dataset_type == "artificial-non-uni"){
  ompGene.list <- seqList("data/artificialSeq-non-unis.fasta", type)
} else {
  print("Incorrect dataset type.")
}
sequences <- unlist(ompGene.list)

#encodingValues <- createRandomSequenceValues(seedValuesList, type) - we don't neet this function
encodingValues <- matrix(0, nrow = length(encodings), ncol = 4)
for (i in 1:nrow(encodingValues)){
  encodingValues[i,] <- encodings[[i]]
}


sizeExp <- length(ompGene.list)
selectedSequences <- c(1:sizeExp)

godelValuePoints <- data.frame(matrix(0, ncol = numberOfEncodings, nrow = sizeExp))
namesList <- list()
for (i in 1:numberOfEncodings) {
  namesList[i] <- paste('godel_log_pos', i, sep = '')
}
godelValuePoints <- setNames(godelValuePoints, namesList)
godelValuePoints$seqNames <- unlist(attributes(ompGene.list)$name)

stringAssignValues <- vector(mode="character", length=numberOfEncodings)

len_of_sequences <- str_length(sequences[1])
primes <- primes[1:len_of_sequences]
primes <- as.numeric(primes)

# main loop - godel number calculation
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
P1 <- sum(log2(primes[1:seqLengthLimit]))
P2 <- sum((log2(primes[1:seqLengthLimit]))^2)
theoreticalMeanEqual <- P1*2.5
theoreticalStdEqual <- sqrt(P2*1.25)
theoretical_dist <- rnorm(90000, theoreticalMeanEqual, theoreticalStdEqual)

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
indexPos = 1 # Specify which encoding to plot
binwidthPlot = 5
this_enc <- paste(encodingValues[indexPos,], collapse = '')

# comparisons file save
filepath <- paste('plots/comparisons_',dataset_type, '.csv', sep = '' )
write.csv(as.data.frame(statsPos), file = filepath)

this_title= expression(paste('Histogram and Density plot of ','log'[2],'(Godel numbers) - encoding 1234', sep = '' ))

# Plots
filepath <- paste('plots/histogram_encoding_', this_enc,'_',dataset_type, '.png', sep = '' )
png(file = filepath, width=1200, height=800)
my_plot <- ggplot(godelValuePoints, aes(x=godelValuePoints[, indexPos]), environment = environment()) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=binwidthPlot)+
  geom_density(alpha=.2, fill="#FF6666") +
  stat_function(fun = dnorm, args = list(mean = theoreticalMeanEqual, sd = theoreticalStdEqual), color = "darkred", size = 2, linetype = "dotdash") +
  labs(title = this_title, x=expression('log'[2]('Godel Numbers')), y = "Density")+
  xlim(min(godelValuePoints,theoretical_dist), max(godelValuePoints,theoretical_dist)) +
  scale_color_brewer(palette="Accent") + 
  theme(title  = element_text(size=20),
        axis.text.x = element_text(size = 20, face = 'bold'), axis.text.y = element_text(size = 15), 
        axis.title.x =  element_text(size=20)) 
print(my_plot)

dev.off()

# ---------------- Other plots -------------------------------------------------

p1 <- ggplot(statsPos, aes(x = reorder(stringAssignValues, fit_estimate), y = fit_estimate, group=1)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title= "Difference between theoretical and observed value of Mean (%)",x="Assignments", y = "Difference in mean (%)")

p2 <- ggplot(statsPos, aes(x = reorder(stringAssignValues, fit_estimate), y = fit_sd, group=1)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title= "Difference between theoretical and observed value of standard deviation (%)",x="Assignments", y = "Difference in  standard deviation (%)")

ggarrange(p1 + rremove("x.text"), p2,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)


# a t g c
lapply(which(statsPos$maxG == min(statsPos$maxG)),
       function(x) assignSets(encodingValues[x,], "DNA"))
# artificial dataset: 4 , 1 , 3 , 2 (min of meanG)
# artificial dataset: 4 , 2 , 1 , 3 (min of minG)
# artificial dataset: 2 , 1 , 3 , 4 (min of maxG)

# a t g c
lapply(which(statsPos$maxG == max(statsPos$maxG)),
       function(x) assignSets(encodingValues[x,], "DNA"))
# artificial dataset: 1 , 4 , 2 , 3 (max of meanG)
# artificial dataset: 4 , 3 , 2 , 1 (max of minG)
# artificial dataset: 1 , 3 , 4 , 2 (max of maxG)


