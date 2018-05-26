seqs <-
  c("CTGAACGAGG",
    "CATTCCTTGG",
    "CGGGCCGGAC",
    "GTGGACCCTC",
    "AAAATGTACC",
    "ATGCCGCCATCG")

# convert to the nucleotide position matrix
seq_mtx <- do.call(rbind, strsplit(seqs, ""))

# tabelize
seq_tbl <- split(seq_mtx, col(seq_mtx))

#count probability of G+C at the each position
unlist(lapply(seq_tbl, function(column) {
  g_prob = length(which(column == "G"))
  c_prob = length(which(column == "C"))
  likelihood = (g_prob + c_prob) / length(column)
}))

thA0 <- 0.4
thB0 <- 0.6


library("stringr")

#' Expectation step of the EM
#'
#'@return matrix with rows A - GC; A - non-GC; the same for B
expect <- function(thA0, thB0, data) {
  E = matrix(ncol = length(data), nrow = 4)
  i = 1
  for (seq in data) {
  
    #probability that sequence is from A
    prob_A <- dbinom(GC_count(seq), nchar(seq), thA0)
    prob_B <- dbinom(GC_count(seq), nchar(seq), thB0)
    
    #A - GC
    E[1, i] = GC_count(seq) * prob_A / (prob_A + prob_B)
    #A - non-GC
    E[2, i] = (nchar(seq) - GC_count(seq)) * prob_A / (prob_A + prob_B)
    
    #B - GC
    E[3, i] = GC_count(seq) * prob_B / (prob_A + prob_B)
    #B - non-GC
    E[4, i] = (nchar(seq) - GC_count(seq)) * prob_B / (prob_A + prob_B)
    
    i = i + 1
  }
  
  return(E)
}

#' Maximization step of the EM
#'
#'@param weighted_mtx matrix of the sequences weights received during E step
#'@return tupel of new ThetaA and ThetaB
maximize <- function(weighted_mtx) {
  # calculate sums of the 4 different counts (see E-step commentary)
  E <- rowSums(weighted_mtx) 
  new_thA <- E[1] / (E[1] + E[2]) # GC-score / sum of all scores in A
  new_thB <- E[3] / (E[3] + E[4])
  return(list(thA = round(new_thA, 3), thB = round(new_thB, 3)))
}

GC_count <- function(seq) str_count(seq, "G") + str_count(seq, "C")

thA <- thA0
thB <- thB0

mtx <- expect(thA, thB, seqs)
mtx

library("gsubfn")

new_th <- maximize(mtx)
thA <- new_th[[1]]
thB <- new_th[[2]]

thA
thB

#2nd run
mtx <- expect(thA, thB, seqs)
mtx

new_th <- maximize(mtx)
thA <- new_th[[1]]
thB <- new_th[[2]]

thA
thB

i = 1
for (seq in seqs) {
  prob_A <- dbinom(GC_count(seq), nchar(seq), thA)
  prob_B <- dbinom(GC_count(seq), nchar(seq), thB)
  
  
  weighted_prbA <- prob_A / (prob_A + prob_B)
  weighted_prbB <- prob_B / (prob_A + prob_B)
  
  if (weighted_prbA > weighted_prbB)
    message(paste(
      "With probability",
      toString(round(weighted_prbA *
                       100, 2)),
      "%, motif",
      toString(i),
      "is expected to be of A!"
    ))
  else
    message(paste(
      "With probability",
      toString(round(weighted_prbB * 100, 2)),
      "%, motif",
      toString(i),
      "is expected to be of B!"
    ))
  i = i + 1
}

#measure the speed of execution
benchmark("expect"= {expect(thA0, thB0, seqs)}, "maximize" = {maximize(mtx)})
