

listfiles  <- c() # creates a vector
viterbi <- list() # creates a list
viterbi[1] <- "empty"
for (i in 1:22){
  listfiles[i] <- paste("/storage/atkinson/home/u242335/lai/bhrcsplit/BHRC_RFMIX_", i, ".Viterbi.txt", sep='')
}

for (k in 1:length(listfiles)){
  viterbi[[k]] <- read.table(listfiles[k], h=F, sep = " ")
  viterbi[[k]] <- viterbi[[k]][ , colSums(is.na(viterbi[[k]])) < nrow(viterbi[[k]])]
}


lastline <- c()
for (k in 1:length(viterbi)){
  lastline[k] <- nrow(viterbi[[k]])
}

for (k in 1:length(viterbi)){
  viterbi[[k]]$epos <- c(viterbi[[k]]$V2[2:lastline[k]] - 1, viterbi[[k]]$V2[lastline[k]] + 1)
}

lastcol <- c()
for (k in 1:length(viterbi)){
  lastcol[k] <- ncol(viterbi[[k]])
}

viterbi[[22]][c(1:10),c((lastcol[22]-10):lastcol[22])]

for (k in 1:length(viterbi)){
viterbi[[k]]$sgpos <- 0
viterbi[[k]]$egpos <- 0
viterbi[[k]]$`n snps` <- 1
}

viterbi_msp <- list() # creates a list
viterbi_msp[1] <- "empty"
for (k in 1:length(viterbi)){
viterbi_msp[[k]] <- viterbi[[k]][,c(1,2,lastcol[k],(lastcol[k]+1),(lastcol[k]+2),(lastcol[k]+3), 557:(lastcol[k]-1))]
}

viterbi_msp[[22]][c(1:10),c(1:10)]

lastcol_msp <- c()
for (k in 1:length(viterbi)){
  lastcol_msp[k] <- ncol(viterbi_msp[[k]])
}

viterbi_msp[[22]][c(1:10),c((lastcol_msp[k]-3):lastcol_msp[k])]

library(readr)
read_lines("/storage/atkinson/home/u242335/lai/BHRC.notref") -> ids
ids_bhrc <- rep(ids,each=2)

for (k in 1:length(viterbi)){
colnames(viterbi_msp[[k]]) <- c("#chm", "spos", "epos", "sgpos", "egpos","n snps", ids_bhrc)
}


viterbi_msp[[22]][c(1:10),c(1:15)]
viterbi_msp[[k]][c(1:10),c((lastcol_msp[k]-7):lastcol_msp[k])]


for (k in 1:22){
  write.table(viterbi_msp[[k]], paste("BHRC_RFMIX_chr",k, ".msp.tsv", sep=""), col.names = T, row.names = F, quote = F, sep = "\t")
}
