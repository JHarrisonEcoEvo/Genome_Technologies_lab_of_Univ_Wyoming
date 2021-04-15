#Determine the extent to which coligos occur where they are not supposed to. 

dat <- read.table("./coligoTable_testRun_iseq10", stringsAsFactors = F, header = T)
dat[1:5,1:5]
dim(dat)

dat <- dat[dat$OTUID != "ISD",]
#modify dat$OTUID so that it matches the names of dat
dat$OTUID <- gsub("_", "", dat$OTUID)

expectedHit <- NA
totalReads <- NA
propExp <- NA
coligo <- NA
for(i in 14:length(names(dat))){
    ptrn <- gsub("c(oligo[A-H]\\d+)_[ABC]","\\1", names(dat)[i])
    #modify pattern so that it doesn't get A10 when it is A1
    ptrn <- paste(ptrn, "$", sep ="")
    expectedHit[i] <- dat[grep(ptrn, dat$OTUID), i]
    totalReads[i] <- sum(dat[, i])
    propExp[i] <- expectedHit[i] / totalReads[i]
    coligo[i] <- names(dat)[i]
}

#Turn into a matrix, or otherwise label by coligo.

propdf <- cbind(propExp, 
                coligo, 
                gsub("coligo([A-H])\\d+_[ABC]", "\\1", coligo),
                gsub("coligo[A-H](\\d+)_[ABC]", "\\1", coligo))
propdf <- propdf[order(propdf[,3], as.numeric(propdf[,4])),]


as <- propdf[grep("coligo([A-H])\\d+_A", propdf[,2]),]
bs <- propdf[grep("coligo([A-H])\\d+_B", propdf[,2]),]
cs <- propdf[grep("coligo([A-H])\\d+_C", propdf[,2]),]

matA <- matrix(nrow = 8, ncol = 12)
matB <- matrix(nrow = 8, ncol = 12)
matC <- matrix(nrow = 8, ncol = 12)

k <- 1
for(i in c("A", "B", "C", "D", "E", "F", "G", "H")){
  matA[k,] <- round(as.numeric(as[as[,3] == i, 1]), 3)
  matB[k,] <- round(as.numeric(bs[bs[,3] == i, 1]), 3)
  matC[k,] <- round(as.numeric(cs[cs[,3] == i, 1]), 3)
  k <- k + 1
}

row.names(matA) <- c("A", "B", "C", "D", "E", "F", "G", "H")
row.names(matB) <- c("A", "B", "C", "D", "E", "F", "G", "H")
row.names(matC) <- c("A", "B", "C", "D", "E", "F", "G", "H")

colnames(matA) <- seq(1,12, 1)
colnames(matB) <- seq(1,12, 1)
colnames(matC) <- seq(1,12, 1)

xtable::xtable(matA, digits = 3)
xtable::xtable(matB, digits = 3)
xtable::xtable(matC, digits = 3)


