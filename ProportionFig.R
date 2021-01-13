rm(list=ls())

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

dat <- read.table("./coligoTable", stringsAsFactors = F, header = T)
head(dat)

its <- read.csv("./ITSmergedstats.csv", stringsAsFactors = F, header = F)
#reformat to faciltiate searching
its$short <- gsub(".*(ITS\\.[ATCG]*\\.[ATCG]*).*","\\1",its$V1)
its$short <- gsub("\\.","_",its$short)

#reorder
its <- its[order(its$short),]

sixteenS <- read.csv("./16Smergedstats.csv", stringsAsFactors = F, header = F)
#reformat to faciltiate searching
sixteenS$short <- gsub(".*(16S\\.[ATCG]*\\.[ATCG]*).*","\\1",sixteenS$V1)
sixteenS$short <- gsub("\\.","_",sixteenS$short)

#reorder
sixteenS <- sixteenS[order(sixteenS$short),]

#now clean up the names of dat and reorder to match
names(dat) <- gsub("^X","",names(dat) )

dat <- dat[,order(names(dat))]
dat <- dat[,-length(dat)]
names(dat) <- gsub("(16S\\_[ATCG]*\\_[ATCG]*).*","\\1",names(dat))
names(dat) <- gsub("(ITS\\_[ATCG]*\\_[ATCG]*).*","\\1",names(dat))


tail(names(dat))
tail(its$short)

sixteenS <- sixteenS[ sixteenS$short %in% names(dat),]
sixteenS <- sixteenS[!is.na(sixteenS$short),]
its <- its[its$short %in% names(dat),]
its <- its[!is.na(its$short),]

#ensure order is correct
dat <- dat[,match(names(dat), c(sixteenS$short, its$short))]
colSums(dat) /c(sixteenS$V4, its$V4)

ISD <- dat[length(dat[,1]),]

#Doesn't exactly match with manual calculations, rounding error maybe
pdf(width = 5, height = 5, file = "ColigoISDProportions1.pdf")

hist(log(as.numeric(ISD[1:length(sixteenS$V4)]/sixteenS$V4)),xlab = "Proportion (log10 scale)", main = "",
     ylim = c(0,1400), xlim = c(-14,0),col = add.alpha("orange",0.5),yaxt = "n", xaxt ="n")
axis(side = 2, at = c(seq(0,1400, by = 200)), label = c(seq(0,1400, by = 200)), las = 2)
axis(side = 1, at = c(seq(-14,0, by = 2)), label =c(seq(-14,0, by = 2)))

hist(log(colSums(dat)[1:length(sixteenS$V4)]/sixteenS$V4),xlab = "", ylab = "", main = "",axes = F, 
     col = add.alpha("aquamarine",0.5), add = T)
legend(c("ISD", "Coligos"), col = c(add.alpha("orange",0.5), add.alpha("aquamarine",0.5)), 
       x="topleft", cex = 2,
       pch = 15, bty = "n")
dev.off()


pdf(width = 5, height = 5, file = "ColigoISDProportions2.pdf")

hist(log(as.numeric(ISD[(length(sixteenS$V4)+1):length(ISD)]/its$V4)),xlab = "Proportion (log10 scale)", main = "",
     ylim = c(0,2500), xlim = c(-14,0),col = add.alpha("orange",0.5), yaxt = "n", xaxt ="n")
axis(side = 2, at = c(seq(0,2500, by = 500)), label = c(seq(0,2500, by = 500)), las = 2)
axis(side = 1, at = c(seq(-14,0, by = 2)), label =c(seq(-14,0, by = 2)))

hist(log(colSums(dat)[(length(sixteenS$V4)+1):length(dat[1,])]/its$V4),xlab = "", ylab = "", main = "",axes = F, 
     col = add.alpha("aquamarine",0.5), add = T)
dev.off()
