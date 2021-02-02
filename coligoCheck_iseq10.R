rm(list=ls())
dat <- read.table("./coligoTable_testRun_iseq10", stringsAsFactors = F, header = T)
dim(dat) #should be 96 long

#remove the synthgene only  and mock community only samples
dat <- dat[,c(1,grep("coligo",names(dat)))]
dat$OTUID <- gsub("^C","c", dat$OTUID)
dat$OTUID <- gsub("_","", dat$OTUID)
dat <- dat[dat$OTUID != "ISD",]

#Visualize technical variation
pdf(width = 11, height= 7, file = "coligoTechVar.pdf")

par(mar=c(6,6,0,3), mfrow = c(1,2), oma = c(2,2,2,2))
plot(xlim = c(0,96),
     ylim = c(0, 30000),
              NULL,
     ylab = "",
     xlab = "Coligo",
     cex.lab = 2,
     xaxt = "n",
     bty = "n",
     las = 2)
mtext(side = 2, "Read count", line = 4, cex = 2)
axis(side = 1, at = c(1,10,20,30,40,50,60,70,80,90,96), labels = c(1,10,20,30,40,50,60,70,80,90,96))

sums <- NA
contams <- NA
for(i in 1:length(dat$OTUID)){
  #subset to the correct replicates for each coligo and omit any of the minor contaminants
  df <- dat[i,grep(paste(dat$OTUID[i],"_[ABC]$",sep = ""), names(dat))]
  points(y = df, x = rep(i, length(df)), col = "darkorchid", pch = 16, type = "o")
  sums[i] <- sum(df)
  contams[i] <- sum(colSums(dat[-i,grep(paste(dat$OTUID[i],"_[ABC]$",sep = ""), names(dat))]))
  }

#Make a histogram
#pdf(7,7,file = "ColigoReadCounts.pdf")
hist(sums, 
     main = "", 
     las = 2,
     xaxt = "n",
     xlab = "Read count (summed)",
     cex.lab = 2,
     col = "darkolivegreen3")
axis(side = 1, at = seq(0,50000, by = 10000), labels = seq(0,50000, by = 10000))

text(x = -12000, 
     y = 35, 
     xpd = NA, 
     labels = "b)", 
     cex = 2)
text(x = -73000, 
     y = 35, 
     xpd = NA, 
     labels = "a)", 
     cex = 2)
dev.off()

options(scipen = 999)
sum(contams) /( sum(sums)+ sum(contams))

#determien proportion of wells with some contam
numContam <- NA
for(i in 2:length(dat)){
  numContam[i] <- length(which(dat[,i] > 0))
}
1 - (table(numContam)[1] / sum(table(numContam)))
