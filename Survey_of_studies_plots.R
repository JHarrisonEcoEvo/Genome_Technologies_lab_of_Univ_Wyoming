dat <- read.csv("Survey_microbiome_methods.csv")
head(dat)

table(dat$robot_TF)
table(dat$sequencingPlatform)
table(dat$pcr_replication)
table(dat$pcrVolume)
table(dat$ISD)
table(dat$cross.contam.checking)
table(dat$positiveControl)
table(dat$blank)
table(dat$multiplexingStrategy)

#Make PCR replication bar plot

pdf(width = 6, height = 6, file = "./visuals/pcr_rep.pdf")
barplot(height = c(1,4,2,10,2,30), las = 2, border = NA, ylab = "Num. of studies")
text(srt = 60, 
     x = c(0.5, 1.75, 2.9, 4.1, 5.5, 6.6), 
     y = -4, c("PCR-free", 
               "Singlet",
               "Duplicate",
               "Triplicate",
               "Quadruplicate",
               "Unclear"), xpd = NA )
dev.off()

#Make pcr volume plot
pdf(width = 6, height = 6, file = "./visuals/pcr_vol.pdf")
barplot(height = c(2,8,2,5,32), las = 2, border = NA, ylab = "Num. of studies")
text(srt = 60, 
     x = c(0.5, 1.75, 2.9, 4.1, 5.5), 
     y = -3, c(expression(paste("20",mu,"l", sep="")),
               expression(paste("25",mu,"l", sep="")),
               expression(paste("30",mu,"l", sep="")),
               expression("">=50),
               "Unclear"), xpd = NA )
dev.off()


#make sideways bar plots, Likert plots

# library
library(likert) 

# make new dataframe to edit for the plot
newdat <- dat[,c(10,9,7,11)]
for(i in 1:length(newdat)){
  newdat[,i] <- as.character(newdat[,i] )
}

newdat$ISD[newdat$ISD == "no, but qPCR"] <- "no"
newdat$blank[is.na(newdat$blank)] <- "no"
newdat$blank[newdat$blank=="unclear"] <- "no"

newdat$multiplexingStrategy[
  newdat$multiplexingStrategy != "unclear"] <- "yes"
newdat$multiplexingStrategy[
  newdat$multiplexingStrategy == "unclear"] <- "no"
names(newdat)[names(newdat) == "multiplexingStrategy"] <- "Library prep. ethods clear"
newdat <- newdat[1:50,]
for(i in 1:length(newdat)){
  newdat[,i] <- as.factor(newdat[,i] )
}

p <- likert(newdat)
pdf(width = 6, height = 6, 
    file = "./visuals/likert.pdf")
plot(p)
dev.off()

#make a scatter plot of how many people use controls for each of the years from
#2015-2020

dat[sample(grep("2019", dat[,1]), 10),]

dat[grep("2020", dat[,1]),]

years <- seq(2015, 2020, by = 1)
prop_Pos_yes <- c(20,50,20,20,40,10)
prop_Neg_yes <- c(20,50,40,30,30,30)

pdf(file = "./visuals/controlPlot.pdf", width = 6, height = 6)
plot(ylim = c(1,100),
     xlim = c(1,6),
     frame.plot = F,
     xlab = "",
     yaxt = "n",
     xaxt = "n",
     ylab = "",
      NULL)
axis(side = 2, lab = seq(0, 100, 20), at = seq(0, 100, 20), las = 2)
points(y = prop_Neg_yes, x = c(1:6), pch = 16, cex = 2)
points(y = prop_Pos_yes, x = c(1:6), pch = 16, col = "red", cex = 2)
axis(side = 1, lab = seq(2015, 2020, 1), 
     at =seq(1, 6, 1), las = 2, las = 1)

legend("topleft", legend = c("Negative control", 
                  "Positive control"),
       col = c("orange","blue"), pch = 19,
       bty = "n")
dev.off()
