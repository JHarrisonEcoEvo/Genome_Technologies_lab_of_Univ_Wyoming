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


