#Determine the extent to which coligos occur where they are not supposed to. 

dat <- read.table("./coligoTable", stringsAsFactors = F, header = T)
dat[1:5,1:5]
dim(dat)

key <- read.csv("./NovaSeq2_DemuxJH.csv", stringsAsFactors = F) #demultiplexing key

dim(key) #some stuff missing. No coligos?

names(dat) <- gsub("X(.*)","\\1",names(dat))

key$forward_barcode <- toupper(key$forward_barcode)
key$reverse_barcode <- toupper(key$reverse_barcode)

key$combo <- paste(key$locus, key$forward_barcode, key$reverse_barcode, key$samplename, sep = "_")
key$combo <- gsub("-","_",key$combo)
          
for(i in 1:length(key$combo)){
  names(dat)[grep(key$combo[i], names(dat))] <- key$wellposition[i]
}

#test to ensure the substitution worked
which(nchar(names(dat)) > 3)


dat$OTUID <- gsub("Coligo_", "",dat$OTUID)

dat <- dat[dat$OTUID != "ISD",]

percent_target <- NA

for(i in 2:length(names(dat))){
  if(dat[dat$OTUID == names(dat)[i],i] > 20){
    percent_target[i] <- dat[dat$OTUID == names(dat)[i],i] / sum(dat[,i])
  }
}
mean(na.omit(percent_target))
summary(percent_target)
table(percent_target > 0.90)
table(percent_target == 0)

pdf(width = 7, height = 5, file = "coligoProportions.pdf")
  hist(percent_target, main = "", xlab = "Proportion of coligo reads that were from expected sequences")
dev.off()

#dat[,1:15]

#Figure out how often contamination was from something adjacent. 
# and How often contamination was two wells or x wells away.

hits <- dat[dat[,3] > 0,c(1,3)]

#use this to determine hits and compare it to the target. Look up the combination of hit to target to find spatial distance. 
hits$OTUID

#create a distance matrix for wells in a plate.

#First, make a cartesian coordinate based matrix of points that represents the 12 x 8 plate with names.
columns <- seq(1,12, 1)
rows <- seq(1,8,1)
rowNames <- c("A","B","C","D","E","F","G","H")
k <- 1

point_coords <- matrix(nrow = 8*12, ncol = 3)
for(i in 1:12){
  for(j in 1:8){
    point_coords[k,1] <- i
    point_coords[k,2] <- j
    point_coords[k,3] <- paste(rowNames[j],i, sep = "")
  k <- k + 1
  }
}

#We can then extract coordinates for each well and calculate their Euclidean distances. 
#e.g., 

# point1 <- point_coords[point_coords[,3] == "A12",]
# point2 <- point_coords[point_coords[,3] == "A1",]
# 
# raster::pointDistance(as.numeric(point1[1:2]), as.numeric(point2[1:2]), lonlat = F)

#We now go through the dataframe and extract the non-target hits and calculate their distance from the target.
hits <- list()

for(i in 2:length(dat)){ #1 is a nonsense index bc is otuid vs. otuid
  hits[[i]] <- dat[dat[,i] > 0,c(1,i)]
}

counter <- 0
numberColigoContam <- vector()
k <- 0
distance <- rep(list(rep(list(),1)),length(dat))
for(i in 2:length(hits)){
#remove the target
  
nontarget <- hits[[i]][!names(hits[[i]])[2] == hits[[i]][,1],]
#calculate distances from target of the misses
numberContams <- 0
if(length(nontarget[,1]) > 0){
  
  for(j in 1:length(nontarget$OTUID)){
    miss <- point_coords[point_coords[,3] == nontarget$OTUID[j] ,]
    target <- point_coords[point_coords[,3] ==  names(hits[[i]])[2] ,]
    
    #sqrt((x-a)^2 + (y-b^2)
    distance[[i]][[j]] <- sqrt((as.numeric(miss[1:2])[1] - as.numeric(target[1:2])[1])^2 +
    (as.numeric(miss[1:2])[2] - as.numeric(target[1:2])[2])^2)
    numberContams <- numberContams + 1
  }
}else{
  counter <-  counter + 1
}
numberColigoContam[k] <- numberContams
k <- k + 1
}

table(unlist(distance) ==1)
summary(unlist(distance))
pdf(width = 8, height = 6, file = "ContamDist.pdf")
plot(sort(table(round(unlist(distance)))),
     ylab = "Incidence",
     xlab = "Euclidean distance (rounded to nearest integer)", las = 2, xaxt = "n")
axis(side = 1, at = seq(1,13,1), labels  = seq(1,13,1))
dev.off()


table(colSums(dat[,2:length(dat)])>0)
sum(table(numberColigoContam)[4:length(table(numberColigoContam))]) / sum(table(numberColigoContam))

