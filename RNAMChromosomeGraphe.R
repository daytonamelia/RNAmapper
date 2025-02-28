### RNAMChromosomeGraphe.R ###
## This script visualizes linked regions in a single chromosome in conjunction with the RNAmapper.py script.
## This script requires 3 inputs:
## argument 1: the path to the atMarkers.vcf file for the chromosome of interest
## argument 2: the path of the stats file to append linkage regions
## argument 3: the path to the output plot (will have chromosome number and .jpg automatically appended to end)
## argument 4: optional change to linkageRatio (default = 0.98)


# Variables
# Get current directory, read in arguments
currentDir <- getwd()
args <- commandArgs(T)
if (length(args) > 4) stop('Error in RNAMChromosomeGraphe.R: Too many arguments.');
mutMarker_path <- args[1]
statsFile <- args[2]
plotOut <- args[3] # Path to plot with plot prefix name as file name (dont include chr or file extension!)
linkedRatio <- if (length(args) < 4) 0.98 else args[4]

# Load in files
mutMarker <- read.delim(mutMarker_path, header=F)
names(mutMarker) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "PL", "FORMAT",
                      "DP", "FREF", "FALT", "RREF", "RALT", "TOTALREF", "TOTALALT",
                      "REFRATIO", "ALTRATIO", "HIGHALLELE", "INDEL", "AVERAGE", "REMOVE")
mutMarker <- subset(mutMarker, select = -REMOVE) # Remove extra column

# Find regions of linked ratio
mutLinked <- mutMarker[mutMarker$AVERAGE >= linkedRatio,]
# get the position of the edges
Ledge <- 0
Ledge <- as.numeric(as.character(mutLinked[1,2]))
Redge <- 0
Redge <- as.numeric(as.character(mutLinked[nrow(mutLinked),2]))
edges <- paste(Ledge,Redge, sep="--")
linkSize <- 0
linkSize <- Redge-Ledge
# append linkage region to stats file
write(paste("Left linkage edge", Ledge),file=statsFile,append=TRUE)
write(paste("Right linkage edge", Redge),file=statsFile,append=TRUE)
write(paste("Linkage size", linkSize),file=statsFile,append=TRUE)

INDELS <- mutMarker[mutMarker$INDEL == "True",]
noINDELS <- mutMarker[mutMarker$INDEL == "False",]
plotName <- head(strsplit(tail(strsplit(plotOut, "/")[[1]], n=1), ".jpg")[[1]], n=1)[[1]]
chrName <- tail(strsplit(tail(strsplit(tail(strsplit(mutMarker_path, "/")[[1]], n=1), "_")[[1]], n=2)[[1]], "chr")[[1]], n=1)

#Plot with INDELS!
jpeg(paste(plotOut, "_chr", chrName, "_indels.jpg", sep=''), width=1000, bg="white")
options(scipen=999)
plot(x=mutMarker$POS, y=mutMarker$AVERAGE,
     pch=16, 
     col="red", 
     cex=4,
     cex.main=2,
     cex.axis=1.1,
     ylim=c(0.5,1),
     xlab="",
     ylab="",
     las=1,
     main = paste(plotName, "Marker Frequency Across Chromosome",chrName, sep=" "))
title(xlab="Position (bp)", ylab="Marker Frequency", cex.lab=1.5, line=2.5)
points(x=mutMarker$POS, y=mutMarker$HIGHALLELE, col="black", pch=20)        # plot raw frequency
points(x=INDELS$POS, y=INDELS$HIGHALLELE, col="blue", pch=20, type="b", lwd=2)      # plot indels
abline(v = Ledge, col = "blue", lwd=2.5)                                 # vertical line at Ledge of homozygosity
abline(v = Redge, col = "blue", lwd=2.5)                                 #  vertical line at Redge of homozygosity
abline(h = 1, col = "black", lty = 3, lwd=2)                           # dashed line at freq=1.0, i.e. homozygosity
mtext(edges, 1, col = "blue", cex=1.2, font=2)              # write position of above lines
dev.off()

#Plot without INDELS!
jpeg(paste(plotOut, "_chr", chrName, ".jpg", sep=''), width=1000, bg="white")
options(scipen=999)
plot(x=noINDELS$POS, y=noINDELS$AVERAGE,
     pch=16, 
     col="red", 
     cex=4,
     cex.main=2,
     cex.axis=1.1,
     ylim=c(0.5,1),
     xlab="",
     ylab="",
     las=1,
     main = paste(plotName, "Marker Frequency Across Chromosome",chrName, sep=" "))
title(xlab="Position (bp)", ylab="Marker Frequency", cex.lab=1.5, line=2.5)
points(x=noINDELS$POS, y=noINDELS$HIGHALLELE, col="black", pch=20)        # plot raw frequency
abline(v = Ledge, col = "blue", lwd=3)                                 # vertical line at Ledge of homozygosity
abline(v = Redge, col = "blue", lwd=3)                                 #  vertical line at Redge of homozygosity
abline(h = 1, col = "black", lty = 3, lwd=2)                           # dashed line at freq=1.0, i.e. homozygosity
mtext(edges, 1, col = "blue", cex=1.2, font=2)              # write position of above lines
dev.off()
