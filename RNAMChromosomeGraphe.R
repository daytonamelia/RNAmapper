### RNAMChromosomeGraphe.R ###
## This script visualizes linked regions in a single chromosome in conjunction with the RNAmapper.py script.
## This script requires 3 inputs:
## argument 1: the atMarkers.vcf file for the chromosome of interest
## argument 2: the name of the stats file to append linkage regions
## argument 3: the name for the output plot (will have .jpg automatically appended to end)
## There is an optional fourth argument to change the linkageRatio (default = 0.98)


# Variables
# Get current directory, read in arguments
currentDir <- getwd()
args <- commandArgs(T)
if (length(args) > 4) stop('Error in RNAMChromosomeGraphe.R: Too many arguments.');
mutMarker_path <- args[1]
statsFile <- args[2]
plotOut <- args[3]
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

#Plot with INDELS!
jpeg(paste(plotOut,"_indels.jpg", sep=""), width=1000, bg="white")
options(scipen=999)
plot(x=mutMarker$POS, y=mutMarker$AVERAGE,
     pch=16, col="red", cex=2,
     ylim=c(0.5,1),
     main = plotOut,
     xlab="Position (bp)",
     ylab="Frequency")
points(x=mutMarker$POS, y=mutMarker$HIGHALLELE, col="black", pch=20)        # plot raw frequency
points(x=INDELS$POS, y=INDELS$HIGHALLELE, col="blue", pch=20, type="b", lwd=2)      # plot indels
abline(v = Ledge, col = "blue")                                 # vertical line at Ledge of homozygosity
abline(v = Redge, col = "blue")                                 #  vertical line at Redge of homozygosity
abline(h = 1, col = "black", lty = 3)                           # dashed line at freq=1.0, i.e. homozygosity
mtext(edges, 1, col = "blue")                                   # write position of above lines
dev.off()

#Plot without INDELS!
jpeg(paste(plotOut,".jpg", sep=""), width=1000, bg="white")
options(scipen=999)
plot(x=noINDELS$POS, y=noINDELS$AVERAGE,
     pch=16, col="red", cex=2,
     ylim=c(0.5,1),
     main = plotOut,
     xlab="Position (bp)",
     ylab="Frequency")
points(x=noINDELS$POS, y=noINDELS$HIGHALLELE, col="black", pch=20)        # plot raw frequency
abline(v = Ledge, col = "blue")                                 # vertical line at Ledge of homozygosity
abline(v = Redge, col = "blue")                                 #  vertical line at Redge of homozygosity
abline(h = 1, col = "black", lty = 3)                           # dashed line at freq=1.0, i.e. homozygosity
mtext(edges, 1, col = "blue")                                   # write position of above lines
dev.off()