# Variables
# Get current directory, read in arguments
currentDir <- getwd()
args <- commandArgs(T)
if (length(args) > 3) stop('Error in RNAMapperGraphe.R: Too many arguments.');
mutMarker_path <- args[1]
plotOut <- args[2]
linkedRatio <- if (length(args) < 3) 0.98 else args[3]

# Load in files
mutMarker <- read.delim(mutMarker_path, header=F)
names(mutMarker) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
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

INDELS <- mutMarker[mutMarker$INDEL == "True",]
noINDELS <- mutMarker[mutMarker$INDEL == "False",]

#Plot with INDELS!
jpeg(paste(plotOut,"_indels.jpg", sep=""), width=1000, bg="white")
options(scipen=999)
plot(mutMarker$POS, mutMarker$AVERAGE,
     pch=16, col="red", cex=2,
     ylim=c(0.5,1),
     main = plotOut,
     xlab="Position (bp)",
     ylab="Frequency")
points(mutMarker$POS, mutMarker$HIGHALLELE, col="black", pch=20)        # plot raw frequency
points(INDELS$POS, INDELS$HIGHALLELE, col="blue", pch=20)      # plot indels
abline(v = Ledge, col = "blue")                                 # vertical line at Ledge of homozygosity
abline(v = Redge, col = "blue")                                 #  vertical line at Redge of homozygosity
abline(h = 1, col = "black", lty = 3)                           # dashed line at freq=1.0, i.e. homozygosity
mtext(edges, 1, col = "blue")                                   # write position of above lines
dev.off()

#Plot without INDELS!
jpeg(paste(plotOut,".jpg", sep=""), width=1000, bg="white")
options(scipen=999)
plot(noINDELS$POS, noINDELS$AVERAGE,
     pch=16, col="red", cex=2,
     ylim=c(0.5,1),
     main = plotOut,
     xlab="Position (bp)",
     ylab="Frequency")
points(noINDELS$POS, noINDELS$HIGHALLELE, col="black", pch=20)        # plot raw frequency
abline(v = Ledge, col = "blue")                                 # vertical line at Ledge of homozygosity
abline(v = Redge, col = "blue")                                 #  vertical line at Redge of homozygosity
abline(h = 1, col = "black", lty = 3)                           # dashed line at freq=1.0, i.e. homozygosity
mtext(edges, 1, col = "blue")                                   # write position of above lines
dev.off()
