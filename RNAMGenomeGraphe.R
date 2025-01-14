### RNAMGenomeGraphe.R ###
## This script visualizes linked regions in a genome in conjunction with the RNAmapper.py script.
## This script requires 3 inputs:
## argument 1: the directory for the atMarkers.vcf files for the genome of interest
## argument 2: prefix for mutant marker files (should be everything before mut_chr#_atMarkers.vcf)
## argument 3: the path and name for the output plot (will have .jpg automatically appended to end)

# Variables
# Get current directory, read in arguments
currentDir <- getwd()
args <- commandArgs(T)
if (length(args) > 4) stop('Error in RNAMGenomeGraphe.R: Too many arguments.');
mutMarker_dir <- args[1] # directory to look for mutant marker files
mutMarker_prefix <- args[2] # everything in filename - including path - before _mut_chr
plotOut <- args[3] # Path to plot with plot name as file name (dont include .jpg)
plotName <- head(strsplit(tail(strsplit(plotOut, "/")[[1]], n=1), ".jpg")[[1]], n=1)[[1]]
writePlot <- if (length(args) < 4) FALSE else TRUE # should you include extra plot info: title, some tick labels, and axis labels

# Get list of chromosomes
setwd(mutMarker_dir)
chr_list <- system(paste(
  'ls',
  'grep -v "/"',
  'grep -E "mut.*_atMarkers\\.vcf"',
  'sed -E "s/.*mut_chr(.*)_atMarkers\\.vcf/\\1/"',
  'sort -n -k1',
  sep = ' | '),
  intern=T)
chr_list <- as.numeric(chr_list)
setwd(currentDir)

# Load in files
markerFiles <- list()
for(chr in chr_list){
  # get name of file
  name <- paste(mutMarker_prefix,"_mut_chr",chr,"_atMarkers.vcf",sep="")
  # save file
  chrfile <- read.delim(name, header=F)
  # save graph to markerFiles
  markerFiles[[chr]] <- chrfile[,c(2,22)]
}

# Put all the mutant zygosity data together into one big graph laid out by chromosome
jpeg(paste(plotOut, ".jpg", sep=""), width=2500, bg="white")   # setup graph environment
# put all graphs into one graphspace
par(mfrow=c(1,length(chr_list)),  # row of graphspace: 1 column, length of chr_list row
    mar=c(5,0.5,6,0.5), # margins of individual plots: c(bottom, left, top, right)
    oma=c(4,8.5,0,0), # margins of outer image: c(bottom, left, top, right)
    cex.lab=4, # magnification to be used for x and y labels relative to the current setting of cex. (default/current cex = 1)
    col.axis="white", # color to be used for axis annotation for each individual graph. Defaults to "black".
    lwd=2) # linewidth (default = 1)
# graph!
for(chr in chr_list){
  plot(markerFiles[[chr]][,1], 
       markerFiles[[chr]][,2], 
       ylim=c(0.5,1), 
       xlim=c(0,max(markerFiles[[chr]][,1])), 
       type="p", 
       xlab=chr)
abline(h = 1, col = "black", lty = 3, lwd=2)                           # dashed line at freq=1.0, i.e. homozygosity
}
# Write extra plot info if asked
if (writePlot){
  # Title
  mtext(paste(plotName, "Marker Frequency Across Genome", sep=" "),  # text
        side = 3, # side (1=bottom, 2=left, 3=top, 4=right)
        cex = 4, # character expansion factor (default/current cex = 1)
        font = 2, # font, 2 = bold
        line = - 4.5, # which margin line
        outer = TRUE) # use outer margins
  # Tick marks
  # Bottom limit (0.5)
  mtext("0.5",  # bottom label
        las = 1,
        adj = 0, # adjustment for string in reading direction (0=left/bottom, 1=right/top)
        side = 2, # side (1=bottom, 2=left, 3=top, 4=right)
        cex = 2, # character expansion factor (default/current cex = 1)
        line = 4, # which margin line
        at=c(0.135), # where to put mtext
        outer = TRUE) # use outer margins
  # Upper limit (1.0)
  mtext("1.0",  # bottom label
        las = 1,
        adj = 0, # adjustment for string in reading direction (0=left/bottom, 1=right/top)
        side = 2, # side (1=bottom, 2=left, 3=top, 4=right)
        cex = 2, # character expansion factor (default/current cex = 1)
        line = 4, # which margin line
        at=c(0.856), # where to put mtext
        outer = TRUE) # use outer margins
  # Axis titles
  # Y axis
  mtext("Marker Frequency",  # bottom label
        las = 0,
        side = 2, # side (1=bottom, 2=left, 3=top, 4=right)
        cex = 3, # character expansion factor (default/current cex = 1)
        line = 4.5, # which margin line
        outer = TRUE) # use outer margins
  # X axis
  mtext("Chromosome",  # bottom label
        las = 0,
        side = 1, # side (1=bottom, 2=left, 3=top, 4=right)
        cex = 3, # character expansion factor (default/current cex = 1)
        line = 2, # which margin line
        outer = TRUE) # use outer margins
}
dev.off()
