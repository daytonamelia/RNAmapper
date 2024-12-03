### RNAMGenomeGraphe.R ###
## This script visualizes linked regions in a genome in conjunction with the RNAmapper.py script.
## This script requires 2 inputs:
## argument 1: the directory for the atMarkers.vcf files for the genome of interest
## argument 2: prefix for mutant marker files (should be everything before mut_chr#_atMarkers.vcf)
## argument 3: the name for the output plot (will have .jpg automatically appended to end)

# Variables
# Get current directory, read in arguments
currentDir <- getwd()
args <- commandArgs(T)
if (length(args) > 3) stop('Error in RNAMGenomeGraphe.R: Too many arguments.');
mutMarker_dir <- args[1]
mutMarker_prefix <- args[2]
plotOut <- args[3]
plotOut <- paste(plotOut, ".jpg", sep="")

# Testing variables DELETE LATER
currentDir <- getwd()
mutMarker_dir <- "~/bioinfo/Bi624/FW_Genetics/GenomeTest/"
mutMarker_prefix <- "Hox40_py"
plotOut <- "Genometest"
plotOut <- paste(plotOut, ".jpg", sep="")

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
  # make graph of file
  graph <- read.delim(name, header=F)
  # save graph to markerFiles
  markerFiles[[chr]] <- graph[,c(2,22)]
}

# Put all the mutant zygosity data together into one big graph laid out by chromosome
jpeg(plotOut, width=2500, bg="white")   # setup graph environment
par(mfrow=c(1,length(chr_list)), mar=c(5,0.5,0,0), cex.lab=4, col.axis="white")   # put all graphs into one graphspace
# graph!
for(chr in chr_list){
  plot(markerFiles[[chr]][,1], markerFiles[[chr]][,2], ylim=c(0.5,1), xlim=c(0,max(markerFiles[[chr]][,1])), type="p", xlab=chr)
}
dev.off()
