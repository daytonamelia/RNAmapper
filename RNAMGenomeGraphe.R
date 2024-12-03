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
mutMarker_dir <- "~/bioinfo/Bi624/FW_Genetics/RNAMapout/"
mutMarker_prefix <- "RNAmapper"
plotOut <- "Genometest"
plotOut <- paste(plotOut, ".jpg", sep="")

# Get list of chromosomes
setwd(mutMarker_dir)
chr_list <- system(paste(
  'ls',
  'grep -v "/"',
  'grep -E "mut.*\\.vcf"',
  'grep -E -v "mut_.*_"',
  'sed -E "s/.*mut_(.*)\\.vcf/\\1/"',
  'sort -n -k1.4',
  sep = ' | '),
  intern=T)
setwd(currentDir)

# Load in files
markerFiles <- list()
for(chr in chr_list){
  # get name of file
  name <- paste(mutMarker_prefix,"_mut_",chr,"_atMarkers.txt",sep="")
  # make graph of file
  graph <- read.delim(name, header=T)
  # save graph to markerFiles
  all[[i]] <- graph[,c(2,21)]
}
##### put all the mutant zygosity data together into one big graph laid out by chromosome
jpeg(paste(), width=2500, bg="white")   # setup graph environment
par(mfrow=c(1,len(chr_list)), mar=c(5,0.5,0,0), cex.lab=4, col.axis="white")   # put all graphs into one graphspace
# graph!
for(chr in chr_list){
  plot(all[[chr]][,1], all[[chr]][,2], ylim=c(0.5,1), xlim=c(0,max(all[[chr]][,1])), type="p", xlab=chr)
}

