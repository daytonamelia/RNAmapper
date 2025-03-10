# RNAmapper

Forward genetics, in which mutations are randomly generated in the genome of an animal model, is a classic method to uncover genes that affect biological processes and provide insight into human genetic disorders. Yet, despite being used as a major tool in genetics for decades, there is still a major challenge in identifying the causative mutations and underlying genes. [Dr. Adam Miller's lab](https://www.adammillerlab.com/) at the University of Oregon has developed an RNA-seq-based method that can identify the causative mutation and examine the related gene regulatory changes in mutants. Implementation of the methodology has been hampered by the lack of a modern, accessible bioinformatic pipeline that would provide broad useability to the community.

This project aims to facilitate the ability of geneticists to gain insight into the genetics of development and disease. This repository is the backbone of scripts that are used for the tool's [web-based interface](https://github.com/ramzymulla/RNAmApp).

## Relevant Papers:

Miller, A. C., Obholzer, N. D., Shah, A. N., Megason, S. G., & Moens, C. B. (2013). RNA-seq-based mapping and candidate identification of mutations from forward genetic screens. Genome research, 23(4), 679–686. [https://doi.org/10.1101/gr.147322.112](https://genome.cshlp.org/content/23/4/679.long)

## Usage:

#### <ins>RNAmapper.py</ins>

RNAmapper.py requires all of your chromosomes split into their own .vcf file (see previous section) and you can run each chromosome file individually or in parallel.

`RNAmapper.py -wt WILDTYPE_FILE -mut MUTANT_FILE -o OUTPUT_DIRECTORY -ch CHROMOSOME_NUMBER`

-wt : .vcf file with SNPs for wildtype.

-mut : .vcf file with SNPs for mutant.

-o : Output directory and output prefix.

-ch : Chromosome number.

-c : Optional argument for coverage cutoff. Default = 20.

-z : Optional argument for zygosity threshold. Default = 25.

-n : Amount of neighbors on either side of SNP for sliding window average calculation. Default = 25.

Once the script is complete, output files will be present in either the specified output directory. These files include a marker file and a stats file for each mutant chromosome.

#### <ins>RNAGenomeGraphe.R</ins>

You can use this Rscript to visualize linked regions in a genome.

`Rscript RNAMGenomeGrapher.R MARKER_DIRECTORY OUTPUT_PREFIX PLOT_DIRECTORY`

-1 : The directory for the marker files for the genome of interest.

-2 : The prefix for the mutant marker files (specified above when creating the files) and consists of everything before “mut_chr#_atMarkers.vcf”

-3 : The path and name for the output plot. This will have .jpg automatically appended to end.

-4 : An optional argument to include extra plot info such as a title, axis labels, and improved tick marks. Can be specified with “TRUE.”

#### <ins>RNAGenomeGraphe.R</ins>

You can use this Rscript to visualize linked regions in a specific chromosome.

`Rscript RNAMChromosomeGraphe.R`

-1 : The path to the marker file of interest.

-2 : The path to the stats file of which to append the region of linkage.

-3 : The path to the output plot which will have chromosome number and .jpg automatically appended to end.

-4 : Optional argument for linkage ratio in a proportion. Default = 0.98.
