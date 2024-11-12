#!/usr/bin/env python

#Imports
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="Parses a WT and Mutant file, outputs a file of all SNPs in wt and mutant files, and outputs a file of all SNPs meeting linkage threshold.\nSee https://github.com/daytonamelia/RNAmapper for more details.")
    parser.add_argument("-wt", "--wtfile", help="Wildtype file to parse.", type=str, required=True)
    parser.add_argument("-mut", "--mutfile", help="Mutant file to parse.", type=str, required=True)
    parser.add_argument("-o", "--out", help="Output files prefix.", type=str, default="../RNAMapout/RNAmapper")
    parser.add_argument("-c", "--coverage", help="Minimum coverage for valid linkage threshold read.", type=int, default=10) # actual default = 25
    parser.add_argument("-z", "--zygosity", help="Minimum zygosity for valid linkage threshold read.", type=int, default=20)
    parser.add_argument("-n", "--neighbors", help="Amount of neighbors for sliding window average.", type=int, default=25) # actual default = 25 ; EITHER SIDE OF SNP!
    parser.add_argument("-lt", "--linkagethreshold", help="Minimum inclusive linkage threshold.", type=float, default=0.98)
    return parser.parse_args()

def vcf_lineparser(vcfline: str) -> list:
    """Parses a vcf line and returns the info of each column as a list."""
    vcfline_return = []
    # First split the read into columns for parsing, and clean out empty elements.
    splitline = re.split(r"\s", vcfline)
    cleanline = []
    for ele in splitline:
        if ele != "":
            cleanline.append(ele)
    vcfline_return.append(int(cleanline[0])) # CHROM
    vcfline_return.append(int(cleanline[1])) # POS
    vcfline_return.append(cleanline[2]) # ID
    vcfline_return.append(cleanline[3]) # REF
    vcfline_return.append(cleanline[4]) # ALT
    vcfline_return.append(int(cleanline[5])) # QUAL
    vcfline_return.append(cleanline[6]) # FILTER
    vcfline_return.append(cleanline[7]) # INFO
    vcfline_return.append(cleanline[7]) # FORMAT
    # Now parse the info column for stats
    splitinfo = re.split(r";", vcfline_return[7])
    infoinfo = [] # the info from the info
    indel = False
    for ele in splitinfo:
        # If read is an INDEL then change "INDEL" to True
        if "INDEL" in ele:
            indel = True
        # Grab the DP element
        if "DP" in ele:
            dpsplit = re.split(r"DP=", ele)
            infoinfo.append(int(dpsplit[1])) # DP
        # Grab the first four numbers in the I16 element
        if "I16" in ele:
            i16split = re.split(r"=|,", ele)
            infoinfo.append(int(i16split[1])) # FREF
            infoinfo.append(int(i16split[2])) # RREF
            infoinfo.append(int(i16split[3])) # FALT
            infoinfo.append(int(i16split[4])) # RALT
    # Calculate totals and ratios for reference and alternate.
    if int(infoinfo[0]) > 0: # if the depth is greater than 0
        infoinfo.append(infoinfo[1] + infoinfo[2]) # total reference
        infoinfo.append(infoinfo[3] + infoinfo[4]) # total alternate
        infoinfo.append(infoinfo[5]/(infoinfo[5] + infoinfo[6])) # reference ratio
        infoinfo.append(infoinfo[6]/(infoinfo[5] + infoinfo[6])) # alternate ratio
        infoinfo.append(round(infoinfo[7], 2) if infoinfo[7] > infoinfo[8] else round(infoinfo[8], 2)) # highest of the ratios!
    else:
        for i in range(0,6): # otherwise theres no reads so just append 0.0 5x times
            infoinfo.append(0.0)
    for ele in infoinfo:
        vcfline_return.append(ele)
    if indel:
        vcfline_return.append(True)
    else:
        vcfline_return.append(False)
    return vcfline_return

def vcf_altcleaner(vcfline: list) -> dict:
    """Given a list of a vcf line, removes ,<*> from the ALT column and keeps just the basepairs."""
    # Split the REF column by ,
    splitref = re.split(r",<\*>", vcfline[3])
    vcfline[3] = splitref[0]
    return vcfline

def vcffileparser(file: str) -> (dict):
    """Given a vcf file, returns a dictionary of all SNPs."""
    snps = {}
    with open(file, "r") as fh:
        for line in fh:
            line = line.strip()
            # Break line into vcfline
            vcfline = vcf_lineparser(line)
            # If ALT column is just <*> its not a SNP/indel and we can toss it
            if vcfline[4] == "<*>":
                continue
            vcfline = vcf_altcleaner(vcfline)
            snps[vcfline[1]] = vcfline
    return snps

def allelefreqcounter(snpdict: dict, zygo: int, coverage: int, wtcheck: bool) -> (set, list):
    """Given a dictionary of snps, filter for 0 reference allele and find the indels."""
    mapsnps = []
    indels = set()
    for pos, snp in snpdict.items():
        # Check for snps with high enough coverage only in wt (not in mutant!)
        if wtcheck:
            if int(snp[9]) < coverage:
                continue
        # Calls with 0 reference allele need to be removed because the vcf doesnt retain information on these calls.
        if int(snp[14]) == 0:
            continue
        # Collect bps of SNPs within 10bp of indels
        if snp[19]: # indel
            indelpos = snp[1]
            for i in range(indelpos-10, indelpos+10):
                # If it exists as a snp, is not an indel, and is not the indel in question, collect it for later deletion.
                if i in snpdict and i != indelpos and not snpdict[i][19]:
                    indels.add(i)
        # Grab high heterozygosity SNPs using the REFratio
        if wtcheck: # need to check zygosity in wt but not in mutant
            highzygo = 1 - (zygo/100)
            lowzygo = (zygo/100)
            if snp[16] > highzygo or snp[16] < lowzygo:
                continue
        if snp[0] not in indels: # quick filter step to save memory
            mapsnps.append(snp[1])
    return indels, mapsnps

def slidingwindowavg(mapsnps: list, reads: dict, neighborn: int) -> dict:
    snpslen = len(mapsnps)
    # If the list is lower than the number of neighbors, just adjust the neighbors to be the length of the list.
    if snpslen < neighborn:
        nieghborn = snpslen - 1
    # Finally calculate the sliding window:
    for i, snppos in enumerate(mapsnps):
        cumsum = 0
        # Before NEIGHBOR does just to the right:
        if i < neighborn:
            for j in range(neighborn):
                neighbor = mapsnps[i+j]
                cumsum += reads[neighbor][18]
            reads[snppos].append(round(cumsum/neighborn, 7))
        # After NEIGHBOR does just to the left:
        elif i > snpslen - neighborn - 1:
            for j in range(neighborn):
                neighbor = mapsnps[i-j]
                cumsum += reads[neighbor][18]
            reads[snppos].append(round(cumsum/neighborn, 7))
        # Within NEIGHBOR does both to left and right
        else:
            for j in range(-1*neighborn+1, neighborn):
                neighbor = mapsnps[i+j]
                cumsum += reads[neighbor][18]
            reads[snppos].append(round(cumsum/(2*neighborn-1), 7))
    return reads

# Read in user-passed arguments
args = get_args()

## STEP 1: ID SNPs in .vcf files
# Read in wt and mut files and get all snps into dictionaries
snps_wt = vcffileparser(args.wtfile)
snps_mut = vcffileparser(args.mutfile)

# Output all SNPs in wt and mut:
with open(f"{args.out}_wt_allALT.vcf", "w") as wtallalt:
    for snp in snps_wt.values():
        for ele in snp:
            wtallalt.write(f"{ele}\t")
        wtallalt.write(f"\n")
with open(f"{args.out}_mut_allALT.vcf", "w") as mutallalt:
    for snp in snps_mut.values():
        for ele in snp:
            mutallalt.write(f"{ele}\t")
        mutallalt.write(f"\n")

## STEP 2: Count allele frequency in wildtype
indels_wt, mapsnps_wt = allelefreqcounter(snps_wt, args.zygosity, args.coverage, True)
indels_mut, mapsnps_mutall = allelefreqcounter(snps_mut, args.zygosity, args.coverage, False)

# Filter indels
for pos in indels_wt:
    if pos in mapsnps_wt:
        mapsnps_wt.remove(pos)
    del snps_wt[pos]

for pos in indels_mut:
    if pos in mapsnps_mutall:
        mapsnps_mutall.remove(pos)
    del snps_mut[pos]

# STEP 3: Take sliding average of highest allele called at any position in the list of mapped snps and append it to the reads
# Make sure mutant snps are also present in wt
mapsnps_mut = []
for pos in mapsnps_mutall:
    if pos in mapsnps_wt:
        mapsnps_mut.append(pos)

snps_wt = slidingwindowavg(mapsnps_wt, snps_wt, args.neighbors)
snps_mut = slidingwindowavg(mapsnps_mut, snps_mut, args.neighbors)

# Write mutant markers to file
with open(f"{args.out}_mut_atMarkers.txt", "w") as mutmarkout:
    for snppos in mapsnps_mut:
        # if snps_mut[snppos][20] >= args.linkagethreshold:
        for ele in snps_mut[snppos]:
            mutmarkout.write(f"{ele}\t")
        mutmarkout.write(f"\n")
