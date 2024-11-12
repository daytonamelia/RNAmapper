#!/usr/bin/env python

#Imports
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-wt", "--wtfile", help="", type=str, required=True)
    parser.add_argument("-mut", "--mutfile", help="", type=str, required=True)
    parser.add_argument("-o", "--out", help="", type=str, default="../RNAMapout/RNAmapper")
    parser.add_argument("-c", "--coverage", help="", type=int, default=10) # actual default = 25
    parser.add_argument("-z", "--zygosity", help="", type=int, default=20)
    parser.add_argument("-n", "--neighbors", help="", type=int, default=3) # actual default = 25 ; EITHER SIDE OF SNP!
    parser.add_argument("-lt", "--linkagethreshold", help="", type=float, default=0.98)
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
        # if (infoinfo[5] + infoinfo[6]) != infoinfo[0]: # if the sum of the totals doesnt equal the depth theres a problem!
        #     print(f"The read totals of line {vcfline} do not match the depth!")
        infoinfo.append(infoinfo[5]/(infoinfo[5] + infoinfo[6])) # reference ratio
        infoinfo.append(infoinfo[6]/(infoinfo[5] + infoinfo[6])) # alternate ratio
        infoinfo.append(round(infoinfo[7], 2) if infoinfo[7] > infoinfo[8] else round(infoinfo[8], 2)) # highest of the ratios!
    else:
        for i in range(0,6): # otherwise theres no reads so just append 0 5 times
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

def vcffileparser(file: str) -> (dict, dict):
    """Given a vcf file, returns a dictionary of all reads with DP > coverage and a dictionary of all SNPs with DP > coverage."""
    reads = {}
    snps = {}
    with open(file, "r") as fh:
        for line in fh:
            line = line.strip()
            # Break line into vcfline
            vcfline = vcf_lineparser(line)
            # If the read was below coverage, it returned 0 and we can toss it
            if vcfline == 0:
                continue
            reads[vcfline[1]] = vcfline
            # If ALT column is just <*> its not a SNP/indel
            if vcfline[4] == "<*>":
                continue
            vcfline = vcf_altcleaner(vcfline)
            snps[vcfline[1]] = vcfline
    return (reads, snps)

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
        print('--')
        print(i, snppos)
        cumsum = 0
        # Before NEIGHBOR does just to the right:
        if i < neighborn:
            print("BEFORE")
            for j in range(neighborn):
                neighbor = mapsnps[i+j]
                print(reads[neighbor])
                cumsum += reads[neighbor][18]
            print(cumsum)
            print(round(cumsum/neighborn, 7))
            reads[snppos].append(round(cumsum/neighborn, 7))
        # After NEIGHBOR does just to the left:
        elif i > snpslen - neighborn - 1:
            print("AFTER")
            for j in range(neighborn):
                neighbor = mapsnps[i-j]
                print(reads[neighbor])
                cumsum += reads[neighbor][18]
            print(cumsum)
            print(round(cumsum/neighborn, 7))
            reads[snppos].append(round(cumsum/neighborn, 7))
        # Within NEIGHBOR does both to left and right
        else:
            print("WITHIN")
            for j in range(-1*neighborn+1, neighborn):
                neighbor = mapsnps[i+j]
                print(reads[neighbor])
                cumsum += reads[neighbor][18]
            print(cumsum)
            print((2*neighborn-1))
            print(round(cumsum/(2*neighborn-1), 7))
            reads[snppos].append(round(cumsum/(2*neighborn-1), 7))
    return reads

# Read in user-passed arguments
args = get_args()

## STEP 1: ID SNPs in .vcf files
# Read in wt and mut files and get all reads and snps into dictionaries
reads_wt, snps_wt = vcffileparser(args.wtfile)
reads_mut, snps_mut = vcffileparser(args.mutfile)

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
print(mapsnps_wt)

mapsnps_mut = []
for pos in mapsnps_mutall:
    if pos in mapsnps_wt:
        mapsnps_mut.append(pos)

print(mapsnps_mut)


reads_wt = slidingwindowavg(mapsnps_wt, reads_wt, args.neighbors)
print('MUTANT')
reads_mut = slidingwindowavg(mapsnps_mut, reads_mut, args.neighbors)

with open(f"{args.out}_mut_atMarkers.txt", "w") as mutmarkout:
    for snppos in mapsnps_mut:
        # if reads_mut[snppos][18] >= args.linkagethreshold:
        for ele in reads_mut[snppos]:
            mutmarkout.write(f"{ele}\t")
        mutmarkout.write(f"\n")
