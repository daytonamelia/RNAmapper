#!/usr/bin/env python

#Imports
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="Parses a WT and Mutant file, outputs a file of all SNPs in wt and mutant files, and outputs a file of all SNPs meeting linkage threshold.\nSee https://github.com/daytonamelia/RNAmapper for more details.")
    parser.add_argument("-wt", "--wtfile", help="Wildtype file to parse.", type=str, required=True)
    parser.add_argument("-mut", "--mutfile", help="Mutant file to parse.", type=str, required=True)
    parser.add_argument("-o", "--out", help="Output files prefix.", type=str, default="../RNAMapout/RNAmapper1")
    parser.add_argument("-c", "--coverage", help="Minimum coverage for valid linkage threshold read.", type=int, default=10) # actual default = 25
    parser.add_argument("-z", "--zygosity", help="Minimum zygosity for valid linkage threshold read.", type=int, default=20)
    parser.add_argument("-n", "--neighbors", help="Amount of neighbors for sliding window average.", type=int, default=25) # actual default = 25 ; EITHER SIDE OF SNP!
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

def vcffileparser(infile: str, outfile: str, zygo: int, coverage: int, wtcheck:bool) -> (dict, set, int, int):
    """Given a vcf file, returns a dictionary of all SNPs and writes mapping SNPs to output file. Also returns a set of indels for later checking, how many total snps, and how many mapping snps."""
    snps = {}
    indels = set()
    snps_num = 0
    mapsnps_num = 0
    with open(infile, "r") as infile, open(outfile, "w") as outfile:
        for line in infile:
            line = line.strip()
            # Break line into vcfline
            vcfline = vcf_lineparser(line)
            # If ALT column is just <*> its not a SNP/indel and we can toss it
            if vcfline[4] == "<*>":
                continue
            vcfline = vcf_altcleaner(vcfline)
            # Write to outfile and count for stats
            for ele in vcfline:
                outfile.write(f"{ele}\t")
            outfile.write(f"\n")
            snps_num +=1
            # Grab indels for later checking
            if vcfline[19]: # indel
                indels.add(vcfline[1])
            # Only check coverage and zygosity in wt
            if wtcheck:
                # Check coverage
                if int(vcfline[9]) < coverage:
                    continue
                highzygo = 1 - (zygo/100) # high zygosity range
                lowzygo = (zygo/100) # low zygosity range
                # Check zygosity of REFratio
                if vcfline[16] > highzygo or vcfline[16] < lowzygo:
                    continue
            # Write to dictionary and stats
            snps[vcfline[1]] = vcfline
            mapsnps_num += 1
    return snps, indels, snps_num, mapsnps_num

def findfilterindels(mappingsnps: dict, indelset: set) -> dict:
    """Given a dictionary of SNPs and a set of indel positions, find positions with 10 bp of indels and delete them."""
    # Collect SNPs
    postodelete = set()
    for indelpos in indelset:
        print('--')
        print(indelpos)
        indel = mappingsnps[indelpos]
        for i in range(indelpos-10, indelpos+10):
            print(i)
            print(mappingsnps[i])
            # If it exists as a SNP, is not an indel, and is not the indel in question, grab it for deletion.
            if i in mappingsnps and i != indelpos and not mappingsnps[i][19]:
                del mappingsnps[i]
    return mappingsnps

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

# # ID SNPs in .vcf files
# Read in wt and mut files, write all SNPs to output file and get all mapping SNPs into a dictionary
mapsnps_wt, indels_wt, snpsnum_wt, mapsnpsnum_wt = vcffileparser(args.wtfile, f"{args.out}_wt_allALT.vcf", args.zygosity, args.coverage, True)
mapsnps_mut, indels_mut, snpsnum_mut, mapsnpsnum_mut = vcffileparser(args.mutfile, f"{args.out}_mut_allALT.vcf", args.zygosity, args.coverage, False)

print(mapsnpsnum_wt)
print(mapsnpsnum_mut)
sleep(10)

# Filter indels
mapsnps_wt = findfilterindels(mapsnps_wt, indels_wt)
mapsnps_mut = findfilterindels(mapsnps_mut, indels_mut)

print(mapsnps_wt.keys())

# # STEP 3: Take sliding average of highest allele called at any position in the list of mapped snps and append it to the reads
# # Make sure mutant snps are also present in wt
# mapsnps_mut = []
# for pos in mapsnps_mut:
#     if pos not in mapsnps_wt:
#         del mapsnps_mut[pos]

# snps_wt = slidingwindowavg(mapsnps_wt, snps_wt, args.neighbors)
# snps_mut = slidingwindowavg(mapsnps_mut, snps_mut, args.neighbors)

# # Write mutant markers to file
# with open(f"{args.out}_mut_atMarkers.vcf", "w") as mutmarkout:
#     for snppos in mapsnps_mut:
#         # if snps_mut[snppos][20] >= args.linkagethreshold:
#         for ele in snps_mut[snppos]:
#             mutmarkout.write(f"{ele}\t")
#         mutmarkout.write(f"\n")

# # Write stats file
# # snps_wt, snps_mut, mapsnps_wt, mapsnps_mutall, mapsnps_mut
# with open(f"{args.out}_stats.txt", "w") as statsfile:
#     statsfile.write(f"Total WT SNPs: {snpsnum_wt}\n")
#     statsfile.write(f"Total MUT SNPs: {snpsnum_mut}\n")
#     statsfile.write(f"Total mapping SNPs in WT: {mapsnpsnum_wt}\n")
#     statsfile.write(f"Total mapping SNPs in MUT pre-filtering: {mapsnpsnum_mut}\n")
#     statsfile.write(f"Total mapping SNPs in MUT post-filtering: {len(mapsnps_mut)}\n")
#     statsfile.write(f"**Filtering for the mapping mutant SNPs finds positions also present in the mapping WT SNPs")
    