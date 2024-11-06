#!/usr/bin/env python

#Imports
import argparse
import re

# Global Variables
COVERAGE = 10 # default = 25
ZYGOSITY = 20 # default = 20 or 25
NEIGHBORS = 3 # default = 50 ; EITHER SIDE OF SNP!
LINKAGE_THRESHOLD = 0.98 # default = 0.98

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-wt", "--wtfile", help="", type=str, required=True)
    parser.add_argument("-mut", "--mutfile", help="", type=str, required=True)
    parser.add_argument("-o", "--out", help="", type=str, default="RNAmapper")
    return parser.parse_args()

def vcf_lineparser(vcfline: str) -> list:
    """Parses a vcf line and returns the relevant info of each column as a list.
    Relevant info: REF, ALT, INFO"""
    # VCF general TSV format: CHROM POS ID  REF ALT QUAL    FILTER  INFO(contains DP and I16)    FORMAT
    vcfline_return = []
    # First split the read into columns for parsing, and clean out empty elements.
    splitline = re.split(r"\s", vcfline)
    cleanline = []
    for ele in splitline:
        if ele != "":
            cleanline.append(ele)
    # Grab POS, REF, ALT, INFO
    vcfline_return.append(int(cleanline[1]))
    vcfline_return.append(cleanline[3])
    vcfline_return.append(cleanline[4])
    vcfline_return.append(cleanline[7])
    return vcfline_return

def vcf_infoparser(vcfline: dict) -> dict:
    """Given a list of relevant vcf line information, takes relevant DP and I16 information from the INFO column and appends it.
    Relevant info includes the reference info, alternate info, total reference, total alternate, reference ratio, and alternate ratio.
    Also marks if the read is an indel."""
    # Split the INFO column by semicolons and add relevant info to infoinfo
    splitinfo = re.split(r";", vcfline[3])
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
    infoinfo.append(infoinfo[1] + infoinfo[2]) # total reference
    infoinfo.append(infoinfo[3] + infoinfo[4]) # total alternate
    # if (infoinfo[5] + infoinfo[6]) != infoinfo[0]: # if the sum of the totals doesnt equal the depth theres a problem!
    #     print(f"The read totals of line {vcfline} do not match the depth!")
    infoinfo.append(infoinfo[5]/(infoinfo[5] + infoinfo[6])) # reference ratio
    infoinfo.append(infoinfo[6]/(infoinfo[5] + infoinfo[6])) # alternate ratio
    # Delete the INFO column as we no longer need it, and add the elements from infoinfo
    vcfline.pop()
    for ele in infoinfo:
        vcfline.append(ele)
    if indel:
        vcfline.append(True)
    else:
        vcfline.append(False)
    return vcfline

def vcf_altcleaner(vcfline: dict) -> dict:
    """Given a list of a vcf line, removes ,<*> from the ALT column and keeps just the basepairs."""
    # Split the REF column by ,
    splitref = re.split(r",<\*>", vcfline[2])
    vcfline[2] = splitref[0]
    return vcfline

# Read in user-passed arguments
args = get_args()

## STEP 1: ID SNPs in .vcf files

# Read in wildtype file and get all snps into snps_wt
snps_wt = {}
with open(args.wtfile, "r") as wtfile:
    for wtline in wtfile:
        # Break each line into a list for easy parsing.
        wtline = wtline.strip()
        # TODO: CHECK FOR HEADER!
        wtvcfline = vcf_lineparser(wtline)
        # If the ALT column is just <*>, it is not of note and we can just continue.
        if wtvcfline[2] == "<*>" or "," not in wtvcfline[2]:
            continue
        # It's a SNP! We need the relevant information (DP, I16) from the INFO column.
        wtvcfline = vcf_infoparser(wtvcfline)
        # Clean the <*> out of the ALT column of the SNP file.
        wtvcfline = vcf_altcleaner(wtvcfline)
        # SNP is ready to be added to dictionary.
        snps_wt[wtvcfline[0]] = wtvcfline

# Next read in mutant file and get all snps into snps_mut
snps_mut = {}
with open(args.mutfile, "r") as mutfile:
    for mutline in mutfile:
        # Break each line into a list for easy parsing.
        mutline = mutline.strip()
        # TODO: CHECK FOR HEADER!
        mutvcfline = vcf_lineparser(mutline)
        # If the ALT column is just <*>, it is not of note and we can just continue.
        if mutvcfline[2] == "<*>":
            continue
        # It's a SNP! We need the relevant information (DP, I16) from the INFO column.
        mutvcfline = vcf_infoparser(mutvcfline)
        # Clean the <*> out of the ALT column of the SNP file.
        mutvcfline = vcf_altcleaner(mutvcfline)
        # SNP is ready to be added to dictionary.
        snps_mut[mutvcfline[0]] = mutvcfline

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

# snps_allele dictionary value has the following layout:
# pos, ref, alt, dp, fREF, rREF, fALT, rALT, REFtotal, ALTtotal, REFratio, ALTratio, INDEL status

## STEP 2: Count allele frequency in wildtype
# Filter the wt snps list
mapsnps_wt = []
for pos, snp in snps_wt.items():
    # If the read does not have enough coverage (the DP) then it should not be counted.
    # However, we do want to output this above! Just in case the causative mutation has very low depth.
    if snp[3] < COVERAGE:
        continue
    # Calls with 0 reference allele need to be removed because the vcf doesnt retain information on these calls.
    if snp[8] == 0:
        continue
    # Collect bps of SNPs within 10bp of indels
    indels = set()
    if snp[12]: # indel
        indelpos = snp[0]
        # for i in range(indelpos-10, indelpos+10):
        #     # If it exists as a snp, is not an indel, and is not the indel in question, collect it for later deletion.
        #     if i in snps_wt and i != indelpos and not snps_wt[i][12]:
        #         indels.add(i)
    # Grab high heterozygosity SNPs using the REFratio
    highzygo = 1 - (ZYGOSITY/100)
    lowzygo = (ZYGOSITY/100)
    if snp[10] <= highzygo or snp[11] >= lowzygo:
        if snp[0] not in indels: # quick filter step to save memory
            mapsnps_wt.append(snp[0])
# Filter out the SNPs that were within 10bp of indels.
for pos in indels:
    mapsnps_wt.remove(pos)
    del snps_wt[pos]

# Filter the mut list
mapsnps_mut = []
for pos, snp in snps_mut.items():
    if snp[3] < COVERAGE:
        continue
    if snp[8] == 0:
        continue
    indels = set()
    if snp[12]: # indel
        indelpos = snp[0]
        # for i in range(indelpos-10, indelpos+10):
        #     if i in snps_mut and i != indelpos and not snps_mut[i][12]:
        #         indels.add(i)
    highzygo = 1 - (ZYGOSITY/100)
    lowzygo = (ZYGOSITY/100)
    if snp[10] <= highzygo or snp[10] >= lowzygo:
        if snp[0] not in indels: # quick filter step to save memory
            mapsnps_mut.append(snp[0])
for pos in indels:
    mapsnps_mut.remove(pos)
    del snps_mut[pos]

# Take sliding average of highest allele called at any position in wt and append it to the SNP call
wtsnpslen = len(mapsnps_wt)
# If the list is lower than the number of neighbors, tell user and just adjust the neighbors to be the length of the list.
if wtsnpslen < NEIGHBORS:
    # print(f"Less than {NEIGHBORS} SNPs in wildtype so sliding window neighbors has been adjusted to {wtsnpslen - 1}!")
    NEIGHBORS = wtsnpslen - 1
# Finally calculate the sliding window:
for i, snppos in enumerate(mapsnps_wt):
    cumsum = 0
    # Before NEIGHBOR does just to the right:
    if i < NEIGHBORS:
        for j in range(i + NEIGHBORS - (i-1)):
            neighbor = mapsnps_wt[j]
            if snps_wt[neighbor][10] > snps_wt[neighbor][11]:
                cumsum += snps_wt[neighbor][10]
            else:
                cumsum += snps_wt[neighbor][11]
        snps_wt[snppos].append((cumsum/(NEIGHBORS+1)))
    # After NEIGHBOR does just to the left:
    elif i >= wtsnpslen - NEIGHBORS - 1:
        window_size = i + NEIGHBORS - (i-1)
        for j in range(i + NEIGHBORS - (i-1)):
            neighbor = mapsnps_wt[i-j]
            if snps_wt[neighbor][10] > snps_wt[neighbor][11]:
                cumsum += snps_wt[neighbor][10]
            else:
                cumsum += snps_wt[neighbor][11]
        snps_wt[snppos].append(cumsum/(NEIGHBORS+1))
    # Within NEIGHBOR does both to left and right
    else:
        for j in range(i-NEIGHBORS, i+NEIGHBORS+1):
            neighbor = mapsnps_wt[j]
            if snps_wt[neighbor][10] > snps_wt[neighbor][11]:
                cumsum += snps_wt[neighbor][10]
            else:
                cumsum += snps_wt[neighbor][11]
        snps_wt[snppos].append(cumsum/(NEIGHBORS*2+1))

# Take sliding average of highest allele called at any position in mut and append it to the SNP call
mutsnpslen = len(mapsnps_mut)
if mutsnpslen < NEIGHBORS:
    # print(f"Less than {NEIGHBORS} SNPs in mutant so sliding window neighbors has been adjusted to {mutsnpslen - 1}!")
    NEIGHBORS = mutsnpslen - 1
for i, snppos in enumerate(mapsnps_mut):
    cumsum = 0
    if snppos not in snps_wt: # if it isnt in the wt, delete it
        continue
    if i < NEIGHBORS:
        for j in range(i + NEIGHBORS - (i-1)):
            neighbor = mapsnps_mut[j]
            if snps_mut[neighbor][10] > snps_mut[neighbor][11]:
                cumsum += snps_mut[neighbor][10]
            else:
                cumsum += snps_mut[neighbor][11]
        snps_mut[snppos].append((cumsum/(NEIGHBORS+1)))
    elif i >= mutsnpslen - NEIGHBORS - 1:
        window_size = i + NEIGHBORS - (i-1)
        for j in range(i + NEIGHBORS - (i-1)):
            neighbor = mapsnps_mut[i-j]
            if snps_mut[neighbor][10] > snps_mut[neighbor][11]:
                cumsum += snps_mut[neighbor][10]
            else:
                cumsum += snps_mut[neighbor][11]
        snps_mut[snppos].append(cumsum/(NEIGHBORS+1))
    else:
        for j in range(i-NEIGHBORS, i+NEIGHBORS+1):
            neighbor = mapsnps_mut[j]
            if snps_mut[neighbor][10] > snps_mut[neighbor][11]:
                cumsum += snps_mut[neighbor][10]
            else:
                cumsum += snps_mut[neighbor][11]
        snps_mut[snppos].append(cumsum/(NEIGHBORS*2+1))

with open(f"{args.out}_mut_atMarkers.txt", "w") as mutmarkout:
    mutmarkout.write(f"#CHROM\tPOS\tREF\tALT\tDP\tFREF\tRREF\tFALT\tRALT\tTOTALREF\tTOTALALT\tREFRATIO\tALTRATIO\tINDEL\tSLIDINGAVG\n")
    for snp in snps_mut.values():
        if float(snp[-1]) <= float(LINKAGE_THRESHOLD): # if the sliding average is above linkage threshold, write
            for ele in snp:
                mutmarkout.write(f"{ele}\t")
            mutmarkout.write(f"\n")

