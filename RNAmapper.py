#!/usr/bin/env python

#Imports
import argparse
import re

# Global Variables
COVERAGE = 10 # default = 25
ZYGOSITY = 20 # default = 20 or 25
WINDOW_NEIGHBORS = 50 # default = 50
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
    vcfline_return.append(cleanline[1])
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
    if (infoinfo[5] + infoinfo[6]) != infoinfo[0]: # if the sum of the totals doesnt equal the depth theres a problem!
        print(f"The read totals of line {vcfline} do not match the depth!")
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
        if wtvcfline[2] == "<*>":
            continue
        # It's a SNP! We need the relevant information (DP, I16) from the INFO column.
        wtvcfline = vcf_infoparser(wtvcfline)
        # Clean the <*> out of the ALT column of the SNP file.
        wtvcfline = vcf_altcleaner(wtvcfline)
        # SNP is ready to be added to dictionary.
        snps_wt[wtvcfline[0]] = wtvcfline

# Next read in mutant file and get all snps into snps_mut
snps_mut = {}
#TODO

# snps_allele dictionary value has the following layout:
# pos, ref, alt, dp, fREF, rREF, fALT, rALT, REFtotal, ALTtotal, REFratio, ALTratio, INDEL status

## STEP 2: Count allele frequency in wildtype
# Filter the wt snps list 
for snp in snps_wt.values():
    # If the read does not have enough coverage (the DP) then it should not be counted.
    # However, we do want to output this above! Just in case the causative mutation has very low depth.
    print('--')
    print(snp)
    if snp[3] < COVERAGE:
        continue
    # Remove SNPs within 10bp of indels
    if snp[12]: # indel
        indelpos = snp[0]

    # TODO: good SNPs will have a high % of heterozygosity in the calls
    # TODO: apparently calls with 0 ref need to removed

    print(snp)
    # TODO: sliding average