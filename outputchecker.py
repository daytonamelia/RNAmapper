#!/usr/bin/env python

#Imports
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-s", "--subject", help="", type=str, required=True)
    parser.add_argument("-q", "--query", help="", type=str, required=True)
    return parser.parse_args()

args = get_args()

subject = []
with open(args.subject, "r") as sf:
    for line in sf:
        line = line.strip()
        splitline = line.split("\t")
        subject.append(splitline[0])

with open(args.query, "r") as qf:
    for line in qf:
        line = line.strip()
        splitline = line.split("\t")
        if splitline[1] not in subject:
            print(splitline)