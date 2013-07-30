#!/usr/bin/python
#
# Prephix to PhenoLink Generator
#
# Usage: pre2phe.py -ref <ref base file> -snp <snp file> -out <PhenoLink output file>
#
# Takes in a one SNP loci output file (possibly containing multiple strains) and one reference base file generated from prephix,
# and generates a PhenoLink export file.
#

VERSION="1.0"

import argparse
import re
import os
import sys

print "\nPrephix to PhenoLink Generator v{0}\n\n".format(VERSION)

# Setup argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('--ref',required=True,type=str, nargs=1,help="Reference base output file from prephix.")
parser.add_argument('--snp',required=True,type=str,nargs=1,help="SNP output file from prephix.")
parser.add_argument('--out',required=True,type=str,nargs=1,help="PhenoLink output file name.")
args = parser.parse_args()

reffile = open(args.ref[0],'r')
snpfile = open(args.snp[0],'r')
outfile = open(args.out[0],'w')

# Read the reference base file into memory (refTable dictionary keyed off of locus with base for value).
refTable = dict()
lineNumber=1
print "Reading base reference file..." ,
try:
    refRegex = re.compile("^(?P<locus>\d+)\t(?P<base>\D)$")
    for line in reffile:
        refMatch = refRegex.match(line)
        if refMatch == None:
            print "*** WARNING: Unrecognized format at line {0}: {1}".format(lineNumber,line)
        else:
            refTable[refMatch.group('locus')] = refMatch.group('base')
        lineNumber += 1

finally:
    reffile.close()

print "done. Read {0} lines.".format(lineNumber)

# Now collect snp data into the snpTable, which is a dict keyed off of StrainID with a value being a list of SNP loci names
# that it shares.  SNP loci names are <refbase>_<locus>_<snpbase>.
# At the same time, collect all known SNP loci names into a unique list.
snpTable = dict()
snpLociList = []
lineNumber=1
print "Processing SNP file..." ,
try:
    snpRegex = re.compile("^(?P<strainid>.+?)\t(?P<locus>\d+)\t(?P<base>.)$")
    for line in snpfile:
        snpMatch = snpRegex.match(line)
        if snpMatch == None:
            print "*** WARNING: Unrecognized SNP format at line {0}: {1}".format(lineNumber,line)
        else:
            strainid = snpMatch.group('strainid')
            strainLocus = snpMatch.group('locus')
            strainBase = snpMatch.group('base')
            snpLociName = "{0}_{1}_{2}".format(refTable[snpMatch.group('locus')],snpMatch.group('locus'),snpMatch.group('base'))

            # Associating snp loci to strain id
            if strainid in snpTable.keys():
                snpTable[strainid].append(snpLociName)
            else:
                snpTable[strainid] = dict()
                snpTable[strainid] = [snpLociName]
                
            # Saving snp loci to unique list.
            if not snpLociName in snpLociList:
                snpLociList.append(snpLociName)

        lineNumber += 1
finally:
    snpfile.close()

print "done. Processed {0} lines.".format(lineNumber)

print "Writing PhenoLink export data..." ,

try:
    # First write out the column headings.  The first is the Strain ID, then all the unique SNP loci names.
    outfile.write("StrainID")
    for lociName in snpLociList:
        outfile.write("\t{0}".format(lociName))

    outfile.write("\n")

    # For each strain, fill out its row values for each snp loci name.
    for strainid in snpTable.keys():
        outfile.write(strainid)
        for lociName in snpLociList:
            if lociName in snpTable[strainid]:
                outfile.write("\t1")
            else:
                outfile.write("\t0")
        outfile.write("\n")

finally:
    outfile.close()

print "done.\n"

print "Final output file is {0}\n".format(args.out[0]);
