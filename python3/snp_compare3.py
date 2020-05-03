#!/usr/bin/python3
#
# SNP Comparator
#
# Takes in an SNP file (formatted like the snp file output from prephix) and a
# list of strain IDs to compare and/or exclude.  SNP Comparator finds all the
# common SNP Loci between the compared Strains and then removes the SNP Loci
# that are common across the exclude Strains.
#
# 7/31/2013 -- V1.1 Added output file name parameter and added more reporting stats.
# 7/31/2013 -- V1.2 Fixed sorting.

import os
import sys
import re

VERSION = "1.1"

print("SNP Comparator v{0}\n".format(VERSION))

if len(sys.argv) < 2:
    print("Usage: {0} --snpfile <snpfile> --compare <compare strain 1> <compare strain 2> ... [--exclude <exclude strain 1> <exclude strain 2> ...] [--outfile <filename>]\n".format(sys.argv[0]))
    sys.exit(1)

mode=""
compareList = []
excludeList = []
snpFilename=""
outFilename=""
for arg in sys.argv:

    if arg == sys.argv[0]:
        continue

    if arg == "--compare":
        mode="compare"
        continue
    elif arg == "--exclude":
        mode="exclude"
        continue
    elif arg == "--snpfile":
        mode="snpfile"
        continue
    elif arg == "--outfile":
        mode="outfile"
        continue

    if mode == "compare":
        compareList.append(arg)
        continue
    elif mode == "exclude":
        excludeList.append(arg)
        continue
    elif mode == "snpfile":
        snpFilename = arg
        continue
    elif mode == "outfile":
        outFilename = arg
        continue

if outFilename == "":
    outFilename = "{0}.compared.txt".format(snpFilename)

# Track which input strains have been included or excluded -- used for warning
# if some strains were not found in snp file
unusedCompareStrains = compareList[:]
unusedExcludeStrains = excludeList[:]

# Read in SNP file and populate the include/exclude dictionaries.  Key is StrainID, value
# is a set of <loci><base>.  Only load the SNP strain data for strains being compared
# or excluded.
print("Reading snp file...")
snpCompareTable = dict()
snpExcludeTable = dict()
try:
    currentStrainID=""
    currentStrainLociList = []
    snpFile = open(snpFilename,"r")
    snpLineRe = re.compile("^(?P<strainID>.+?)\t(?P<locus>\d+)\t(?P<base>.)$")
    lineNumber = 0

    for line in snpFile:
        snpLineMatch = snpLineRe.search(line)

        if snpLineMatch == None:
            print("*** ERROR: Bad formatted line at line {0}: {1}".format(lineNumber,line))
            sys.exit(1)

        else:
            matchStrainID = snpLineMatch.group('strainID')
            matchLocus = snpLineMatch.group('locus')
            matchBase = snpLineMatch.group('base')

            if  not (matchStrainID in compareList) and not (matchStrainID in excludeList):
                continue

            if currentStrainID != matchStrainID:

                if (currentStrainID != "") and (currentStrainLociList != ""):

                    if currentStrainID in compareList:
                        snpCompareTable[currentStrainID] = set(currentStrainLociList)

                        if currentStrainID in unusedCompareStrains:
                            unusedCompareStrains.remove(currentStrainID)

                    else:
                        snpExcludeTable[currentStrainID] = set(currentStrainLociList)

                        if currentStrainID in unusedExcludeStrains:
                            unusedExcludeStrains.remove(currentStrainID)

                currentStrainID = matchStrainID
                currentStrainLociList = ["{0}\t{1}".format(matchLocus,matchBase)]

            else:
                currentStrainLociList.append("{0}\t{1}".format(matchLocus,matchBase))


    # Process last strain from loop.
    if (currentStrainID != "") and (currentStrainLociList != ""):

        if currentStrainID in compareList:
            snpCompareTable[currentStrainID] = set(currentStrainLociList)

            if currentStrainID in unusedCompareStrains:
                unusedCompareStrains.remove(currentStrainID)
        else:
            snpExcludeTable[currentStrainID] = set(currentStrainLociList)

            if currentStrainID in unusedExcludeStrains:
                unusedExcludeStrains.remove(currentStrainID)

finally:
    snpFile.close()

print("Comparison Strains: {0}\n".format(', '.join(compareList)))
print("Exclusion Strains: {0}\n\n".format(', '.join(excludeList)))

if len(unusedCompareStrains) > 0:
    print("*** WARNING: Some comparison strains were not found in snp file: {0}".format(', '.join(unusedCompareStrains)))

if len(unusedExcludeStrains) > 0:
    print("*** WARNING: Some exclusion strains were not found in snp file: {0}".format(', '.join(unusedExcludeStrains)))

print("Generating comparison report...")
try:
    outFile = open(outFilename,"w")
    outFile.write("SNP Comparator Report\n")
    outFile.write("Comparison Strains: {0}\n".format(', '.join(compareList)))
    outFile.write("Exclusion Strains: {0}\n".format(', '.join(excludeList)))
    outFile.write("\n")
    outFile.write("Common SNP Loci (after exclusion, if any):\n")

    # Create the set of intersections of SNP loci from the comparison strains.
    commonSet = set()
    exclusionSet = set()
    for strainSet in snpCompareTable.values():
        if len(commonSet) == 0:
            commonSet = strainSet
        else:
            commonSet = commonSet.intersection(strainSet)

    # Create the set of intersections of SNP loci from the exclusion strains.
    for strainSet in snpExcludeTable.values():
        if len(exclusionSet) == 0:
            exclusionSet = strainSet
        else:
            exclusionSet = exclusionSet.intersection(strainSet)

    # Create the set of the difference between the common set and exclusion set
    # -- i.e. exclude the exclusion set snp loci from the common set.
    finalSet = commonSet.difference(exclusionSet)

    # Write out the final common, filtered SNP loci.
    # Split the locus and base up into distinct elements for sorting. Note we get a list of lists like [locus, base].
    finalList = [ [x.split("\t")[0],x.split("\t")[1]] for x in finalSet ]
    finalList.sort()

    for snpLoci in finalList:
        outFile.write("{0}\t{1}\n".format(snpLoci[0],snpLoci[1]))

    outFile.write("\nCount of common SNP loci from comparison set: {0}\n".format(len(commonSet)))
    outFile.write("Count of common SNP loci from exclusion set: {0}\n".format(len(exclusionSet)))
    outFile.write("FINAL count of common SNP loci (after exclusions): {0}\n".format(len(finalSet)))

    if len(unusedCompareStrains) > 0:
        outFile.write("\n*** WARNING: Some comparison strains were not found in snp file: {0}".format(', '.join(unusedCompareStrains)))

    if len(unusedExcludeStrains) > 0:
        outFile.write("\n*** WARNING: Some exclusion strains were not found in snp file: {0}".format(', '.join(unusedExcludeStrains)))

    outFile.write("\n")
finally:
    outFile.close()


print("\nDone.\n")
print("Count of common SNP loci from comparison set: {0}".format(len(commonSet)))
print("Count of common SNP loci from exclusion set: {0}".format(len(exclusionSet)))
print("FINAL count of common SNP loci (after exclusions): {0}\n".format(len(finalSet)))


if len(unusedCompareStrains) > 0:
    print("*** WARNING: Some comparison strains were not found in snp file: {0}".format(', '.join(unusedCompareStrains)))

if len(unusedExcludeStrains) > 0:
    print("*** WARNING: Some exclusion strains were not found in snp file: {0}".format(', '.join(unusedExcludeStrains)))

print("Output file is {0}\n".format(outFilename))

