#!/usr/bin/python
#
# Prephix (Pre-Phrecon Input fiXer)
#
# Takes in a list of input files and generates one SNP loci output file (containing multiple strains), one reference base file,
# and a file containing the excluded indel entries.  It generates the reference base sequence based on the SNP loci information.
# Duplicate SNP loci with different bases is considered an error.
#
# The indel file is formated as "STRAIN ID" [TAB] k28|nuc|vcf [TAB] ...excluded VAAL K28, VCF, or NUCMER line...
# So prephix basically dumps the excluded lines as-is into the indel file, but prefixes each line
# with the strain id and then either k28 for VAAL formatted line, nuc for NUCMER formatted line, or vcf for VCF formatted line.
#
# The SNP loci and reference base files are suitable as input to the phrecon (Phylo Reconstructor) utility.
# SNP loci file has the tab-delimited columnar format: STRAIN_ID [TAB] LOCI [TAB] BASE
# Reference base file has the tab-delimited column format: LOCI [TAB] BASE
#
# The SNP loci file and indel files are suitable input for the SNP Swapper utility.
#
# Recognized input file formats: k28.out (VAAL), NUCMER, VCF
#
# You may mix input file types within a single run/batch.
#

import sys
import os
import re
import stat
import argparse
import logging

import statprof

# Custom include
import SNPInputReader

VERSION = '3.5.0'



#####################
# UTILITY FUNCTIONS
#####################

def print_all(msg):
    '''Prints output to STDOUT and logfile in one call.
       If the quiet funnction is set, then only log to file.

       Parameters:
       msg = Text to print/log.
    '''
    if not quietMode:
        print msg
    logfile.write("{}\n".format(msg))


########
# MAIN #
########
if __name__ == '__main__':
    print "\nPrephix (Pre-Phrecon Input fiXer) v{0}\n".format(VERSION)

    # Define our parameters and parse with argparse module.
    argParser = argparse.ArgumentParser()

    progname = os.path.abspath(sys.argv[0])

    # Define arguments.
    argParser.add_argument("input_files",nargs='*',help="One or more input files. Supported input file formats are k28.out (VAAL), NUCMER, VCF.")
    argParser.add_argument("-batchid","--batchid",required=True,help="Used to create the output filenames for the combined SNP loci and ref files.")
    argParser.add_argument("-exclude","--exclude",metavar='loci_exclusion_file',help="Exclude any bases within the loci ranges listed in the loci exclusion file.  The file should contain lines of the format 'label,start_loci,end_loci'")
    argParser.add_argument("-filter_quality","--filter_quality",action="store_true",help="Exclude any bases not passing the quality filter in VCF files.",default=False)
    argParser.add_argument("-tablog","--tablog",action="store_true",help="Print out the summary in tabular format.")
    argParser.add_argument("-export_phenolink","--export_phenolink",action="store_true",help="Additionally write out a Pheonolink compatible output file.")
    argParser.add_argument("-debug","--debug",action="store_true",help="Debug mode.")
    argParser.add_argument("-multichrom","--multichrom",action="store_true",help="VCF and NUCMER files have multi-chrom values.")
    argParser.add_argument("-quiet","--quiet",action="store_true",help="Quiet mode.")


    # Process arguments.
    args = argParser.parse_args()

    # Set argument variables.
    debugMode = args.debug
    if debugMode:
        print "Producing debug output."
        print "Argparse results: "
        print args
    
    quietMode = args.quiet
    if quietMode:
        print "Producing quiet (no stdout) output. Log file is still generated."

    excludeFileName = args.exclude
    if excludeFileName:
        print "Exclusion file is {}".format(excludeFileName)

    batchid = args.batchid
    print "Batch id is {}".format(batchid)

    filterQuality = args.filter_quality
    if not filterQuality:
        print "Will process all lines, ignoring quality value (only applicable for VCF input files)."
    else:
        print "Will processes only lines passing quality value (for VCF input files)."

    writeTablog = args.tablog
    if writeTablog:
        print "Will write summary in tabular format."

    exportPhenoLink = args.export_phenolink
    if exportPhenoLink:
        print "Will export a PhenoLink file."

    multiChrom = args.multichrom
    if multiChrom:
        print "Will use multi-FASTA (multi-chrom) parsing on VCF and NUCMER input files."

    inputFileList = args.input_files
    fileCount = 0
    for inputFilename in inputFileList:
        fileCount += 1 
        if os.path.isfile(inputFilename):
            print "Found {} file to process.".format(inputFilename)
        else:
            print "*** ERROR: File does not exist: {}".format(inputFilename)
            sys.exit(1)
            


    # Prep files for read/write.

    logFilename = "{}.log".format(batchid)
    logfile = open(logFilename,"w")

    print_all("Found {} files to process.".format(fileCount))

    debugLogFilename = "{}.debug.log".format(batchid)
    if debugMode:
        debugLevel = logging.DEBUG
    else:
        debugLevel = logging.WARNING
    logging.basicConfig(filename=debugLogFilename,filemode='w',level=debugLevel)

    # SNP output file -- SNP loci file format is (StrainId [TAB] Loci [TAB] Base)
    outFilename = "{}.snp".format(batchid)
    outfile = open(outFilename,"w")

    # Reference base file -- file format is (Loci [TAB] Base)
    refFilename = "{}.ref".format(batchid)
    reffile = open(refFilename,"w")

    # The indel file format is a tab-delimted line of strain id, file format, the raw line, and indel type.
    #
    # NOTE: Compared to version 2 of prephix (perl), this has an additional field at the end, which is INS or DEL to
    # indicate if the indel entry is an insertion or a deletion.
    indelFilename = "{}.indel".format(batchid)
    indelfile = open(indelFilename,"w")

    # Setup data dictionaries

    # Reference Base Data Table
    #
    # Key is locus integer having a value of a 3-tuple in the order (base, source file name, line number)
    # So {locus: (base, source filename, line number) }
    # The source file name and line number are for tracking the original source of the base.
    refDataTable = {} 
    
    # Exclusion Reporting Data Table
    #
    # Key is strain id, with a value of another dictionary with key exclude_label having a value of integer count.
    # So {strainID: {exclude_label: count} }
    excludeDataTable = {}

    # Overall Stats Report Data Table
    #
    # Key is strain id having value of a list of SNP count, insert count, delete counts, exclusion count.
    # So {strainID: [# SNPs, # insertions, # deletions, # exclusions] }
    statsTable = {}

    # Read exclusion file, if any.
    # Expect exclusion file to contain one loci range per line, formatted: label,start_loci,end_loci (inclusive)
    exclusionTable = {}
    if excludeFileName:
        excludefile = open(excludeFileName,"r")
        excludeCount = 0
        excludeRe = re.compile("^(?P<label>[^,]+),(?P<start_loci>[0-9]+),(?P<end_loci>[0-9]+)$")

        for excludeLine in excludefile:
            excludeCount += 1
            # Store exclusions in hash with label as key and start,end loci as array values.
            excludeMatch = excludeRe.match(excludeLine)

            if excludeMatch:
                excludeLabel = excludeMatch.group('label')
                excludeStartLoci = excludeMatch.group('start_loci')
                excludeEndLoci = excludeMatch.group('end_loci')

                if excludeLabel in exclusionTable:
                    print_all("Duplicate exclusion label {} encountered on line {}. Quitting!".format(excludeLabel,excludeCount))
                    sys.exit(1)
                            
                exclusionTable[excludeLabel] = [excludeStartLoci,excludeEndLoci]
                logging.debug("Adding exclusion: {},{},{}".format(excludeLabel,excludeStartLoci,excludeEndLoci)) 
            else:
                print_all("Badly formatted line: {} at line {}. Quitting!".format(excludeLine,excludeCount))
                sys.exit(1)

        excludefile.close()
        print_all("{} exclusions read.".format(excludeCount))

    print_all("\n***");
    print_all("*** REMINDER: This program assumes that all input files refer to the same reference sequence.");
    print_all("***           Ref file will be generated from consolidated snp loci information of ALL input files.");
    print_all("***");


    fileCount = 0
    skippedFileCount = 0
    # Process the snp files.
    while len(inputFileList) > 0:

        inputFilename = inputFileList.pop()

        fileCount += 1
        try:
	    snpFileReader = SNPInputReader.getSNPFileReader(inputFilename,filterQuality, multiChrom)
        except SNPInputReader.EmptyFileError as efe:
            print_all("*** WARNING: Skipping empty file: {0}".format(inputFilename))
            logging.debug("*** WARNING: Skipping empty file: {0}".format(inputFilename))
            skippedFileCount = skippedFileCount + 1
            continue

        strainid = snpFileReader.strainID
        fileFormat = snpFileReader.fileFormat
        shortFilename = os.path.basename(inputFilename)

        print_all("Processing {} file {}...".format(fileFormat,shortFilename))
        sys.stdout.flush()

        logging.debug("Strain ID is %s",strainid)

        # Initialize stats entry for this strain.
        statsTable[strainid] = [0,0,0,0]

        for snpData in snpFileReader:
            rawLine = snpData[0]
            lineNumber = snpData[1]
            locus = snpData[2]
            snpBase = snpData[3]
            refBase = snpData[4]
            isIndel = snpData[5]
            isInsert = snpData[6]
            isDelete = snpData[7]

            logging.debug("At line: %s",lineNumber)

            # Record indels into their table, and skip futher processing.
            if isIndel:
                if isInsert:
                    logging.debug("Found indel line (insertion): %s (%s)",rawLine,shortFilename)
                    indelType = "INS"

                    # Increment insertion count in report
                    statsTable[strainid][1] += 1
                elif isDelete:
                    logging.debug("Found indel line (deletion): %s (%s)",rawLine,shortFilename)
                    indelType = "DEL"

                    # Increment deletion count in report
                    statsTable[strainid][2] += 1
                else:
                    print_all("*** ERROR: Unknown indel type found at line {}: {}".format(lineNumber,rawLine))
                    sys.exit(1)
        
                # Write entry to indel file. Format is tab-delimted line of strain id, file format, the raw line, and indel type.
                # Compared to version 2 of prephix (perl), this has an additional field at the end, which is INS or DEL to
                # indicate if the indel entry is an insertion or a deletion.
                indelfile.write("{}\t{}\t{}\t{}\n".format(strainid,fileFormat,rawLine,indelType))

            # Check exclusion of locus.
            excluded = False
            for excludeLabel in exclusionTable:
                excludeStartLocus = int(exclusionTable[excludeLabel][0])
                excludeEndLocus = int(exclusionTable[excludeLabel][1])
                #logging.debug("Exclusion test: is {} between {} and {}".format(locus,excludeStartLocus,excludeEndLocus))
                if locus >= excludeStartLocus and locus <= excludeEndLocus:
                    excluded = True
                    logging.debug("Excluded loci %s",locus)

                    # Report exclusion count.
                    statsTable[strainid][3] += 1

                    if not strainid in excludeDataTable:
                        excludeDataTable[strainid] = {excludeLabel: 1}
                    elif not excludeLabel in excludeDataTable[strainid]:
                        excludeDataTable[strainid][excludeLabel] = 1
                    else:
                        excludeDataTable[strainid][excludeLabel] += 1

            # If SNP data was an indel or excluded, skip further processing. 
            if excluded or isIndel:
                logging.debug("Data was excluded or indel, so skipping...")
                continue


            # Write data to SNP output file.
            # SNP loci file format is (StrainId [TAB] Loci [TAB] Base)
            outfile.write("{}\t{}\t{}\n".format(strainid,locus,snpBase))

            # Record SNP count.
            statsTable[strainid][0] += 1

            # Insert data into REF dictionary
            # First check if this is a collision with existing reference base data at same locus (but different base).
            if locus in refDataTable:
                if refDataTable[locus][0] != refBase:
                    print_all("*** ERROR: Reference base mismatch at loci {}! Input file {} line {} has ref={}, but ref base at this loci was already recorded as {} while processing file {} line {}!".format(locus,shortFilename,lineNumber,refBase,refDataTable[locus][0],refDataTable[locus][1],refDataTable[locus][2]))
                    print_all("*** Are you sure all input files are from the same reference?")
                    print_all("Failed.")
                    sys.exit(1)
                else:
                    logging.debug("Duplicate (but identical - so this is OK) ref base found at locus %s: %s",locus,refBase)
            else:
                # If no collision, add it to the reference data table.
                refDataTable[locus] = (refBase,shortFilename,lineNumber)
                logging.debug("Added ref locus %s: (%s,%s,%s)",locus,refBase,shortFilename,lineNumber)

        # If snp input file was empty, write a placeholder for it.
        if statsTable[strainid][0] == 0:
            logging.debug("Empty SNP file for strain {}".format(strainid))
            outfile.write("{}\t{}\t{}\n".format(strainid,-1,"-"))


    # Write out the reference file from the table of merged ref loci bases from the input file.
    # Output format is Loci [TAB] Base
    print_all("")
    print_all("{} files were processed.".format(fileCount))
    print_all("")
    print_all("WARNING --> {} files were SKIPPED (empty?).".format(skippedFileCount))
    print_all("")
    print_all("Merging and generating reference file from input file data....")
    refDataTableKeyList = refDataTable.keys()
    refDataTableKeyList.sort()
    for refLocus in refDataTableKeyList:
        reffile.write("{}\t{}\n".format(refLocus,refDataTable[refLocus][0]))


    print_all("Done.")

    # Cleanup
    outfile.close()
    reffile.close()

    if exportPhenoLink:
        print_all("Exporting PhenoLink file...")
        # For now, assuming pre2phe.py is in the same working directory.
        if os.path.exists("./pre2phe.py"):
            scriptStat = os.stat("./pre2phe.py")
            if (stat.S_IXUSR & scriptStat.st_mode) or (stat.S_IXGRP & scriptStat.st_mode) or (stat.S_IXOTH & scriptStat.st_mode):
                os.system("./pre2phe.py --ref {} --snp {} --out {}.phenonlink.txt".format(refFilename,outFilename,batchid))
            else:
                print_all("***\n*** ERROR: Cannot find or execute helper script pre2phe.py in working directory!\n***")
                print_all("Unable to generate PhenoLink file!")


    print_all("Merged SNP loci file is {}".format(outFilename))
    print_all("Merged reference base file from this run is {}".format(refFilename))
    print_all("Merged indel file from this run is {}".format(indelFilename))

    # Write out stats report.
    print_all("\n=== Final Report ===\n")

    if not writeTablog:
        for strainid in statsTable:
            print_all("Strain: {}".format(strainid))
            print_all("SNPs: {}".format(statsTable[strainid][0]))
            print_all("Insertions: {}".format(statsTable[strainid][1]))
            print_all("Deletions: {}".format(statsTable[strainid][2]))

            totalIndels = statsTable[strainid][1] + statsTable[strainid][2]

            print_all("Total indels: {}".format(totalIndels))
            print_all("Loci excluded: {}".format(statsTable[strainid][3]))

            if strainid in excludeDataTable:
                for excludeLabel in excludeDataTable[strainid]:
                    print_all("* Loci excluded from {}: {}".format(excludeLabel,excludeDataTable[strainid][excludeLabel]))
            print_all("")

        print_all("")

    else:
        # Tabular format requsted by -tablog flag.
        print_all("Strain ID\tSNPs\tInserts\tDeletes\tTotal Indels\tLoci Excluded\n")
        for strainid in statsTable:
            totalIndels = statsTable[strainid][1] + statsTable[strainid][2]
            print_all("{}\t{}\t{}\t{}\t{}\t{}".format(strainid,statsTable[strainid][0],statsTable[strainid][1],statsTable[strainid][2],totalIndels,statsTable[strainid][3]))

        print_all("")

        for strainid in statsTable:
            if strainid in excludeDataTable:
                print_all("Loci exclusion summary for {}:".format(strainid))
                for excludeLabel in excludeDataTable[strainid]:
                    print_all("* Loci excluded from {}: {}".format(excludeLabel,excludeDataTable[strainid][excludeLabel]))

    print_all("\nDone.\n")
