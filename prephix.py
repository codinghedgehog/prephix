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
# 8/08/2013 - Andrew Pann - Initial port of perl-version Prephix (v2.5.2) to 3.0.


import sys
import os
import re
import argparse
import string
import math
import sqlite3

# Custom include
import SNPInputReader

VERSION = '3.0.0'

def main():
    print "\nPrephix (Pre-Phrecon Input fiXer) v{0}\n".format(VERSION)

    # Define our parameters and parse with argparse module.
    argParser = argparse.ArgumentParser()

    progname = os.path.abspath(sys.argv[0])

    # Define arguments.
    argParser.add_argument("input_files",nargs='*'help="One or more input files. Supported input file formats are k28.out (VAAL), NUCMER, VCF.")
    argParser.add_argument("-batchid","--batchid",required=True,help="Used to create the output filenames for the combined SNP loci and ref files.")
    argParser.add_argument("--dbfile",nargs=1,help="The database filename.")
    argParser.add_argument("-exclude","--exclude",metavar='loci_exclusion_file',nargs=1,help="Exclude any bases within the loci ranges listed in the loci exclusion file.  The file should contain lines of the format 'label,start_loci,end_loci'")
    argParser.add_argument("-ignore_quality","--ignore_quality",action="store_true",help="Exclude any bases within the loci ranges listed in the loci exclusion file.  The file should contain lines of the format 'label,start_loci,end_loci'")
    argParser.add_argument("-tablog","--tablog",action="store_true",help="Print out the summary in tabular format.")
    argParser.add_argument("-export_phenolink","--export_phenolink",action="store_true",help="Additionally write out a Pheonolink compatible output file.")
    argParser.add_argument("-debug","--debug",action="store_true",help="Debug mode.")
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
        print "Exclusion file is {}".format(excludeFilename)

    batchid = args.batchid
    print "Batch id is {}".format(batchid)

    filterQuality = not args.ignore_quality
    if not filterQuality:
        print "Will process all lines, ignoring quality value (only applicable for VCF input files)."

    writeTablog = args.tablog
    if writeTablog:
        print "Will write summary in tabular format."

    exportPhenoLink = args.export_phenolink
    if exportPhenoLink:
        print "Wille xport a PhenoLink file."

    inputFileList = args.input_files
    for inputFilename in inputFileList:
        if os.path.isfile(inputFilename):
            print "Found {} file to process.".format(inputFilename)
        else:
            print "*** ERROR: File does not exist: {}".format(inputFilename)
            sys.exit(1)
            

    # Prep files for read/write.
    logFilename = "{}.log".format(batchid)
    logfile = open(logFilename,"w")

    outFilename = "{}.snp".format(batchid)
    outfile = open(outFilename,"w")

    refFilename = "{}.ref".format(batchid)
    reffile = open(refFilename,"w")

    indelFilename = "{}.indel".format(batchid)
    indelfile = open(indelFilename,"w")

    # Setup in-memory database.
    if args.dbfile == None:
        dbconn = sqlite3.connect(':memory:')
    else:
        dbfilename = args.dbfile[0]
        if os.path.isfile(dbfilename):
            print "*** ERROR database file {0} already exists.  Please remove or use a different filename.\n".format(dbfilename)
            sys.exit(1)
        else:
            dbconn = sqlite3.connect(dbfilename)
    dbcursor = dbconn.cursor()

    # Create the table for the ref file data (COLUMNS: locus, base, source_file, line_number)
    dbcursor.execute('''CREATE TABLE REF_DATA (locus integer PRIMARY KEY ASC, base text NOT NULL, source_file text, line_number integer)''')

    # Create the table for the snp file data (COLUMNS: strainid, locus, base)
    dbcursor.execute('''CREATE TABLE SNP_DATA (strainid text, locus integer, base text NOT NULL, PRIMARY KEY(strainid,locus))''')

    # Create the table for the indel file data (COLUMNS: strain,format, indel_type, raw line).
    # indel_type values should be one of INS or DEL.
    dbcursor.execute('''CREATE TABLE INDEL_DATA (strainid text NOT NULL, format text NOT NULL, indel_type text NOT NULL, rawline text NOT NULL)''')

    # Create the table for the exclusion file report (COLUMNS: strainid, exclude_label, locus, raw line).
    # indel_type values should be one of INS or DEL.
    dbcursor.execute('''CREATE TABLE EXCLUSION_DATA (strainid text NOT NULL, exclude_label text NOT NULL, locus integer NOT NULL, rawline text NOT NULL)''')

    # Create indexes.
    dbcursor.execute('''CREATE INDEX SNP_DATA_LOCUS_IDX ON SNP_DATA (locus)''')
    dbcursor.execute('''CREATE INDEX INDEL_DATA_STRAINID_IDX ON SNP_DATA (strainid)''')
    dbcursor.execute('''CREATE INDEX EXCLUSION_DATA_STRAINID_LABEL_IDX ON SNP_DATA (strainid,exclude_label)''')

    dbconn.commit()

    # Read exclusion file, if any.
    # Expect exclusion file to contain one loci range per line, formatted: label,start_loci,end_loci (inclusive)
    if excludeFileName:
        exclusionTable = {}
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

                if excludeLabel in excludeTable:
                    print_all("Duplicate exclusion label {} encountered on line {}. Quitting!".format(excludeLabel,excludeCount))
                    sys.exit(1)
                            
                exclusionTable[excludeLabel] = [excludeStartLoci,excludeEndLoci]
            else
                print_all("Badly formatted line: {} at line {}. Quitting!".format(excludeLine,excludeCount))
                sys.exit(1)

        excludefile.close()
        print_all("{} exclusions read.".format(excludeCount))

    print_all("\n***");
    print_all("*** REMINDER: This program assumes that all input files refer to the same reference sequence.");
    print_all("***           Ref file will be generated from consolidated snp loci information of ALL input files.");
    print_all("***");


    fileCount = 0
    # Process the snp files.
    for inputFilename in inputFileList:
        fileCount += 1
        snpFileReader = SNPInputReader.getSNPFileReader(inputFilename,filterQuality)
        strainid = snpFileReader.strainID
        fileFormat = snpFileReader.format
        shortFilename = os.path.basename(inputFilename)

        print_all("Processing {} file {}...".format(fileFormat,shortFilename))
        sys.stdout.flush()

        print_debug("Strain ID is {}".format(strainid))

        for snpData in snpFileReader:
            locus = snpData.locus
            snpBase = snpData.snpBase
            refBase = snpData.refBase

            # Record indels into their table, and skip futher processing.
            if snpData.isIndel:
                if snpData.isInsert:
                    print_debug("Found indel line (insertion): {} ({})".format(snpData.rawLine,shortFilename))
                    indelType = "INS"
                elif snpData.isDelete:
                    print_debug("Found indel line (deletion): {} ({})".format(snpData.rawLine,shortFilename))
                    indelType = "DEL"
                else:
                    print_all("*** ERROR: Unknown indel type found at line {}: {}".format(snpData.lineNumber,snpData.rawLine))
                    sys.exit(1)
        
                try:
                    dbconn.execute('''INSERT INTO INDEL_DATA (strainid,format,indel_type,rawline) VALUES (?,?,?,?)''',(strainid,fileFormat,indelType,snpData.rawLine))
                except:
                    print_all("*** DATABASE error inserting into INDEL_DATA {} {} {} {}".format(strainid,fileFormat,indelType,snpData.rawLine))
                    raise 

            # Check exclusion of locus.
            excluded = False
            for excludeLabel in excludeTable:
                excludeStartLocus = excludeTable[excludeLabel][0]
                excludeEndLocus = excludeTable[excludeLabel][1]
                if locus >= excludeStartLocus and locus <= excludeEndLocus:
                    excluded = True
                    print_debug("Excluded loci {}".format(str(locus)))
                    try:
                        dbconn.execute('''INSERT INTO EXCLUSION_DATA (strainid,exclude_label,locus,rawline) VALUES (?,?,?,?)''',(strainid,excludeLabel,locus,snpData.rawLine))
                    except:
                        print_all("*** DATABASE error inserting into EXCLUSION_DATA {} {} {} {}".format(strainid,excludeLabel,locus,snpData.rawLine))
                        raise 
                    break

           
            # If SNP data was an indel or excluded, skip further processing. 
            if excluded or snpData.isIndel:
                continue
            else:
                # Insert data into SNP and REF tables.
                try:
                    dbconn.execute('''INSERT INTO SNP_DATA (strainid,locus,base) VALUES (?,?,?)''',(strainid,locus,snpBase))
                except:
                    print_all("*** DATABASE error inserting into SNP_DATA {} {} {}".format(strainid,locus,snpBase))
                    raise 

                try:
                    dbconn.execute('''INSERT INTO REF_DATA (locus,base,source_file,line_number) VALUES (?,?,?)''',(locus,refBase,shortFilename,snpData.lineNumber))
                except:
                    print_all("*** DATABASE error inserting into REF_DATA {} {} {}".format(locus,refBase,shortFilename,snpData.lineNumber))
                    print_all("*** Are you sure all input files are using the same reference?")
                    print_all("Failed.")
                    raise 

#####################
# UTILITY FUNCTIONS
#####################

def print_debug(msg):
    '''Prints output to STDOUT and logfile if debug flag is set.
       If the quiet funnction is also set, then only log to file.

       Parameters:
       msg = Text to print/log.
    '''
    if debugMode:
        if not quietMode:
            print msg
        logfile.write("{}\n".format(msg))

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
    main()
