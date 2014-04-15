#!/usr/bin/python
#
# Prephix 2 SNP Effect Converter
# by Andrew Pann
#
# This script takes the ref and snp output files from a prephix run and generates an SNP Effect compliant input file.
# SNP Effect file format is: SNPlocus<tab>sample=snpbase<tab>ref=refbase
# SNP Effect doesn't care about strain ids so if two strains share a locus but have different bases, then two lines will
# be generated, with identical locus and ref values but different sample values.
#
# Usage: prephix2snpeffect.py ref_file snp_file
#
# 7/31/2013 -- VERSION 2.0.0 - Moved to Sqlite3 datastore and query processing, instead of dictionaries.
# 7/31/2013 -- VERSION 2.0.1 - Add better error output on database insert failure.
# 8/21/2013 -- VERSION 2.1.0 - Removed reporting of duplicate bases at same locus.
# 12/28/2013 -- VERSION 2.2.0 - Added logic to ignore inserts of duplicate strain id/locus/base during SNP processing.
# 4/13/2014 -- VERSION 2.2.1 - Fixed minor bug in reporting of empty strain inputs.

import sys
import os
import re
import argparse
import string
import cStringIO
import math
import sqlite3

VERSION = '2.3.1'

# MAIN #

print "\nPrephix to SNP Effect Converter v{0}\n".format(VERSION)

# Define our parameters and parse with argparse module.
argParser = argparse.ArgumentParser()

progname = os.path.abspath(sys.argv[0])

# Define arguments.
argParser.add_argument("ref_file",help="The prephix written reference output file.")
argParser.add_argument("snp_file",help="The prephix written snp file output.")
argParser.add_argument("--dbfile",nargs=1,help="The SQLite database filename. Optional and useful for debugging.")
argParser.add_argument("-d","--debug",action="store_true",help="Debug mode.")

args = argParser.parse_args()

# Set debug mode.
debugMode = args.debug

# Working variables.
refFilename = args.ref_file
snpFilename = args.snp_file
outputFilename = refFilename + ".snpeffect"

if debugMode:
    print "Argparse results: "
    print args


try:
    snpFile = open(snpFilename,"r")
except IOError as e:
    print "Unable to open file " + snpFilename
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit(1)
except:
    print "Unexpected error while opening SNP file:", sys.exc_info()[0]
    raise

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

# Create the table for the snp file data (COLUMNS: strainid, locus, base)
dbcursor.execute('''CREATE TABLE SNP_DATA (strainid text, locus integer, base text NOT NULL, PRIMARY KEY(strainid,locus))''')

# Create indexes.
dbcursor.execute('''CREATE INDEX SNP_DATA_LOCUS_IDX ON SNP_DATA (locus)''')

dbconn.commit()

print "Reading SNP file data: {0}...".format(snpFilename)
sys.stdout.flush()

with dbconn:
    for snpLine in snpFile:
        # Skip strains with no SNPs (will have locus value -1)
        snpBadLocusMatch = re.match("^(?P<strainid>[^\t]+)\t-1",snpLine)
        if snpBadLocusMatch:
            print "Ignoring entry for strain {} due to locus value of -1 (No SNPs).".format(str(snpBadLocusMatch.group("strainid")))
            continue

        # Expect the prephix snp output file to be in format: STRAIN_ID\tLOCUS\tSNP_BASE
        snpLineMatch = re.match("^(?P<strainid>[^\t]+)\t(?P<locus>\d+)\t(?P<snpBase>[ACGT])$",snpLine)
        if snpLineMatch:
            # Database implementation here - store the SNP information in the SNP_DATA table.
            strainid = snpLineMatch.group("strainid")
            locus = snpLineMatch.group("locus")
            snpBase = snpLineMatch.group("snpBase")

            try:
                dbconn.execute('''INSERT INTO SNP_DATA (strainid,locus,base) VALUES (?,?,?)''',(strainid,locus,snpBase))
            except sqlite3.IntegrityError as e:
                # IntegrityError can be caused by non-unique straid and locus.  This is ok if the BASE is the same, however.
                # Ignore if this is a true duplicate.  Else it is an error.
                dbcursor.execute('''SELECT base FROM SNP_DATA WHERE STRAINID= ? AND LOCUS = ?''',(strainid,locus))
                for row in dbcursor.fetchall():
                    selectBase = str(row[0])

                if selectBase == snpBase:
                    print "Ignoring duplicate base {} at locus {} for strain {}".format(snpBase,locus,strainid)
                else:
                    print "*** DATABASE error inserting {} {} {}".format(strainid,locus,snpBase)
                    print "Existing base at this locus for this strain is {}.  Tried to insert different base {} at same locus for same strain id!".format(selectBase,snpBase)
                    raise 
            except:
                print "*** DATABASE error inserting {} {} {}".format(strainid,locus,snpBase)
                raise 

        else:
            print "*** ERROR: Bad line in {0}! Not a prephix SNP file?".format(snpFilename)
            print "Cannot parse line: {0}".format(snpLine)
            sys.exit(1)

snpFile.close()


print "Reading ref file data and generating output file: {0}".format(outputFilename)
sys.stdout.flush()

try:
    refFile = open(refFilename,"r")
except IOError as e:
    print "Unable to open file " + refFilename
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit(1)
except:
    print "Unexpected error while opening ref file:", sys.exc_info()[0]
    raise

try:
    outFile = open(outputFilename,"w")
except IOError as e:
    print "Unable to open file " + outputFilename
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit(1)
except:
    print "Unexpected error while opening snp effect output file:", sys.exc_info()[0]
    raise

dupTracker = {} # Variable to track locus/base pairs to skip duplicate values.

for refLine in refFile:
    # Reference file generated by prephix should be in format LOCUS\tREF_BASE.
    refLineMatch = re.match("^(?P<locus>\d+)\t(?P<refBase>[ACGT])$",refLine)
    if refLineMatch:
        locus = refLineMatch.group("locus")
        refBase = refLineMatch.group("refBase")
        dbcursor.execute('''SELECT base FROM SNP_DATA WHERE LOCUS = ?''',(locus,))

        got_results = False

        for row in dbcursor.fetchall():
            got_results = True
            snpBase = str(row[0])

            # SNP Effect file format has a -1 offset from actual locus value.
            snpEffectLocus = int(locus) - 1
            
            dupKey = "{0}{1}".format(snpEffectLocus,snpBase)
            if not dupKey in dupTracker:
                # SNP Effect file format is: SNPlocus<tab>sample=snpbase<tab>ref=refbase
                outFile.write("{0}\tsample={1}\tref={2}\n".format(snpEffectLocus,snpBase,refBase))
                dupTracker[dupKey] = "Y"


        if not got_results:
            print "*** ERROR: No SNP base found for locus {0} in ref file!".format(locus)
            sys.exit(1)

    else:
        print "*** ERROR: Bad line in {0}! Not a prephix ref file?".format(refFilename)
        print "Cannot parse line: {0}".format(refLine)
        sys.exit(1)


refFile.close()
outFile.close()

dbconn.close()

dupCount = len(dupTracker)
print "Duplicate locus/pair entries removed: {0}".format(dupCount)
print "Done.  SNP Effect output file is {0}".format(outputFilename)
