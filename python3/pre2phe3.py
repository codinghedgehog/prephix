#!/usr/bin/python3
#
# Prephix to PhenoLink Generator
#
# Usage: pre2phe.py -ref <ref base file> -snp <snp file> -out <PhenoLink output file>
#
# Takes in a one SNP loci output file (possibly containing multiple strains) and one reference base file generated from prephix,
# and generates a PhenoLink export file.
#
# The PhenoLink export file is a tab-delimited file, with columns of all Loci Names (in this case RefBase_Locus_SnpBase)
# and rows of StrainIDs with 1's indicating the presence of that Locus Name in its strain, and 0 otherwise.
#
# 7/31/2013 - VERSION 2.0 - Moved to a database backed data store and query process.

VERSION="2.0"

import argparse
import re
import os
import sys
import sqlite3

print("\nPrephix to PhenoLink Generator v{0}\n\n".format(VERSION))

# Setup argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('--ref',required=True,type=str, nargs=1,help="Reference base output file from prephix.")
parser.add_argument('--snp',required=True,type=str,nargs=1,help="SNP output file from prephix.")
parser.add_argument('--out',required=True,type=str,nargs=1,help="PhenoLink output file name.")
parser.add_argument('--dbfile',type=str,nargs=1,help="PhenoLink database name (created new for each run).")
args = parser.parse_args()

reffile = open(args.ref[0],'r')
snpfile = open(args.snp[0],'r')
outfile = open(args.out[0],'w')

# Setup in-memory (or disk, if --dbfile is used) database.
if args.dbfile != None:
    dbfilename = args.dbfile[0]
    if os.path.isfile(dbfilename):
        print("*** ERROR database file {0} already exists.  Please remove or use a different filename.\n".format(dbfilename))
        sys.exit(1)
    else:
        dbconn = sqlite3.connect(dbfilename)
else:
    dbconn = sqlite3.connect(':memory:')

dbcursor = dbconn.cursor()

# Create the table for the ref file data (COLUMNS: locus, base)
dbcursor.execute('''CREATE TABLE REF_DATA (locus integer PRIMARY KEY ASC, base text NOT NULL)''')

# Create the table for the snp file data (COLUMNS: strainid, locus, base)
dbcursor.execute('''CREATE TABLE SNP_DATA (strainid text, locus integer, base text NOT NULL, PRIMARY KEY(strainid,locus))''')

# Create the table for the PhenoLink file data (COLUMNS: locus_name,strainid)
dbcursor.execute('''CREATE TABLE PHENOLINK_DATA (strainid text, locus_name text, PRIMARY KEY(strainid,locus_name))''')

# Create indexes.
dbcursor.execute('''CREATE INDEX SNP_DATA_LOCUS_IDX ON SNP_DATA (locus)''')
dbcursor.execute('''CREATE INDEX PHENOLINK_DATA_LOCUS_NAME_IDX ON PHENOLINK_DATA (locus_name)''')
dbcursor.execute('''CREATE INDEX PHENOLINK_DATA_STRAINID_IDX ON PHENOLINK_DATA (strainid)''')

dbconn.commit()


# Read the reference base file into the database.
refTable = dict()
lineNumber=1
print("Reading base reference file..." ,)
sys.stdout.flush()

try:
    with dbconn:
        refRegex = re.compile("^(?P<locus>\d+)\t(?P<base>\D)$")
        for line in reffile:
            refMatch = refRegex.match(line)
            if refMatch == None:
                print("*** WARNING: Unrecognized format at line {0}: {1}".format(lineNumber,line))
            else:
                locus = refMatch.group('locus')
                base = refMatch.group('base')
                dbconn.execute('''INSERT INTO REF_DATA (locus,base) VALUES (?,?)''',(locus,base))
            lineNumber += 1

finally:
    reffile.close()

print("done. Read {0} lines.".format(lineNumber))

# Now read the snp data into the database.
snpTable = dict()
snpLociList = []
lineNumber=1
print("Reading SNP file..." ,)
sys.stdout.flush()
try:
    with dbconn:
        snpRegex = re.compile("^(?P<strainid>.+?)\t(?P<locus>\d+)\t(?P<base>.)$")
        for line in snpfile:
            snpMatch = snpRegex.match(line)
            if snpMatch == None:
                print("*** WARNING: Unrecognized SNP format at line {0}: {1}".format(lineNumber,line))
            else:
                strainid = snpMatch.group('strainid')
                strainLocus = snpMatch.group('locus')
                strainBase = snpMatch.group('base')
                dbconn.execute('''INSERT INTO SNP_DATA (strainid,locus,base) VALUES (?,?,?)''',(strainid,strainLocus,strainBase))

            lineNumber += 1
finally:
    snpfile.close()

print("done. Read {0} lines.".format(lineNumber))

print("Processing..." ,)
sys.stdout.flush()

# Generate the PhenoLink data by loading the data table.  Locus name is <ref_base>_<locus>_<snp_base>.
# The existance of a strainid,locus_name entry in the data table indicates that the strain has a base
# at that locus.
with dbconn:
    dbconn.execute('''INSERT INTO PHENOLINK_DATA (strainid,locus_name) SELECT snp.strainid, ref.base || '_' || snp.locus || '_' || snp.base as locus_name from SNP_DATA snp inner join REF_DATA ref on snp.locus = ref.locus''')

print("Done.")


print("Writing PhenoLink export data..." ,)
sys.stdout.flush()

try:
    # First write out the column headings.  The first is the Strain ID, then all the unique SNP loci names.
    outfile.write("StrainID")

    dbcursor = dbconn.cursor()
    dbcursor.execute('select distinct locus_name from phenolink_data order by locus_name asc')
    for lociNameResult in dbcursor.fetchall():
        outfile.write("\t{0}".format(lociNameResult[0]))

    outfile.write("\n")

    straincursor = dbconn.cursor()
    locuscursor = dbconn.cursor()
    # For each strain, fill out its row values for each snp loci name.
    straincursor.execute('select distinct strainid from snp_data order by strainid asc')
    for strainIDResult in straincursor.fetchall():
        strainID = strainIDResult[0]
        outfile.write(strainID)

        dbcursor.execute('select distinct locus_name from phenolink_data order by locus_name asc')
        for locusNameResult in dbcursor.fetchall():
            locusName = locusNameResult[0]
            locuscursor.execute('SELECT 1 from phenolink_data where strainid = ? and locus_name = ? LIMIT 1',(strainID,locusName))
            locusCheckResult = locuscursor.fetchone()
            if locusCheckResult == None:
                outfile.write("\t0")
            else:
                outfile.write("\t1")

        outfile.write("\n")

finally:
    outfile.close()

print("done.\n")

print("Final output file is {0}\n".format(args.out[0]))
