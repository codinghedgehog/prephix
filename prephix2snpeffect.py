#!/usr/bin/python
#
# Prephix 2 SNP Effect Converter
#
# This script takes the ref and snp output files from a prephix run and generates an SNP Effect compliant input file.
#
# Usage: prephix2snpeffect.py ref_file snp_file

import sys
import os
import re
import argparse
import string
import cStringIO
import math

VERSION = '1.0.0'

# MAIN #

print "\nPrephix to SNP Effect Converter v{0}\n".format(VERSION)

# Define our parameters and parse with argparse module.
argParser = argparse.ArgumentParser()

progname = os.path.abspath(sys.argv[0])

# Define arguments.
argParser.add_argument("ref_file",help="The prephix written reference output file.")
argParser.add_argument("snp_file",help="The prephix written snp file output.")
argParser.add_argument("-d","--debug",action="store_true",help="Debug mode.")

args = argParser.parse_args()

# Set debug mode.
debugMode = args.debug

# Working variables.
refFilename = args.ref_file
snpFilename = args.snp_file

if debugMode:
    print "Argparse results: "
    print args

