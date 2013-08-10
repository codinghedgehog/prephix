# This is a helper module for Prephix. It uses a pseudo-factory pattern to
# handle the different SNP input file formats that Prephix expects to support.
#

import re
import os
import sys

#
# Public Methods
#
def getSNPFileReader(fileName,filterQuality=True):
'''
This is basically an factory method that takes
in an SNP input file name for Prephix, determines
the file format, and returns an appropriate 
reader/iterator object for that file.

SNP File Reader/Iterator objects are expected
to support iteration behavior (i.e. "for...in" statements)
that returns an SNPDataLine object with each call. 

Parameters:
    fileName = The SNP input file name.  Assumes data for one strain per file.

    filterQuality = (OPTIONAL) Used by VCF file reader for filtering on quality (default is to return only PASS quality).
        Setting filterQuality to False will have the reader return all lines without regard to quality.


'''
    fileFormat = "unknown"

    # Define regexes for recognizing the file type based on its contents.
    k28Re = re.compile("^#(.+?)\/")
    nucmerRe = re.compile("^NUCMER$")
    vcfRe = re.compile("^##fileformat=VCF")

    with open(fileName,"r") as inputFile:
        for line in inputfile:
            if k28Re.match(line):
                fileFormat = "k28"
                break
            elif nucmerRe.match(line):
                fileFormat = "nucmer"
                break
            elif vcfRe.match(line):
                fileFormat = "vcf"

    if fileFormat == "k28":
        return K28FileReader(fileName)
    elif fileFormat == "nucmer":
        return NucmerFileReader(fileName)
    elif fileFormat = "vcf":
        return VCFFileReader(fileName,filterQuality)
    else:
        raise NotImplementedError

#
# Custom Exceptions
#

class SNPFileReadError(Exception):
    '''
    This is a custom exception thrown when an error has occurred
    while reading an SNP data file.

    It requires a message, with optional line number and the
    line text itself where the error occurred.
    '''
    def __init__(self,message,lineNumber=None,line=None):
        super.__init__(self,message)
        self.lineNumber = lineNumber
        self.lineText = line


#
# Classes
#
class SNPFileReader:
    '''
    This is an abstract base class that defines the interface
    that SNPFileReader objects are expected to support.

    Subclasses should be iterable (i.e. override the __iter__ method).
    '''

    # Public Attributes.
    strainID = None # The strain name associatd with the SNP data
    fileName = None # The data source (typically filename)
    fileFormat = None # The file format

    def __init__(self,fileName):
        '''
        Default constructor.
        '''
        self.fileName = fileName
    
    def __iter__():
        raise NotImplementedError

class SNPDataLine:
    '''
    This is primary a data structure used as the return value from
    SNPFileReader's GetNextData method.
    '''
    def __init__(self,rawLine,locus,snpBase,refBase=None,isIndel=False,isInsert=False,isDelete=False):
        self.rawLine = rawLine       # Actual text line of file.
        self.locus = locus           # Locus number
        self.snpBase = snpBase       # The SNP (sample) base at the locus number for the SNP.
        self.refBase = refBase       # The reference base at the locus number for the SNP.
        self.isIndel = isIndel       # Is this line an indel?
        self.isInsert = isInsert     # If this is an indel, is it an insertion?
        self.isDelete = isDelete     # If this is an indel, is it a deletion?
    

#
# Private Classes
#

class K28FileReader(SNPFileReader):
    '''
    This is a concrete implementation of a file reader
    for k28.out (VAAL) data.
    '''

    k28lineRe = re.compile("^[0-9]+\s+(?P<locus>[0-9]+)\s+left=[ATCG]*\s+sample=(?P<sample_base>[ATCG]*)\s+ref=(?P<ref_base>[ATCG]*)\s+right=[ATCG]*$")
    strainRe = re.compile ("^#(?P<strainid>.+?)\/")

    def __init__(self,fileName):
        super(K28FileReader,self).__init__(self,fileName)
        self.fileFormat = "k28"
        self.lineNumber = 0

        # Find the strainID
        fh = open(fileName,"r")

        # Get strain ID from header comments: #<strain_id>/<reference_genome_filename>
        for line in fh:
            strainMatch = strainRe.match(line)
            if strainMatch:
                self.strainID = str(strainMatch.group('strainid'))
                break


        fh.close()

    def __iter__(self):
        '''
        Concrete implementation of iterator protocol (as a generator) to get the next line of data.
        '''
        # Open the file and skip to the first line of data.
        fh = open(fileName,"r")
        for line in fh:
            self.lineNumber += 1
            # Ignore other comments in the file (lines starting with #).  This includes the header comments.
            # Also skip the > line (don't care about genbank_id_from_ref_genome_file).
            if re.match("^(#)|(>)",line):
                continue

            # At this point, should be at a data line.

            # Assuming k28 input body data to be in the format:
            # 0 <snp_locus> <left flank seq> <sample> <ref> <right flank seq>
            # Only care about the locus, sample, and ref columns.
            #
            # Regex note:
            #
            # Some times sample= and ref= may have no value, so match for [ATCG] and check for length 1.
            # If it is not length 1, then it is either blank or have more than one base, so skip as indel.
            lineMatch = k28lineRe.match(line)
            if lineMatch:
                # VAAL k28.out file loci is offset by +1
                realLocus = int(lineMatch.group(locus)) + 1
                snpBase = lineMatch.group('sample_base')
                refBase = lineMatch.group('ref_base')

                isIndel = False
                isDelete = False
                isInsert = False

                # Check for indels.
                if len(snpBase) != 1:
                    if len(snpBase) == 0:
                        # Deletion found.
                        isIndel = True
                        isDelete = True
                    else:
                        # Insertion found.
                        isIndel = True
                        isInsert = True
                elif len(refBase) != 1:
                    if len(refBase) == 0:
                        # Insertion found.
                        isIndel = True
                        isInsert = True
                    else:
                        # Deletion found.
                        isIndel = True
                        isInsert = True


                yield SNPDataLine(line,realLocus,snpBase,refBase,isIndel,isInsert,isDelete)
            else:
                raise SNPFileReadError("Unrecognized line at {}: {}".format(self.lineNumber,line),self.lineNumber,line)

class NucmerFileReader(SNPFileReader):
    '''
    This is a concrete implementation of a file reader
    for nucmer data.
    '''

    strainRe = re.compile("\s.+\/(?P<strainid>[^\/]+)$")
    nucmerlineRe = re.compile("^(?P<locus>[0-9]+)\t(?P<ref_base>[ATCG]*)\t(?P<sample_base>[ATCG]*)\t[0-9]+")

    def __init__(self,fileName):
        super(NucmerFileReader,self).__init__(self,fileName)
        self.fileFormat = "nucmer"
        self.lineNumber = 0

        # Find the strainID
        fh = open(fileName,"r")

        # Get strain ID from header line: /path/to/reference/file /path/to/query/file
        # Assuming Strain ID is the query file name
        for line in fh:
            strainMatch = strainRe.match(line)
            if strainMatch:
                self.strainID = str(strainMatch.group('strainid'))
                break

        fh.close()

    def __iter__(self):
        '''
        Concrete implementation of iterator protocol (as a generator) to get the next line of data.
        '''
        # Open the file and skip to the first line of data.
        fh = open(fileName,"r")
        foundHeader = False
        for line in fh:
            self.lineNumber += 1

            # Keep skipping lines until we reach the data portion.  This should occur after the data header line:
            # [P1]  [SUB] [SUB] [P2]  [BUFF]  [DIST]  [LEN R] [LEN Q] [FRM] [TAGS]
            # So look for [P1]
            if not foundHeader:
                if re.match("^\[P1\]",line):
                    # Found header.  Advance to next line (a valid data line) to begin processing.
                    foundHeader = True
                continue

            # At this point, should be at a data line.

            # Assuming NUCMER input body data to be in the format:
            # [P1]  [SUB] [SUB] [P2]  [BUFF]  [DIST]  [LEN R] [LEN Q] [FRM] [TAGS]
            # Only care about the locus, ref, and sample base columns.  So, P1 and the first two SUBs, respectively.
            #
            # Assuming fields are tab-delimited.
            #
            # Regex note:
            #
            # Some times one or another SUB may have no value, so match for [ATCG] and check for length 1.
            # If it is not length 1, then it is either blank or have more than one base, so skip as indel.

            lineMatch = nucmerlineRe.match(line)
            if lineMatch:
                realLocus = int(lineMatch.group(locus))
                snpBase = lineMatch.group('sample_base')
                refBase = lineMatch.group('ref_base')

                isIndel = False
                isDelete = False
                isInsert = False

                # Check for indels.
                if len(snpBase) != 1:
                    if len(snpBase) == 0:
                        # Deletion found.
                        isIndel = True
                        isDelete = True
                    else:
                        # Insertion found.
                        isIndel = True
                        isInsert = True
                elif len(refBase) != 1:
                    if len(refBase) == 0:
                        # Insertion found.
                        isIndel = True
                        isInsert = True
                    else:
                        # Deletion found.
                        isIndel = True
                        isInsert = True


                yield SNPDataLine(line,realLocus,snpBase,refBase,isIndel,isInsert,isDelete)
            else:
                raise SNPFileReadError("Unrecognized line at {}: {}".format(self.lineNumber,line),self.lineNumber,line)

class VCFFileReader(SNPFileReader):
    '''
    This is a concrete implementation of a file reader
    for VCF data.
    '''

    vcflineRe = re.compile("^[^\t]+\t(?P<locus>[0-9]+)\t[^\t]+\t(?P<ref_base>[ATCGN,]+)\t(?P<sample_base>[ATCGN,]+)\t[^\t]+\t(?P<filter>[^\t]+)\t")

    def __init__(self,fileName,filterQuality=True):
        '''
        filterQuality parameter is a flag to indicate if it should filter data lines having only PASS quality filter values.
        '''
        super(NucmerFileReader,self).__init__(self,fileName)
        self.fileFormat = "vcf"
        self.lineNumber = 0
        self.filterQuality = filterQuality

        # Set the strainID to the filename for now.
        self.strainID = os.file.basename(fileName)


    def __iter__(self):
        '''
        Concrete implementation of iterator protocol (as a generator) to get the next line of data.
        '''
        # Open the file and skip to the first line of data.
        fh = open(fileName,"r")
        foundHeader = False
        for line in fh:
            self.lineNumber += 1

            # Keep skipping lines until we reach the data portion.  This should occur after the data header line:
            #CHROM  POS ID  REF ALT QUAL  FILTER  INFO
            # So look for #CHROM
            if not foundHeader:
                if re.match("^#CHROM\s+POS\s+ID",line):
                    # Found header!  Advance one more line to actual data line.
                    foundHeader = True
                continue


            # At this point, should be at a data line.

            # Assuming VCF input body data to be in the format:
            # CHROM  POS ID  REF ALT QUAL  FILTER  INFO
            # Only care about the pos (loci), ref, alt (sample), and filter columns.
            #
            # Assuming fields are tab-delimited.
            #

            lineMatch = vcflineRe.match(line)
            if lineMatch:
                realLocus = int(lineMatch.group(locus))
                snpBase = lineMatch.group('sample_base')
                refBase = lineMatch.group('ref_base')
                filter = lineMatch.group('filter')

                isIndel = False
                isDelete = False
                isInsert = False

                if filterQuality and filter != "PASS":
                    # Ignore low quality line.
                    continue

                # Check for indels.
                if len(refBase) != 1:
                    # Deletion found.
                    isIndel = True
                    isDelete = True
                elif len(snpBase) != 1:
                    # Insertion found.
                    isIndel = True
                    isInsert = True

                yield SNPDataLine(line,realLocus,snpBase,refBase,isIndel,isInsert,isDelete)
            else:
                raise SNPFileReadError("Unrecognized line at {}: {}".format(self.lineNumber,line),self.lineNumber,line)
