# This is a helper module for Prephix. It uses a pseudo-factory pattern to
# handle the different SNP input file formats that Prephix expects to support.
#
# Each SNP input file reader class is expected to be ITERABLE, and returns
# a tupe of (line,lineNumber,realLocus,snpBase,refBase,isIndel,isInsert,isDelete)
#
# Where:
#
# line = The raw line as read from the input file.
# lineNumber = The line number of the file line.
# realLocus = Adjusted locus value (some formats have an offset relative to the base reference)
# snpBase = The sample base at the given locus
# refBase = The reference base at the given locus
# isIndel = Boolean to indicate if this line is an indel.
# isInsert = Boolean to indicate if this line is an indel insertion.
# isDelete  = Boolean to indicate if this line is an indel deletion.
#

import re
import os
import sys

#
# Public Methods
#
def getSNPFileReader(fileName,filterQuality=True,multiChrom=False):
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

        multiChrom = (OPTIONAL) Used by VCF file read for indicating  multi-Chrom parsing or not.


    '''
    fileFormat = "unknown"

    # Define regexes for recognizing the file type based on its contents.
    k28Re = re.compile("^#(.+?)\/")
    nucmerRe = re.compile("^NUCMER$")
    vcfRe = re.compile("^##fileformat=VCF")

    with open(fileName,"r") as inputFile:
        for line in inputFile:
            if k28Re.match(line):
                fileFormat = "k28"
                break
            elif nucmerRe.match(line):
                fileFormat = "nucmer"
                break
            elif vcfRe.match(line):
                fileFormat = "vcf"
                break

    if fileFormat == "k28":
        return K28FileReader(fileName)
    elif fileFormat == "nucmer":
        return NucmerFileReader(fileName)
    elif fileFormat == "vcf":
        return VCFFileReader(fileName,filterQuality,multiChrom)
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
        super(SNPFileReadError,self).__init__(message)
        self.lineNumber = lineNumber
        self.lineText = line

class SNPFileUnrecognizedLineError(Exception):
    '''
    This is a custom exception thrown when an unrecognized line
    is read from an SNP data file.

    It requires a message, with optional line number and the
    line text itself where the error occurred.
    '''
    def __init__(self,message,lineNumber=None,line=None):
        super(SNPFileUnrecognizedLineError,self).__init__(message)
        self.lineNumber = lineNumber
        self.lineText = line


#
# Classes
#
class SNPFileReader(object):
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

#
# Private Classes
#

class K28FileReader(SNPFileReader):
    '''
    This is a concrete implementation of a file reader
    for k28.out (VAAL) data.
    '''


    def __init__(self,fileName):
        super(K28FileReader,self).__init__(fileName)
        self.fileFormat = "k28"
        self.lineNumber = 0

        self.k28lineRe = re.compile("^[0-9]+\s+(?P<locus>[0-9]+)\s+left=[ATCG]*\s+sample=(?P<sample_base>[ATCG]*)\s+ref=(?P<ref_base>[ATCG]*)\s+right=[ATCG]*$")

        # Find the strainID
        fh = open(self.fileName,"r")
        self._fileLines = fh.readlines()
        fh.close()

        # Get strain ID from header comments: #<strain_id>/<reference_genome_filename>
        strainRe = re.compile ("^#(?P<strainid>.+?)\/")
        for line in self._fileLines:
            self.lineNumber +=1
            strainMatch = strainRe.match(line)
            if strainMatch:
                self.strainID = str(strainMatch.group('strainid'))
                break

        # Skip to the first line of data.
        while self.lineNumber <= len(self._fileLines):
            line = self._fileLines[self.lineNumber - 1]
            # Ignore other comments in the file (lines starting with #).  This includes the header comments.
            # Also skip the > line (don't care about genbank_id_from_ref_genome_file).
            if not re.match("^(#|>)",line):
                # Exit now since this IS a data line. Want to end reading right here.
                break
            else:
                self.lineNumber += 1


    def __iter__(self):
        '''
        Concrete implementation of iterator protocol (as a generator) to get the next line of data.
        '''

        while self.lineNumber <= len(self._fileLines):

            # At this point, should be at a data line.

            # Assuming k28 input body data to be in the format:
            # 0 <snp_locus> <left flank seq> <sample> <ref> <right flank seq>
            # Only care about the locus, sample, and ref columns.
            #
            # Regex note:
            #
            # Some times sample= and ref= may have no value, so match for [ATCG] and check for length 1.
            # If it is not length 1, then it is either blank or have more than one base, so skip as indel.
            lineNumber = self.lineNumber
            line = self._fileLines[lineNumber - 1].strip(os.linesep)
            self.lineNumber += 1

            # Skip comments (for K28/VAAL files, there are some comments at the end, not just the beginning).
            if line.startswith("#"):
                continue

            # Else parse the line.
            lineMatch = self.k28lineRe.match(line)
            if lineMatch:
                # VAAL k28.out file loci is offset by +1
                realLocus = int(lineMatch.group('locus')) + 1
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
                        isDelete = True


                yield (line,lineNumber,realLocus,snpBase,refBase,isIndel,isInsert,isDelete)
            else:
                raise SNPFileUnrecognizedLineError("Unrecognized line at {}: {}".format(lineNumber,line),lineNumber,line)

class NucmerFileReader(SNPFileReader):
    '''
    This is a concrete implementation of a file reader
    for nucmer data.
    '''


    def __init__(self,fileName):
        super(NucmerFileReader,self).__init__(fileName)
        self.fileFormat = "nucmer"
        self.lineNumber = 0

        self.nucmerlineRe = re.compile("^(?P<locus>[0-9]+)\t(?P<ref_base>[ATCG]*)\t(?P<sample_base>[ATCG]*)\t[0-9]+")

        # Find the strainID
        fh = open(fileName,"r")
        self._fileLines = fh.readlines()
        fh.close()

        # Get strain ID from header line: /path/to/reference/file /path/to/query/file
        strainRe = re.compile("\s.+\/(?P<strainid>[^\/]+)$")

        # Assuming Strain ID is the query file name
        for line in self._fileLines:
            self.lineNumber += 1
            strainMatch = strainRe.search(line)
            if strainMatch:
                self.strainID = str(strainMatch.group('strainid')).strip(os.linesep)
                break

        # Open the file and skip to the first line of data.
        while self.lineNumber <= len(self._fileLines):
            line = self._fileLines[self.lineNumber - 1].strip(os.linesep)

            # Keep skipping lines until we reach the data portion.  This should occur after the data header line:
            # [P1]  [SUB] [SUB] [P2]  [BUFF]  [DIST]  [LEN R] [LEN Q] [FRM] [TAGS]
            # So look for [P1]
            if re.match("^\[P1\]",line):
                # Found header. Advance line number to the first data line (first line after this header) and stop.
                self.lineNumber += 1
                break

            self.lineNumber += 1

    def __iter__(self):
        '''
        Concrete implementation of iterator protocol (as a generator) to get the next line of data.
        '''
        while self.lineNumber <= len(self._fileLines):
            lineNumber = self.lineNumber
            line = self._fileLines[lineNumber - 1].strip(os.linesep)
            self.lineNumber += 1

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

            lineMatch = self.nucmerlineRe.search(line)
            if lineMatch:
                realLocus = int(lineMatch.group('locus'))
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
                        isDelete = True


                yield (line,lineNumber,realLocus,snpBase,refBase,isIndel,isInsert,isDelete)
            else:
                raise SNPFileUnrecognizedLineError("Unrecognized line at {}: {}".format(lineNumber,line),lineNumber,line)

class VCFFileReader(SNPFileReader):
    '''
    This is a concrete implementation of a file reader
    for VCF data.
    '''


    def __init__(self,fileName,filterQuality=True, multiChrom=False):
        '''
        filterQuality parameter is a flag to indicate if it should filter data lines having only PASS quality filter values.
        '''
        super(VCFFileReader,self).__init__(fileName)
        self.fileFormat = "vcf"
        self.filterQuality = filterQuality
        self.lineNumber = 0
        self._fileLines = []
        self.multiChrom = multiChrom

        # NOTE: If processing multi-chrom FASTA files, Locus becomes CHROM + POS, instead of just POS.
        self.vcflineRe = re.compile("^(?P<chrom>[^\t]+)\t(?P<locus>[0-9]+)\t[^\t]+\t(?P<ref_base>[ATCGN,]+)\t(?P<sample_base>[ATCGN,]+)\t[^\t]+\t(?P<filter>[^\t]+)\t")

        # Set the strainID to the filename for now.
        self.strainID = os.path.basename(fileName)

        # Open the file and skip to the first line of data.
        # Move past the header lines (not needed since it's been determined elsewhere what this file type is).
        fh = open(fileName,"r")
        self._fileLines = fh.readlines()
        fh.close()
        for line in self._fileLines:
            self.lineNumber += 1

            # Keep skipping lines until we reach the data portion.  This should occur after the data header line:
            #CHROM  POS ID  REF ALT QUAL  FILTER  INFO
            # So look for #CHROM
            if re.match("^#CHROM\s+POS\s+ID",line):
                # Found header!  Advance one more line to actual data line and stop looking.
                self.lineNumber += 1
                #print "lineNumber is " + str(self.lineNumber) + " and filesLines is "  + str(len(self._fileLines))
                break


    def __iter__(self):
        '''
        Concrete implementation of iterator protocol (as a generator) to get the next line of data.
        '''

        while self.lineNumber <= len(self._fileLines):
            lineNumber = self.lineNumber
            line = self._fileLines[lineNumber - 1].strip(os.linesep)
            self.lineNumber += 1
            #print "lineNumber is " + str(self.lineNumber) + " and filesLines is "  + str(len(self._fileLines))
            # At this point, should be at a data line.

            # Assuming VCF input body data to be in the format:
            # CHROM  POS ID  REF ALT QUAL  FILTER  INFO
            # Only care about the pos (loci), ref, alt (sample), and filter columns.
            #
            # Assuming fields are tab-delimited.
            #

            lineMatch = self.vcflineRe.match(line)
            if lineMatch:
                if self.multiChrom:
                    realLocus = str(lineMatch.group('chrom')) + '-' + str(lineMatch.group('locus')) 
                else:
                    realLocus = int(lineMatch.group('locus'))
                snpBase = lineMatch.group('sample_base')
                refBase = lineMatch.group('ref_base')
                filter = lineMatch.group('filter')

                isIndel = False
                isDelete = False
                isInsert = False

                if self.filterQuality and filter != "PASS":
                    # Ignore low quality line.
                    continue

                # Check for indels.
                if len(snpBase) != len(refBase):
                    if len(snpBase) < len(refBase):
                        # Deletion found.
                        isIndel = True
                        isDelete = True
                    else:
                        # Insertion found.
                        isIndel = True
                        isInsert = True

                yield (line,lineNumber,realLocus,snpBase,refBase,isIndel,isInsert,isDelete)
            else:
                raise SNPFileUnrecognizedLineError("Unrecognized line at {}: {}".format(lineNumber,line),lineNumber,line)
