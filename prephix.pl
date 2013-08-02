#!/usr/bin/perl
#
# Prephix (Pre-Phrecon Input fiXer) 
#
# Usage: prephix.pl k28.out_file1 [k28.out_file2 ...] -batchid <id> [-debug] [-quiet] [-exclude <loci_file>] [-ignore_quality]
#
# Takes in a list of input files and generates one SNP loci output file (containing multiple strains), one reference base file,
# and a file containing the excluded indel entries.  
#
# The indel file is formated as "STRAIN ID" [TAB] k28|nuc|vcf [TAB] ...excluded VAAL K28, VCF, or NUCMER line...
# So prephix basically dumps the excluded lines as-is into the indel file, but prefixes each line
# with the strain id and then either k28 for VAAL formatted line, nuc for NUCMER formatted line, or vcf for VCF formatted line.
#
# The SNP loci and reference base files are suitable as input to the phrecon (Phylo Reconstructor) utility (i.e. tab 
# delimited columnar format: STRAIN_ID [TAB] LOCI [TAB] BASE).
#
# The SNP loci file and indel files are suitable input for the SNP Swapper utility.
#
# Recognized input file formats: k28.out (VAAL), NUCMER, VCF
#
# You may mix input file types within a single run/batch.
#
# 2/01/2012 - Andrew Pann - Initial development.
# 2/02/2012 - Andrew Pann - v2.1.  Minor textual changes, and more debugging reporting.
# 2/03/2012 - Andrew Pann - v2.1.2 Minor text changes.  Reports proper log file name if it can't open it for writing.  Streamlined regex expression (removed nested captures for sample and ref base).  Fixed loci offset from VAAL input.
# 5/04/2012 - Andrew Pann - v2.2.0 Added support for NUCMER files.
# 5/04/2012 - Andrew Pann - v2.2.1 Added support for loci exclusions.
# 5/18/2012 - Andrew Pann - v2.3.0 Updated to generate just one indel file, with mixed VAAL and NUCMER indel lines with strain id in first column.
# 7/07/2012 - Andrew Pann - v2.3.1 Added support for VCF files.
# 11/14/2012 - Andrew Pann - v2.4.0 Added stats reporting.
# 11/14/2012 - Andrew Pann - v2.4.1 Added total indel count and loci exclusions to stats.  Fixed stats logic.
# 7/19/2013 - Andrew Pann v2.5.0 Added reporting breakdown for loci exclusions per strain.
# 7/20/2013 - Andrew Pann v2.5.1 Added PhenoLink export option. This uses the pre2phe.py script, called externally.
# 7/31/2013 - Andrew Pann v2.5.2 Fixed insertion/deletion reporting for k28 files.
# 
#

use strict;

my $VERSION="2.5.2";

print "\nPrephix (Pre-Phrecon Input fiXer) v$VERSION\n\n";

if ($#ARGV < 0){
  print "Usage: $0 input_file1 [input_file2 ...] -batchid <id> [-debug] [-quiet] [-exclude <loci_file>] [-ignore_quality] [-tablog]\n";
	print "-batch <id> is used to create the output filenames for the combined SNP loci and ref files.\n";
	print "-exclude <loci_file> indicates to prephix to exclude any bases within the loci ranges in the\n";
	print "<loci_file>.  The file should contain lines of the format 'label,start_loci,end_loci' \n";
	print "-debug and -quiet flags do what you expect them to do (i.e. verbose output and surpressed output).\n";
	print "-ignore_quality (currently only applies to VCF files) will process all lines, not just lines with QUALITY=PASS.\n";
	print "-tablog instructs prephix to print out the summary in tabular format.\n";
	print "-export_phenolink instructs prephix to also write out a PhenoLink compatible output file.\n";
  exit 1;
}

my $batchid="";
my $inFilename;
my $reffile;
my $refFilename;
my $outfile;
my $outFilename;
my $excludeFile;
my $excludeFilename;
my $exportPhenoLink="N";
my $logfile;
my $i=0;
my $fileCount=0;
my $loci;
my $charCount=0;
my $warnings=0;
my $debug="N";
my $quiet="N";
my $arg_num=0;
my @inputFileList;
my %refBaseTable;
my %exclusionTable;
my $realLoci;
my $input_type="k28";
my $indelfile;
my $indelFilename;
my $ignore_quality="N";
my $tablog = "N";
my $phenoStrain;


my %reportHash;  # Format is {strainId => [# SNPs, # insertions, # deletions, # exclusions]} (key is strainID, value is an array).
my %exclusionReportHash; # Format is {strainID => {Exclusion strainID => count}}
my %phenoLinkStrainReportHash; # Format is {strainID => {Loci => Base}}

# Process parameters
$arg_num=0;
while ($arg_num <= $#ARGV){
  # Using if-else logic because given/when isn't compatible with< 5.10 and Switch module was removed from >5.14.
  if ($ARGV[$arg_num] eq "-debug"){
      $debug="Y";
      print "Producing debug output.\n";
  }
  elsif ($ARGV[$arg_num] eq "-quiet") { 
      print "Producing quiet (no stdout) output.  Log file is still generated.\n";
      $quiet="Y";
  }
	elsif ($ARGV[$arg_num] eq "-exclude" ){
			$arg_num++;
			$excludeFilename="$ARGV[$arg_num]";
			print "Exclusion file is $excludeFilename\n";
  }
	elsif ($ARGV[$arg_num] eq "-batchid"){
			$arg_num++;
			$batchid="$ARGV[$arg_num]";
			print "Batch id is $batchid\n";
  }
	elsif ($ARGV[$arg_num] eq "-ignore_quality"){
			$ignore_quality="Y";
			print "Will process all lines, ignoring quality value (only applicable for vcf files).\n";
  }
	elsif ($ARGV[$arg_num] eq "-tablog"){
		  print "Will format summary in tabular format.\n";
		  $tablog = "Y";
  }
	elsif ($ARGV[$arg_num] eq "-export_phenolink"){
		  print "Will export a PhenoLink file.\n";
		  $exportPhenoLink = "Y";
  }
  elsif ( -r $ARGV[$arg_num] ){
			push(@inputFileList, $ARGV[$arg_num]);
			print "Found $ARGV[$arg_num] file to process.\n";
			$fileCount++;
  }
  else{
      print "*** ERROR: Unknown command line argument or file $ARGV[$arg_num].\n";
      exit 1;
  }
	$arg_num++;
}

if ($batchid eq ""){
  print "Missing required parameter: -batchid <id>\n\n";
  print "Usage: $0 -batchid <id> k28.out_file1 [k28.out_file2 ...3 ...] [-debug] [-quiet]\n";
	print "-batch <id> is used to generate the output filenames for the combined SNP loci and ref files.\n";
  exit 1;
}

if ($#inputFileList < 0 ){
  print "*** ERROR: No input files found.  Quitting...\n";
	exit 1;
}

open($logfile,">","$batchid.log") or die "unable to open log file $batchid.log for writing! $!\n";

$outFilename = "$batchid.snp";
open($outfile,">","$outFilename") or die "unable to open SNP loci output file $outFilename for writing! $!\n";

$refFilename = "$batchid.ref";
open($reffile,">","$refFilename") or die "Unable to open output reference base file $refFilename for writing!  $!\n";

$indelFilename = "$batchid.indel";
open($indelfile,">","$indelFilename") or die "Unable to open output indel file $indelFilename for writing!  $!\n";

if ($excludeFilename ne ""){
  print_all("Reading exclusion list from $excludeFilename.\n");

  open($excludeFile,"<",$excludeFilename) or die "Unable to read exclusion file $excludeFilename! $!\n";

  my $excludeCount=0;
  # Expect exclusion file to contain one loci range per line, formatted: label,start_loci,end_loci (inclusive)
	while (<$excludeFile>){
		chomp;
		# Store exclusions in hash with label as key and start,end loci as array values.
		if (/^([^,]+),([0-9]+),([0-9]+)$/){
			if (exists($exclusionTable{"$1"})){
				print_all("Duplicate exclusion label encountered: $_.  Quitting!\n");
				exit 1;
			}
			else{
				$exclusionTable{"$1"} = [$2,$3];
				$excludeCount++;
			}	
		}
		else{
			print_all("Badly formatted line: $_.  Quitting!\n");
			exit 1;
		}
	}
	close($excludeFile);
	print_all("$excludeCount exclusions read.\n");
}

print_all("\n***\n");
print_all("*** REMINDER: This program assumes that all input files refer to the same reference.\n");
print_all("***           Ref file will be generated from consolidated ref loci information of ALL input files.\n");
print_all("***\n");

print_all("\nFound $fileCount input files...\n\n");


$fileCount=0;

# Read in each input file and generate a SNP loci file and base reference file
while ($inFilename = pop(@inputFileList)){

  $input_type = get_input_type($inFilename);

  # Process each file based on its type.  Note it fills and uses the same global hash variables (to consolidate across allinput files).
  if ($input_type eq "k28"){
	  do_k28_file($inFilename);
  }
	elsif ($input_type eq "nucmer"){
	  do_nucmer_file($inFilename);
  }
	elsif ($input_type eq "vcf"){
	  do_vcf_file($inFilename);
  }
  else{
	  print_all("*** ERROR: Unknown or unset input file type: $input_type\n");
		exit 1;
  }
	$fileCount++;
}

print_all("$fileCount files were processed.\n");

print_all("\nMerging and generating reference file from input file data...");

# Write out the reference file from the table of merged ref loci bases from the input file.
# Output format is Loci [TAB] Base
$charCount = 0;
foreach $loci (sort {$a <=> $b} keys %refBaseTable){

  print $reffile "$loci\t$refBaseTable{$loci}[0]\n";

}

print_all("Done.\n");

if ($exportPhenoLink eq "Y"){
	print_all("Exporting PhenoLink file...");

	if ( ! -x "./pre2phe.py" ){
		print_all("*** ERROR: Cannot find external script pre2phe.py to generate PhenoLink output!\n");
	}
  else{	
		system("./pre2phe.py --ref $refFilename --snp $outFilename --out $batchid.phenolink.txt");
		print_all("PhenoLink export file is $batchid.phenolink.txt\n");
	}
}

close($outfile);
close($reffile);
close($indelfile);
print_all("\nMerged SNP loci file is $outFilename\n");
print_all("Merged reference base file from this run is $refFilename\n");
print_all("Merged indel file from this run is $indelFilename\n");


# Write out stats report.
print_all("\n=== Final Report ===\n\n");

if ($tablog eq "N"){

	foreach my $strain (keys %reportHash){
		print_all("Strain: $strain\n");
		print_all("SNPs: $reportHash{$strain}[0]\n");
		print_all("Insertions: $reportHash{$strain}[1]\n");
		print_all("Deletions: $reportHash{$strain}[2]\n");
		my $totalIndels =  $reportHash{$strain}[1] + $reportHash{$strain}[2];
		print_all("Total indels: $totalIndels\n");
		print_all("Loci excluded: $reportHash{$strain}[3]\n");

		foreach my $excludeStrainName (keys %{$exclusionReportHash{"$strain"}}){
			print_all("* Loci excluded from $excludeStrainName: " . $exclusionReportHash{"$strain"}{"$excludeStrainName"} . "\n");
		}
		print_all("\n");
	}


}
else{
	# Tabular format requested by -tablog flag.

	print_all("Strain ID\tSNPs\tInserts\tDeletes\tTotal Indels\tLoci Excluded\n");
	foreach my $strain (keys %reportHash){
		print_all("$strain\t$reportHash{$strain}[0]\t$reportHash{$strain}[1]\t$reportHash{$strain}[2]\t");
		my $totalIndels =  $reportHash{$strain}[1] + $reportHash{$strain}[2];
		print_all("$totalIndels\t$reportHash{$strain}[3]\n");
	}

	print "\n";

	foreach my $strain (keys %reportHash){
		print_all("Loci exclusion summary for $strain:\n");
		foreach my $excludeStrainName (keys %{$exclusionReportHash{"$strain"}}){
			print_all("* Loci excluded from $excludeStrainName: " . $exclusionReportHash{"$strain"}{"$excludeStrainName"} . "\n");
		}
	}
}
print_all("\nDone.\n");

sub get_exclusion_loci
{
  my $loci=$_[0];
  my $label;

	foreach $label (keys %exclusionTable){
	  if (($loci >= $exclusionTable{$label}[0]) && ($loci <= $exclusionTable{$label}[1])){
		  return $label;
    }
	}
	return "NO_EXCLUSION";
}

sub get_input_type
{
  my $fileName = $_[0];
	my $infile;
	my $my_file_type="unknown";

  open($infile,"<",$fileName) or die "Unable to open input file $inFilename for reading!  $!\n";
  while (<$infile>){
    $i++;
    chomp;

    # See if header matches k28 (VAAL) file format. Get strain ID from header comments: #<strain_id>/<reference_genome_filename>
    if (/^#(.+?)\//){
		  $my_file_type = "k28";
		  last;
    }
		# See if it is a NUCMER file.
		elsif (/^NUCMER$/){
		  $my_file_type = "nucmer";
			last;
    }
		# See if it is a VCF file.
		elsif (/^##fileformat=VCF/){
		  $my_file_type = "vcf";
			last;
	  }
	}
	close($infile);
	return $my_file_type;
}

sub do_vcf_file
{
  my $inFilename = $_[0];
	my $infile;
	my $currentStrain="";
	my $exclude_name="";
	my $quality="";

  open($infile,"<",$inFilename) or die "Unable to open input file $inFilename for reading!  $!\n";

  print_all("Processing VCF file $inFilename...");

	# Set strainID to file name for now.  Strip the path, of course.
	if ($inFilename =~ /\/([^\/]+)$/){
    $currentStrain=$1;
	}
	else{
    $currentStrain=$inFilename;
	}
	print_debug("Strain ID is $currentStrain");

	# Initialize report entry for this strain.
	$reportHash{$currentStrain} = [0,0,0,0];

  $i=0;
  my $dataStart="N";
  while (<$infile>){

    $i++;
    chomp;


    # Keep skipping lines until we reach the data portion.  This should occur after the data header line:
		#CHROM  POS ID  REF ALT QUAL  FILTER  INFO
		# So look for #CHROM
		if (/^#CHROM/){
		  $dataStart="Y";
      next;
    }	
    elsif ($dataStart eq "N"){
			print_debug("Ignoring line $1: $_ (in $inFilename)\n");
		  next;
    }

		if (/^#/){
		  # Ignore other comments in the file (lines starting with #).
			print_debug("Ignoring comment line $1: $_ (in $inFilename)\n");
		  next;
    }

    # Assuming VCF input body data to be in the format:
		# CHROM  POS ID  REF ALT QUAL  FILTER  INFO
		# Only care about the pos (loci), ref, and alt (sample) base columns.  
		#
		# Assuming fields are tab-delimited.
    #
    if (/^[^\t]+\t([0-9]+)\t[^\t]+\t([ATCGN,]+)\t([ATCGN,]+)\t[^\t]+\t([^\t]+)\t/){
			$realLoci=$1;
      $quality=$4;
			if (($quality ne "PASS") && ($ignore_quality eq "N")){
        print_debug("Skipping low-quality line: $_ ($inFilename)\n");
				next;
			}
      if (length($2) != 1 ){
        print_debug("Found indel line (deletion): $_ ($inFilename)\n");

        print_debug("Skipping indel line (deletion): $_ ($inFilename)\n");
        print $indelfile "$currentStrain\tvcf\t$_\n";

				# Increment deletion count in report
				$reportHash{"$currentStrain"}[2]++;
				next;
      }

      if (length($3) != 1){
        print_debug("Found indel line (insertion): $_ ($inFilename)\n");

        print_debug("Skipping indel line (insertion): $_ ($inFilename)\n");
        print $indelfile "$currentStrain\tvcf\t$_\n";

				# Increment insertion count in report
				$reportHash{"$currentStrain"}[1]++;
				next;
      }

			$exclude_name = get_exclusion_loci($realLoci);
      if ($exclude_name eq "NO_EXCLUSION"){
        print $outfile "$currentStrain\t$realLoci\t$3\n"; # SNP loci file format is (StrainId [TAB] Loci [TAB] Base)

				if ( (exists($refBaseTable{$realLoci})) && ($refBaseTable{$realLoci}[0] ne "$2") ){
					print_all("*** ERROR: Reference base mismatch at same loci! Input file $inFilename line $i has ref=$2, but ref base at this loci was already recorded as $refBaseTable{$realLoci}[0] while processing file $refBaseTable{$realLoci}[1] line $refBaseTable{$realLoci}[2]!\n"); 
		      print_all("\n*** Are you sure all input files are using the same reference?\n");
					print "Failed.\n";
					exit 1;
				}
				else{
					# Storing ref base in hash table; use loci as key, and value is array of base, input file name, and file line number.
					# The file name and line number help us with reporting mismatched loci base collisions.
					$refBaseTable{$realLoci}=["$2","$inFilename","$i"];

					# Do report - increment SNP count.
					$reportHash{"$currentStrain"}[0]++;

				}
      }
			else{
        print_debug("Excluded loci $realLoci.\n");

				# Do report - increment exclusion count.
				$reportHash{"$currentStrain"}[3]++;

				# Do tally in the exclusion report.
				if ( exists($exclusionReportHash{"$currentStrain"}) ){
					$exclusionReportHash{"$currentStrain"}{"$exclude_name"}++;
				}
				else{
					$exclusionReportHash{"$currentStrain"}{"$exclude_name"} = 1;
				}	
			}
    }
    else{
      print_all("*** ERROR: VCF file format not recognized. Got \"$_\" at line $i.\n");
			print "Failed.\n";
      exit 1;
    }
  }
	$currentStrain="";
  close($infile);

  print_all("Done.\n\n");


}

sub do_nucmer_file
{
  my $inFilename = $_[0];
	my $infile;
	my $currentStrain="";
	my $exclude_name="";

  open($infile,"<",$inFilename) or die "Unable to open input file $inFilename for reading!  $!\n";

  print_all("Processing NUCMER file $inFilename...");

  $i=0;
  my $dataStart="N";
  while (<$infile>){

    $i++;
    chomp;

    # Get strain ID from header line: /path/to/referenc/file /path/to/query/file
		# Assuming Strain ID is the query file
    if (/\s.+\/([^\/]+)$/){
		  if ($currentStrain eq "" ){
				$currentStrain=$1;
				print_debug("$inFilename appears to be NUCMER format.\n");
				print_debug("Found strain id of $currentStrain in file $inFilename\n");

				# Initialize report entry for this strain.
				$reportHash{"$currentStrain"} = [0,0,0,0];
				next;
			}
			else{
				print_all("***ERROR: Encountered a second strain id in same file?  On line $i: $_\n");
				print "Failed.\n";
				exit 1;
			}
    }
		elsif (/^#/){
		  # Ignore other comments in the file (lines starting with #).
			print_debug("Ignoring comment line $1: $_ (in $inFilename)\n");
		  next;
    }

    # Keep skipping lines until we reach the data portion.  This should occur after the data header line:
		# [P1]  [SUB] [SUB] [P2]  [BUFF]  [DIST]  [LEN R] [LEN Q] [FRM] [TAGS]
		# So look for [P1]
		if (/^\[P1\]/){
		  $dataStart="Y";
      next;
    }	
    elsif ($dataStart eq "N"){
			print_debug("Ignoring line $1: $_ (in $inFilename)\n");
		  next;
    }

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
    if (/^([0-9]+)\t([ATCG]*)\t([ATCG]*)\t[0-9]+/){
			$realLoci=$1;
      if (length($2) != 1){
				if (length($2) == 0){
					print_debug("Found indel line (insertion): $_ ($inFilename)\n");

					# Increment insertion count in report
					$reportHash{"$currentStrain"}[1]++;
				}
        else{
					print_debug("Found indel line (deletion): $_ ($inFilename)\n");

					# Increment deletion count in report
					$reportHash{"$currentStrain"}[2]++;
				}
        print_debug("Skipping indel line: $_ ($inFilename)\n");
        print $indelfile "$currentStrain\tnuc\t$_\n";
				next;
      }

      if (length($3) != 1){
				if (length($3) == 0){
					print_debug("Found indel line (deletion): $_ ($inFilename)\n");

					# Increment deletion count in report
					$reportHash{"$currentStrain"}[2]++;
				}
        else{
					print_debug("Found indel line (insertion): $_ ($inFilename)\n");

					# Increment insertion count in report
					$reportHash{"$currentStrain"}[1]++;
				}
        print_debug("Skipping indel line: $_ ($inFilename)\n");
        print $indelfile "$currentStrain\tnuc\t$_\n";
				next;
      }

			$exclude_name = get_exclusion_loci($realLoci);
      if($exclude_name eq "NO_EXCLUSION"){
        print $outfile "$currentStrain\t$realLoci\t$3\n"; # SNP loci file format is (StrainId [TAB] Loci [TAB] Base)

				if ( (exists($refBaseTable{$realLoci})) && ($refBaseTable{$realLoci}[0] ne "$2") ){
					print_all("*** ERROR: Reference base mismatch at same loci! Input file $inFilename line $i has ref=$2, but ref base at this loci was already recorded as $refBaseTable{$realLoci}[0] while processing file $refBaseTable{$realLoci}[1] line $refBaseTable{$realLoci}[2]!\n"); 
		      print_all("\n*** Are you sure all input files are using the same reference?\n");
					print "Failed.\n";
					exit 1;
				}
				else{
					# Storing ref base in hash table; use loci as key, and value is array of base, input file name, and file line number.
					# The file name and line number help us with reporting mismatched loci base collisions.
					$refBaseTable{$realLoci}=["$2","$inFilename","$i"];

					# Do report - increment SNP count.
					$reportHash{"$currentStrain"}[0]++;

				}
      }
			else{
        print_debug("Excluded loci $realLoci.\n");

				# Do report - increment exclusion count.
				$reportHash{"$currentStrain"}[3]++;

				# Do tally in the exclusion report.
				if ( exists($exclusionReportHash{"$currentStrain"}) ){
					$exclusionReportHash{"$currentStrain"}{"$exclude_name"}++;
				}
				else{
					$exclusionReportHash{"$currentStrain"}{"$exclude_name"} = 1;
				}	
			}
    }
    else{
      print_all("*** ERROR: NUCMER file format not recognized. Got \"$_\" at line $i.\n");
			print "Failed.\n";
      exit 1;
    }
  }
	$currentStrain="";
  close($infile);

  print_all("Done.\n\n");

}

sub do_k28_file
{
  my $inFilename = $_[0];
	my $infile;
	my $currentStrain;
	my $exclude_name="";

  open($infile,"<",$inFilename) or die "Unable to open input file $inFilename for reading!  $!\n";

  print_all("Processing k28.out (VAAL) file $inFilename...");

  $i=0;
  while (<$infile>){

    $i++;
    chomp;

    # Get strain ID from header comments: #<strain_id>/<reference_genome_filename>
    if (/^#(.+?)\//){
		  if ($currentStrain eq "" ){
				$currentStrain=$1;
				print_debug("$inFilename appears to be k28.out/VAAL format.\n");
				print_debug("Found strain id of $currentStrain in file $inFilename\n");

				# Initialize report entry for this strain.
				$reportHash{"$currentStrain"} = [0,0,0,0];
				next;
			}
			else{
				print_all("***ERROR: Encountered a second strain id in same file?  On line $i: $_\n");
				print "Failed.\n";
				exit 1;
			}
    }
		elsif (/^#/){
		  # Ignore other comments in the file (lines starting with #).
			print_debug("Ignoring comment line $1: $_ (in $inFilename)\n");
		  next;
    }

    # Skip the > line (don't care about genbank_id_from_ref_genome_file).
    if (/^>/){
			print_debug("Ignoring line $1: $_ (in $inFilename)\n");
		  next;
    }

    # Assuming k28 input body data to be in the format:
		# 0 <snp_locus> <left flank seq> <sample> <ref> <right flank seq>
		# Only care about the locus, sample, and ref columns.
    #
		# Regex note:
    #
    # Some times sample= and ref= may have no value, so match for [ATCG] and check for length 1.
    # If it is not length 1, then it is either blank or have more than one base, so skip as indel.
    if (/^[0-9]+\s+([0-9]+)\s+left=[ATCG]*\s+sample=([ATCG]*)\s+ref=([ATCG]*)\s+right=[ATCG]*$/){
			$realLoci=$1;
			$realLoci++; # VAAL k28.out file loci is offset by +1
      if (length($2) != 1){
				if (length($2) == 0){
					print_debug("Found indel line (deletion): $_ ($inFilename)\n");

					# Increment deletion count in report
					$reportHash{"$currentStrain"}[2]++;
				}
        else{
					print_debug("Found indel line (insertion): $_ ($inFilename)\n");

					# Increment insertion count in report
					$reportHash{"$currentStrain"}[1]++;
				}
        print_debug("Skipping indel line: $_ ($inFilename)\n");
        print $indelfile "$currentStrain\tk28\t$_\n";
				next;
      }

      if (length($3) != 1){
				if (length($3) == 0){
					print_debug("Found indel line (insertion): $_ ($inFilename)\n");

					# Increment deletion count in report
					$reportHash{"$currentStrain"}[1]++;
				}
        else{
					print_debug("Found indel line (deletion): $_ ($inFilename)\n");

					# Increment insertion count in report
					$reportHash{"$currentStrain"}[2]++;
				}
        print_debug("Skipping indel line: $_ ($inFilename)\n");
        print $indelfile "$currentStrain\tk28\t$_\n";
				next;
      }

			$exclude_name = get_exclusion_loci($realLoci);
      if ($exclude_name eq "NO_EXCLUSION"){
        print $outfile "$currentStrain\t$realLoci\t$2\n"; # SNP loci file format is (StrainId [TAB] Loci [TAB] Base)

				if ( (exists($refBaseTable{$realLoci})) && ($refBaseTable{$realLoci}[0] ne "$3") ){
					print_all("*** ERROR: Reference base mismatch at same loci! Input file $inFilename line $i has ref=$3, but ref base at this loci was already recorded as $refBaseTable{$realLoci}[0] while processing file $refBaseTable{$realLoci}[1] line $refBaseTable{$realLoci}[2]!\n"); 
		      print_all("\n*** Are you sure all input files are using the same reference?\n");
					print "Failed.\n";
					exit 1;
				}
				else{
					# Storing ref base in hash table; use loci as key, and value is array of base, input file name, and file line number.
					# The file name and line number help us with reporting mismatched loci base collisions.
					$refBaseTable{$realLoci}=["$3","$inFilename","$i"];

					# Do report - increment SNP count.
					$reportHash{"$currentStrain"}[0]++;

				}
      }
			else{
        print_debug("Excluded loci $realLoci.\n");

				# Do report - increment exclusion count.
				$reportHash{"$currentStrain"}[3]++;

				# Do tally in the exclusion report.
				if ( exists($exclusionReportHash{"$currentStrain"}) ){
					$exclusionReportHash{"$currentStrain"}{"$exclude_name"}++;
				}
				else{
					$exclusionReportHash{"$currentStrain"}{"$exclude_name"} = 1;
				}	
			}
    }
    else{
      print_all("*** ERROR: k28.out file format not recognized. Got \"$_\" at line $i.\n");
			print "Failed.\n";
      exit 1;
    }
  }
	$currentStrain="";
  close($infile);

  print_all("Done.\n\n");
}



#############
# FUNCTIONS
#############
sub print_debug(){
  # Prints output to STDOUT and log file if debug flag is set.  Otherwise nothing.
  # If the quiet function is also set, then only log to file.
  #
  # Parameters:
  # $1 = Text to print.
  #
  if ($debug eq "Y"){
    if ($quiet eq "N"){
      print "$_[0]";
    }
    print $logfile "$_[0]";
  }
}

sub print_all(){
  # Silly MUX'ed function to write output to both standard out and logfile in one call.
  # If the quiet function is set, then only log to file.
  #
  # Parameters:
  # $1 = Text to print.
  #
  if ($quiet eq "N"){
    print "$_[0]";
  }
  print $logfile "$_[0]";
}

