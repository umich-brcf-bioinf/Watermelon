#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

### Manjusha Pande, mpande@umich.edu, 10/31/2013.
### Original code is from http://cpansearch.perl.org/src/JMCNAMARA/Spreadsheet-WriteExcel-2.39/examples/tab2xls.pl modified to use Excel::Writer::XLSX instead of Spreadsheet::WriteExcel
### Usage:  perl tab2xls.pl <path/to/inDir> <path/to/outDir> [runInfoFile]
### Example:  perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/tab2xls.pl /ccmb/home/mpande/RNA-seq_pipeline/Sample_cuffdiff /ccmb/home/mpande/RNA-seq_pipeline/Sample_cuffdiff/report
### Example:  perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/tab2xls.pl /ccmb/BioinfCore/Projects/Manjusha/Rsq_test_v2/noNovelTranscripts/cuffdiff/E7.5_WT_v_E7.5_mutant /ccmb/BioinfCore/Projects/Manjusha/Rsq_test_v2/noNovelTranscripts/report 

use Excel::Writer::XLSX;

my $reportDir = $ARGV[1];
unless (-e $reportDir || mkdir ($reportDir, 0775)) {
	die "Unable to create $reportDir \n";
	};
my $inDir = $ARGV[0];
#if ($ARGV[0] eq ".") {
#	$inDir = `pwd`;
#}	
#my @letters = 'A'..'XFD';
#print Dumper @letters[0..26];

my @inFiles = glob($ARGV[0] . "/*{gene,isoform}*foldchange*_annotated.txt");
push (@inFiles, glob($ARGV[0] . "/*DAVIDchartReport.txt"));
my $runInfoFile;
my $xlFileName = '';

my @tmp = split("/", $inDir); 
if ($tmp[$#tmp] =~ /_v_/)
{
	$xlFileName = $tmp[$#tmp] . "_DE.xlsx";
}
else { 
	$xlFileName = $tmp[$#tmp] . ".xlsx";
}

if ($ARGV[2]) {
	$runInfoFile =  $ARGV[2];
	push (@inFiles, $runInfoFile);
	my @tmp2 = split(/\//, $runInfoFile);
	my $filename = $tmp2[$#tmp2];
	$filename =~ s/.txt//;
	$xlFileName = $filename . ".xlsx";
} 

#print Dumper @inFiles;

my $xlFile = $reportDir . "/" . $xlFileName;

#print $xlFile . "\n";
# Create a new Excel workbook
my $workbook = Excel::Writer::XLSX->new($xlFile);


for my $f (0..$#inFiles) 
{
	my $gene = 0;
	my $enrichment = 0;
	my @header = ();
	my $headerRow = 0;
	my $geneIdCol = 0;
	my $GenesCol = 0;
	my $hyperlink_format = $workbook->add_format(
   	 color     => 'blue',
   	 underline => 1,
	);
	my $worksheet;
	#my @tmp = split("/", $inFiles[$f]);
	#my $filename = $tmp[$#tmp];
	# Open the tab delimited file
	open (TABFILE, $inFiles[$f]) or die "Can not read $inFiles[$f] $!";
	if ($inFiles[$f] =~ /gene/) { 
		$worksheet = $workbook->add_worksheet('genes');
		$gene = 1;
		#print "Processing gene worksheet..\n";

	} elsif ($inFiles[$f] =~ /isoform/) { 
		$worksheet = $workbook->add_worksheet('isoforms');
	} elsif ($inFiles[$f] =~ /DAVIDchartReport/) { 
		$worksheet = $workbook->add_worksheet('enrichment');
		$enrichment = 1;
	} elsif ($inFiles[$f] =~ /runinfo/) { 
		$worksheet = $workbook->add_worksheet('runinfo');
	}

	# Row and column are zero indexed
	my $row = 0;
	my $col;

	while (<TABFILE>) {
		my $header = 0;		
		my $link; 
		my $formula;
    		chomp;	
		if ($_ =~ /Entrez GeneId|#Category/) {
			$header = 1;
			$headerRow = $row+1;
			$worksheet->freeze_panes($headerRow, 0);
			#print "Found header\n";
		}	
    		# Split on single tab
   		my @Fld = split('\t', $_);
		if ($header) {
			@header = @Fld; 
			for my $i (0..$#header) {
				if ($header[$i] =~ /Entrez GeneId/) {
					$geneIdCol = $i;
					#print "Gene id column is $geneIdCol \n";
				} elsif ($header[$i] =~ /Genes/) {
					$GenesCol = $i;
				}
			}
			#if ($gene) { splice (@Fld, $geneIdCol+2, 0, "NCBI Gene Record");}
			if ($gene) { $Fld[$geneIdCol] = "NCBI Gene";}
		}
		elsif ( $gene && scalar (@header) && $Fld[$geneIdCol]) {
			#print "Converting gene ids to hyperlinks..\n";		
			if (scalar (@Fld) &&  scalar (@Fld) >= scalar (@header) -1) {
				$link = 'http://www.ncbi.nlm.nih.gov/gene/?term=' . $Fld[$geneIdCol];  ### url
				#splice (@Fld, $geneIdCol+2, 0, $link);
				#$formula =  "=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/gene/?term=" . $Fld[$geneIdCol] . "[uid]\", $Fld[$geneIdCol])"; ### formula
				#$Fld[$geneIdCol] = $formula;
			}	
		} #elsif ( $enrichment && scalar (@header) && $Fld[$GenesCol]) { 
			#$link = 'http://www.ncbi.nlm.nih.gov/gene/?term=' . $Fld[$GenesCol];
		#}
		
    		$col = 0;
    		foreach my $token (@Fld) { 	
			$token =~ s/^\s+|\s+$//; ## remove spaces at the beginning or end
			if($gene && scalar (@header)  && !$header && scalar (@Fld) && scalar (@Fld) >= (scalar (@header)-1) && $col == $geneIdCol && $Fld[$geneIdCol]) { ## sheet is gene, header is defined, line is not a header or blank or a comment, gene id is not blank
				#$worksheet->write_url($row, $col, $link);
				#$worksheet->write_formula($row, $col, $token, undef, $hyperlink_format);
				#$Fld[$geneIdCol] =~ s/^\s+|\s+$//; ## remove spaces at the beginning or end
				$Fld[$geneIdCol] = int($Fld[$geneIdCol]);
				$worksheet->write($row, $col, $link, $hyperlink_format, $Fld[$geneIdCol]);
			} 
			#if ($enrichment && scalar (@header)  && !$header && scalar (@Fld) && scalar (@Fld) >= (scalar (@header)-1) && $col == $GenesCol && $Fld[$GenesCol]) {			
				#print "$link\n";
				#$worksheet->write($row, $col, $link, $hyperlink_format, $Fld[$GenesCol]);
			#} 
			else {
        			$worksheet->write($row, $col, $token);
			}
       		$col++;
    		}
    	$row++;
	}
	close TABFILE;

	#my $worksheet_name = $worksheet->get_name();
	#print "The worksheet is $worksheet_name, ";
	
### Following code attempts to get number of differentially expessed genes, does not work at this time
=pod
	my $col_letter = $letters[$col];
	
	my $formula =  "7, ".$col.", \'=COUNTIF(".${letters[$col-1]}."7:".${letters[$col-1]}.$row.", \"=YES \")\'\n";
	my $formula =  "\'".$col_letter."7\', \'=COUNTIF(".${letters[$col-1]}."7:".${letters[$col-1]}.$row.", \"YES \")\'\n";
	print "Formula is $formula\n";
	$worksheet->write_formula($formula);
	$worksheet->write($formula);
=cut
}
$workbook->close();




