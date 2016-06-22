#!/usr/bin/perl

use strict;
use warnings;

#
#  Maintainers: Rich McEachin <mceachin@umich.edu>; Ashwini Bhasi <abhasi@umich.edu>
#  Date: 8 Feb 2012
#  Version: 6
#  Status: Development 
#  This version automatically creates the output file name and creates a .txt (tab separated text) file as output
#  Compatible with Cuffdiff Version: v1.1.x
#  Description:
#  This program takes cuffdiff output file: gene_exp.diff or isoform_exp.diff, and adds the following columns:
#  (1) Fold Change (Calculated by 2 ** log2fold change (column log2fold change in gene_exp.diff)
#  (2)  Diff Exp (YES/NO) : Calculated by checking IF [(FDR <=0.05) AND (Fold Change >= Threshold or Fold Change <= 1/Threshold) AND( Test Status eq OK)]
#  FDR is the column "qvalue" and Test Status is column "Status" in gene_exp.diff; Threshold = the user input $threshold, and it is assumed to be >= 1.0.  
#  By default, the program checks for both up and down regulated transcripts.
#  It also sorts the results on Diff Exp (YES before NO) & Fold Change (highest to lowest)
#  Import from the command line: cuffdiff output file ($infile1), the text (no spaces) you'd like to start the filename with ($text), 
#  and the threshold for defining differential expression ($threshold)
#  Command: perl Cuffdiff_out_format_v5.pl $infile $text $threhsold
#  For example: "perl Cuffdiff_out_format_v5.pl gene_exp.diff goldstrohm_puf 2.0"
#

my ($infile, $text, $threshold) = @ARGV;
my $outfile;

if (($infile eq "") or ($threshold eq "") ){
	print "\nWARNING! Please enter input parameters in the following order:\n\n";
	print "\"Cuffdiff_out_format.pl Input_File_name Identifying_text Threshold_Value\"\n";
	exit;
}

#Open in and out files
open IN1, "$infile" or die "Cannot open $infile $!\n";
open OUT, ">$text.foldchange.$threshold.txt" or die "Cannot open $outfile $!\n";

#Set up header for out file
print OUT "Select Differentially Expressed Transcripts Based on Three Criteria\n";
print OUT "Test Status OK\n";
print OUT "FDR (q-value) <= 0.05\n";
print OUT "Fold Change >= Threshold or <= 1/Threshold, where Threshold = $threshold \n\n";
print OUT "test_id\t gene_id\t gene\t locus\t sample_1\t sample_2\t status\t value_1\t value_2\t log2(fold_change)\t test_stat\t p_value\t q_value\t significant\t Fold Change\t Diff Exp\n";
my $i=0;
my @AOA;

 while (my $line1 = <IN1>) {
 unless ($line1=~/^test_id/)	{
 
		($AOA[$i][0], $AOA[$i][1], $AOA[$i][2],$AOA[$i][3],$AOA[$i][4],$AOA[$i][5],$AOA[$i][6], $AOA[$i][7], $AOA[$i][8],$AOA[$i][9],$AOA[$i][10],$AOA[$i][11],$AOA[$i][12],$AOA[$i][13])= split(/\t/,$line1);
			chomp $AOA[$i][13];
			
			#calculating Fold Change
			$AOA[$i][14] = 2**($AOA[$i][9]);
			
			#calculating Diff Exp 
			if (($AOA[$i][14] >= $threshold or $AOA[$i][14] <= 1/$threshold) && ($AOA[$i][12]<= 0.05)&& ($AOA[$i][6] eq "OK")){
				$AOA[$i][15] = "YES";
			} else {
				$AOA[$i][15] = "NO";
			}
			
			$i++;
		}
}

#sorting results based on Diff Exp & Fold Change
my @AOA_sort = sort{
				$b->[6] cmp $a->[6] || $b->[13] cmp $a->[13] || $b->[15] cmp $a->[15] || $b->[14] <=> $a->[14];
            } @AOA; 
			
foreach my $item1 (@AOA_sort){
  foreach my $item2 (@{$item1}){
    print OUT "$item2 \t";  
  }
  print OUT "\n";
}                   

print "done\n";
close IN1;
close OUT;