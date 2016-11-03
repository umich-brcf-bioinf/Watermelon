#!/usr/bin/perl
### Manjusha Pande, mpande@umich.edu, 11/11/13
### Get alignment rates for each sample
### Usage: perl getQCMetrics.pl <path/to/AlignSummaryFile> <path/to/OutputFile> <sampleName>


## abhasi, 06/16
## made it compatible with Watermelon (snakemake pipeline)

use warnings;
use strict;
use Data::Dumper;

my $inFile = $ARGV[0];  #   03-tophat/Sample_01/Sample_01_align_summary.txt
my $outFile = $ARGV[1]; #outputfile
my $sample= $ARGV[2];

if (-e $inFile) {
    open (QCFILE, ">>$outFile") || die "Can not open $outFile $!\n";
     print QCFILE "$sample\t";
    getAlignmentRate ($inFile);
    print QCFILE "\n";
	}
	else {
	print QCFILE "Alignment summary file: $inFile is not available.\n";
    die;
    }
    

sub getAlignmentRate {
	my $file = $_[0];
	my $isLeft = 0;
	my $isRight = 0;
	open (INFILE, "<$file") || die "Can not open $file. $!\n";
	while (<INFILE>) {
		chomp $_;
		if ($_ =~ /Left reads/) {
			$isLeft = 1;
		} elsif ($_ =~ /Right reads/) {
			$isRight = 1;
		}
		if ($isLeft && !$isRight && $_ =~ /Input/) {
			$_ =~ m/(\d+)/; 
			print QCFILE "$1 \t";
		}
		#elsif ($isRight && $_ =~ /Input/) {
			#$_ =~ m/(\d+)/;
			#print QCFILE "$1 \t";
		#}
		elsif (!$isLeft && !$isRight && $_ =~ /Input/){
			$_ =~ m/(\d+)/;
			print QCFILE "$1 \t";
		}
		if ($_ =~/% overall/) { #abhasi: 5/9/16
			$_ =~ s/\soverall.*$//;
			print QCFILE $_;
		}
	}
}
close QCFILE;