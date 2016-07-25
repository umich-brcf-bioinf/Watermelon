#!/usr/bin/perl
### Manjusha Pande, mpande@umich.edu, 11/11/13
### Get alignment rates for each sample
### Usage: perl getQCMetrics.pl <path/to/alnDir> <path/to/sample_names_file> <path/to/outDir> [runInfoFile]
### Example: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/getQCMetrics.pl /ccmb/BioinfCore/Projects/Hu_826_pathu_mceachin_mpande_RS2/tophat /ccmb/BioinfCore/Projects/Hu_826_pathu_mceachin_mpande_RS2/runInfo/Sample_list.txt /ccmb/BioinfCore/Projects/Hu_826_pathu_mceachin_mpande_RS2/runInfo


## abhasi, 06/16
## made it compatible with Watermelon (snakemake pipeline)

use warnings;
use strict;
use Data::Dumper;

my $alnDir = $ARGV[0];
my $sampleList = $ARGV[1];
my $runInfoDir = $ARGV[2];

#unless (-e $runInfoDir || mkdir ($runInfoDir, 0775)) {
#	die "Unable to create $runInfoDir \n";
#	};

my @samples = ();
open (SAMPLELIST, "<$sampleList") || die "Can not open $sampleList $!\n";
while (<SAMPLELIST>) {
	unless ($_ =~/^#/) {
		chomp $_;
		push (@samples, $_);
	}
}	
close SAMPLELIST;
my $outFile;
if ($ARGV[3]) {
	$outFile =  $ARGV[3];
}
else {$outFile = $runInfoDir . "/Align_summary_all.txt";}

open (QCFILE, ">>$outFile") || die "Can not open $outFile $!\n";
print QCFILE "\nPer sample read counts and alignment rates:\n";
print QCFILE "Sample Name\tRead Count\tOverall Alignment Rate\n"; 

my $inFile;
foreach my $s (@samples) {
	print QCFILE "$s\t";
	$inFile = "$alnDir/$s/". $s."_align_summary.txt";
	if (-e $inFile) {
		getAlignmentRate ($inFile);
		print QCFILE "\n";
	}
	else { print QCFILE "Alignment summary file $alnDir/$s/".$s."_align_summary.txt not available.\n";}
}
print QCFILE "\n";

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
	
