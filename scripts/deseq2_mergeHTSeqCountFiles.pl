#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

#######################################################################################################################
### Manjusha Pande, mpande@umich.edu, 01/15/14
### The script merges second columns of all files based on values in col 1, designed to merge the HTSeq counts from individual samples
### Usage: perl mergeCountFiles.pl <inDir>
### Example: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/mergeCountFiles.pl /ccmb/BioinfCore/Projects/Manjusha/Rsq_test_v3/SEnoNovelTranscripts/HTSeq
#######################################################################################################################
my $inDir = $ARGV[0];

my @inFiles = glob($inDir . "/Sample*_counts.txt");
my %gene2count = ();
my %countStats = ();
my @header = ("GeneId");
my @header2 = ("Readcounts");
my @feature_readcounts = ();
my @total_readcounts = ();

#print Dumper @inFiles;

foreach my $f (@inFiles) {
	my $feature_readcount = 0;
	my $nonfeature_readcount = 0;
	my $total_readcount = 0;
	my @tmp1 = split ("/", $f);
	my $sampleName = $tmp1[$#tmp1];
	$sampleName =~ s/_counts.txt//;
	push (@header, $sampleName);
	push (@header2, $sampleName);
	open (INFILE, "<".$f) or die "Cannot open $f $!\n";
	print "Processing file $f... \n";
	my @stats = (qr/no_feature/, qr/ambiguous/, qr/too_low_aQual/, qr/not_aligned/, qr/alignment_not_unique/);
	my $count = 1;
	while (<INFILE>) {
		chomp $_;
		my @tmp = split (/\t/, $_);
		if ($tmp[0] ~~ @stats) {
			push (@{$countStats{$tmp[0]}}, $tmp[1]); 
			$nonfeature_readcount += $tmp[1];
	
		} else {
			if ($tmp[0] eq "") { $tmp[0] = "-"; }
			push (@{$gene2count{$tmp[0]}}, $tmp[1]);
			$feature_readcount += $tmp[1];
			#if ($count <= 5 ) {print "Gene id is $tmp[0], read count is $tmp[1]\n";}
		}
		$count++;
	}
	close INFILE;
	$total_readcount = $feature_readcount + $nonfeature_readcount;
	push (@total_readcounts, $total_readcount);
	#print "Total feature readcount is $feature_readcount \n";
	push (@feature_readcounts, $feature_readcount);
}

#print Dumper @feature_readcounts;
push(@{$countStats{"feature_assigned"}},@feature_readcounts);
push(@{$countStats{"Total"}},@total_readcounts);

#print Dumper @header;
#print Dumper %gene2count;
#print Dumper %countStats;

open (OUTFILE, "> $inDir/HTSeq_counts.txt") or die "Cannot open $inDir/HTSeq_counts.txt $!\n";
my $header = "";
foreach my $h (@header) {
    $header .= "$h\t"
}
$header =~s/\t$//;
print OUTFILE "$header\n";

my @k = sort(keys%gene2count);
foreach my $k (@k) {
	print OUTFILE "$k";
	foreach my $i (@{$gene2count{$k}}) {
		print OUTFILE "\t$i";
	}
	print OUTFILE "\n";
}
close OUTFILE;

open (OUTFILE2, "> $inDir/HTSeq_counts_stats.txt") or die "Cannot open $inDir/HTSeq_counts_stats.txt $!\n";

$header = "";
foreach my $h (@header2) {
    $header .= "$h\t"
}
$header =~s/\t$//;
print OUTFILE2 "$header\n";

my @keys = sort(keys%countStats);
foreach my $k (@keys) {
	if($k ne "Total") {
		print OUTFILE2 "$k";
		foreach my $i (@{$countStats{$k}}) {
			print OUTFILE2 "\t$i";
		}
		print OUTFILE2 "\n";
	}
}
print OUTFILE2 "\n";
print OUTFILE2 "Total";
foreach my $i (@{$countStats{"Total"}}) {
	print OUTFILE2 "\t$i";
}
print OUTFILE2 "\n";

close OUTFILE2;

print "All done!\n";








	


	
	

