### Manjusha Pande, mpande@umich.edu, 12/05/13
### Getting inner mate distance for tophat parameter -r
### Usage: $ perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/GetMateDistance.pl <readsDir> <runInfoDir> 
### Example: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/GetMateDistance.pl /ccmb/BioinfCore/Projects/Manjusha/PE_tophat_test/reads /ccmb/BioinfCore/Projects/Manjusha/PE_tophat_test/runInfo
#########################################################################################################################
use strict;
use warnings;
use Data::Dumper;

my $readsDir = $ARGV[0]; 
my $runInfoDir = $ARGV[1]; ## this dir contains a readCountandLength.txt for the dataset 
my $outDir = $runInfoDir;

my $inDir = "$ARGV[0]/subSample/bowtie2"; ## this dir contains one Sample_<id>_insert_dist_file per sample
### Getting inner mate distance per sample....
my %samples;
open (READLENGTH, "<$runInfoDir/readCountandLength.txt") || die "Cannot open $runInfoDir/readCountandLength.txt for reading $!\n";
while (<READLENGTH>) {
	chomp $_;
	unless ($_ =~ /^#/) {
		my @tmp = split (/\t/, $_);
		push (@{$samples{$tmp[0]}}, @tmp[2..5]);
	}
}
#print Dumper %samples;

my @samples = sort (keys %samples);
#print Dumper @samples;
#=pod

foreach my $s (@samples) {
	#print "The sample is $s\n";
	my $line_count = 0;
	my @header;
	my $median_col;
	my $mean_col;
	my $sd_col;
	my $insertFile = "$inDir/" . $s . "_insert_dist_file";
	open (INFILE, "<$insertFile") || die "Cannot open $insertFile for reading. $!\n";
	while (<INFILE>) {
		if ($_ =~ /^#|^\s/) { next;}
		else {
			chomp $_;
			if ($_ =~ /^MEDIAN/) {
				#print "The line is $_\n";
				@header = split (/\t/, $_);
				#print "size of header is " . scalar (@header) . "\n";
				for my $i (0..$#header) {
					if ($header[$i] eq "MEDIAN_INSERT_SIZE") {
						$median_col = $i;
					} elsif ($header[$i] eq "MEAN_INSERT_SIZE") {
						$mean_col = $i;
					} elsif ($header[$i] eq "STANDARD_DEVIATION") {
						$sd_col = $i;
					}
				}

				$line_count++;
			}

			elsif ($line_count == 1 ){
				#print "The line is $_\n";
				#print "median_col is $median_col\tmean_col is $mean_col\tsd_col is $sd_col\n";
				my @tmp = split (/\t/, $_);
				if (!defined (@{$samples{$s}}[4]) ){
					${$samples{$s}}[4] = $tmp[$median_col];
				}
				if (!defined (@{$samples{$s}}[5]) ){
					${$samples{$s}}[5] = $tmp[$mean_col];
					${$samples{$s}}[5] = sprintf ("%.2f", ${$samples{$s}}[5]);
				}
				if (!defined (@{$samples{$s}}[6]) ){
					${$samples{$s}}[6] = $tmp[$sd_col];
					${$samples{$s}}[6] = sprintf ("%.2f", ${$samples{$s}}[6]);
				}
				#push (@{$samples{$s}}, @tmp[$median_col, $mean_col, $sd_col]);
				$line_count--;		
			}		
			else {next;}
		}
	}			
}
#=cut
#print Dumper %samples;
open (OUTFILE, ">$outDir/innerMateDist.txt") || die "Cannot open $outDir/innerMateDist.txt for writing. $! \n";
print OUTFILE "#Sample\tMean Read Length R1\tSD Read Length R1\tMean Read Length R2\tSD Read Length R2\tMedian Insert Length\tMean Insert Length\tSD Insert Length\tMate Inner Distance\tMate SD\n";
foreach my $s (@samples) {
	print OUTFILE "$s\t"; 
	my $mateInnerDist = ${$samples{$s}}[5] - (${$samples{$s}}[0] + ${$samples{$s}}[2]);
	$mateInnerDist = sprintf ("%.2f", $mateInnerDist);
	#my $mateStdDev = ${$samples{$s}}[6] - (${$samples{$s}}[1] + ${$samples{$s}}[3]);
	my $mateStdDev = ${$samples{$s}}[6]; ## ignoring std dev in read lengths
	$mateStdDev = sprintf ("%.2f", $mateStdDev);
	foreach my $i (@{$samples{$s}}) {
		print OUTFILE "$i\t";
	}
	print OUTFILE "$mateInnerDist\t$mateStdDev\n"; 
}
















