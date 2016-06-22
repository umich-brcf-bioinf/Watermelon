#! /usr/bin/perl 

### Manjusha Pande, mpande@umich.edu, 12/05/13
### This script calulates distribution of read lengths for each fastq.gz file and and outputs that to a <Sample>_R*.fastq.gz_readCountandLength.txt file 
### Usage: $ perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/CalculateReadLength.pl <readsDir> <nProc> <runInfoDir> <readFiles>[cmdFile]
### Example: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/CalculateReadLength.pl /ccmb/home/mpande/RNA-seq_pipeline/PE_inner_dist_tophat_test/reads 1 /ccmb/home/mpande/RNA-seq_pipeline/PE_inner_dist_tophat_test/runInfo
### Example: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/CalculateReadLength.pl /ccmb/BioinfCore/Projects/Hu_826_pathu_mceachin_mpande_RS2/reads 20 /ccmb/BioinfCore/Projects/Hu_826_pathu_mceachin_mpande_RS2/runInfo
#########################################################################################################################
use strict;
use warnings;
use Data::Dumper;

my $readsDir = $ARGV[0]; ### this dir has .fastq.gz files
my $nProc = $ARGV[1];
my $runInfoDir = $ARGV[2];
my $readFiles = $ARGV[3];
my $RCL;
my $outDir = "$ARGV[0]/readCountLength"; ### this dir is inside readsDir to store the read length and read count per fastq.gz file
unless (-e $outDir || mkdir ($outDir, 0775)) {
	die "Unable to create $outDir \n";
};

if ($ARGV[4]) {
	$RCL = $ARGV[4];
	open (CMDFILE, ">>" . $RCL . ".sh");	
} else {
	my $libDir =  "/ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts";
	my $runErrorLog = "$readsDir/RCL_error.log";
	$RCL = $outDir . "/GetReadCountandLength";
	open (CMDFILE, ">" . $RCL . ".sh");
	print CMDFILE "#!/bin/sh\n";	
	print CMDFILE "\nPROGNAME=\$(basename \$0)
	#echo \$PROGNAME
	function error_exit
	{
		echo \"\${PROGNAME}: \${1:-\"Unknown Error\"}\" >> $runErrorLog 2>&1
      		ERROR=\"\${1:-\"Unknown Error\"} in \${PROGNAME}\" 
		DATE=\$(date +\"%a, %B %e, %Y\")
		TIME=\$(date +\"%I:%M %p\")
		perl ${libDir}/send_email.pl \"Your RNA-seq run terminated (exit status = 1) on \${DATE} at \${TIME}. 
\${ERROR}.\" \"\" $runErrorLog
		exit 1
	}
	function error_log
	{
		echo \"\${PROGNAME}: \${1:-\"Unknown Error\"}\" >> $runErrorLog 2>&1	
	}\n";
}
print CMDFILE "\necho \"Calculating read count and mean read length for fastq.gz files..\";\n";
my $n_bg_processes = 0;
if ($nProc >=3 ) {
	$nProc = $nProc / 3;
}

open (READFILES, "<$readFiles") || die "Cannot open $readFiles for reading $!\n";
### Getting subset of fastq.gz files...
my @inFiles = ();
while (<READFILES>) {
	chomp $_;
	push (@inFiles, $_);
}

#print Dumper @inFiles;
#my %readCounts = ();

my $sample_count = scalar(@inFiles);
#print "Sample count is $sample_count\n";
my $count = 1;
my @subsetInFiles = ();
while (scalar(@inFiles)){
	#print "New round\n";
	if (scalar(@inFiles) >= $nProc) {
		@subsetInFiles = splice (@inFiles, 0, $nProc);
	} else {
		@subsetInFiles = splice (@inFiles, 0, scalar(@inFiles));
	}
	#print Dumper @subsetInFiles;
#=pod
	foreach my $f (@subsetInFiles) {
		my @tmp = split (/\//, $f);
		my $sample = $tmp[$#tmp];
		#print "Processing $sample [$count of $sample_count]... \n";
		#print "`zcat $f | awk '{if(NR%4==2) print length(\$1)}' | sort | uniq -c &`\n";
		#push (@{$readCounts{$sample}}, `zcat $f | awk '{if(NR%4==2) print length(\$1)}' | sort | uniq -c`); 
		my $cmd = "zcat $f | awk '{if(NR%4==2) print length(\$1)}' | sort | uniq -c > $outDir/${sample}_readCountandLength.txt &";
		#system ("zcat $f | awk '{if(NR%4==2) print length(\$1)}' | sort | uniq -c > $outDir/${sample}_readCountandLength.txt &");
		print CMDFILE "$cmd\n";
		$n_bg_processes++;
		print CMDFILE "p" . $n_bg_processes . "=\$!\n";
		$count++;
	}
	#print CMDFILE "\nwait;\n\n";
	print CMDFILE "wait";
	for my $bg (1..$n_bg_processes) { print CMDFILE " \$p" . $bg};
		print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
	$n_bg_processes = 0;
#=cut
}
my $cmd = "perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/GetReadLength.pl $outDir $runInfoDir";
print CMDFILE "$cmd\n";
print CMDFILE "status=\$? 
if [ \$status != 0 ] ; then error_exit \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n";
close CMDFILE;
#print "nohup sh $RCL.sh 2>&1 > ${RCL}.log &\n"
#system("nohup sh $RCL.sh 2>&1 > ${RCL}.log &");


##############################################################################################################################
=pod
### The following code works when one file is processed at a time using the awk command in backticks
my $pe = 0;
open (OUTFILE, ">${runInfoDir}/readCountandLength.txt") || die "Cannot open ${runInfoDir}/readCountandLength.txt for writing\n";
#print Dumper %readCounts;
#print Dumper keys(%readCounts);
my %sample_pairs = ();
for my $k (keys %readCounts) {
	my $sample_name = $k;
	$sample_name =~ s/_R.*//;
	if ($k =~ /R1/) {
		@{$sample_pairs{$sample_name}}[0] = $k;
	} elsif ($k =~ /R2/) {
		@{$sample_pairs{$sample_name}}[1] = $k;
		$pe = 1;
	}		
}
if ($pe) {
	print OUTFILE "#Sample\tRead Count\tMean Read Length R1\tMean Read Length R2\n";
} else { print OUTFILE "#Sample\tRead Count\tRead Length\n";
}
#print Dumper %sample_pairs;
#=pod
for my $k (keys %sample_pairs) {	
	print OUTFILE "$k\t";
	my @readLengths1 = @{$readCounts{@{$sample_pairs{$k}}[0]}};
	#print "@readLengths1[0..$#readLengths1]\n";
	my $sum1 = 0;
	my $readCount1 = 0;
	foreach my $r1 (@readLengths1) {
		my @r1 = $r1 =~ m/(\d+)/g;
		$sum1 = $sum1 + ($r1[0] * $r1[1]);
		$readCount1 = $readCount1 + $r1[0];
	}
	my $mean_readLength1 = $sum1 / $readCount1;
	print OUTFILE "$readCount1\t$mean_readLength1";

	if (defined (@{$sample_pairs{$k}}[1])) {
		my @readLengths2 = @{$readCounts{@{$sample_pairs{$k}}[1]}};
		#print "@readLengths2[0..$#readLengths2]\n";
		my $sum2 = 0;
		my $readCount2 = 0;
		foreach my $r2 (@readLengths2) {
			my @r2 = $r2 =~ m/(\d+)/g;
			$sum2 = $sum2 + ($r2[0] * $r2[1]);
			$readCount2 = $readCount2 + $r2[0];
		}
		my $mean_readLength2 = $sum2 / $readCount2;
		print OUTFILE "\t$mean_readLength2";		
	}
	print OUTFILE "\n";
}
=cut

#print "All done!\n";
