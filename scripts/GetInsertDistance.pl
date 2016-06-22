#! /usr/bin/perl

### Manjusha Pande, mpande@umich.edu, 12/05/13
### This script gets a random subsample of the paired .fastq files
###             aligns the reads to reference genome, sorts the aligned files and 
###             gets the insert size distribution and writes to <Sample>_insert_dist_file
### Usage: $ perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/GetInsertDistance.pl <readsDir> <fraction>  <refGenomeBTindex> <nProc> <readFiles> <runInfoDir> [cmdFile]
### Example: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/GetInsertDistance.pl /ccmb/home/mpande/RNA-seq_pipeline/PE_inner_dist_tophat_test 0.01 /ccmb/BioinfCore/Projects/Manjusha/Genomes/iGenomes/Caenorhabditis_elegans_UCSC_ce10/UCSC/ce10/Sequence/Bowtie2Index/genome 20
### Example: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/GetInsertDistance.pl /ccmb/BioinfCore/Projects/Manjusha/PE_tophat_test/reads 0.001 /ccmb/BioinfCore/Projects/Manjusha/Genomes/iGenomes/Caenorhabditis_elegans_UCSC_ce10/UCSC/ce10/Sequence/Bowtie2Index/genome 20
#########################################################################################################################
use strict;
use warnings;
use Data::Dumper;

my $projectDir = $ARGV[0];
my $fraction = $ARGV[1];
#my $refGenome =  $ARGV[2];
my $refGenomeBTindex = $ARGV[2];
my $nProc = $ARGV[3];
my $readFiles = $ARGV[4];
my $runInfoDir = $ARGV[5];

my $n_bg_processes = 0;

my $readsDir = "$projectDir/reads";
my $outDir = "$readsDir/subSample";
unless (-e $outDir || mkdir ($outDir, 0775)) {
	die "Unable to create $outDir \n";
};

my $GMD;

if ($ARGV[6]) {
	$GMD = $ARGV[6];
	open (CMDFILE, ">>" . $GMD . ".sh");	
} else { 
	my $libDir =  "/ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts";
	my $runErrorLog = "$readsDir/GMD_error.log";
	my $GMD = $outDir . "/GetInsertDist";
	open (CMDFILE, ">" . $GMD . ".sh");
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

open (READFILES, "<$readFiles") || die "Cannot open $readFiles for reading $!\n";
### Getting subset of fastq.gz files...
my @inFiles1 = ();
while (<READFILES>) {
	chomp $_;
	push (@inFiles1, $_);
}
#print Dumper @inFiles1;
print CMDFILE "\necho \"Getting subsample of fastq.gz files...\";\n";
#print "Getting subset of fastq.gz files...\n";
my %sample_pairs = ();

foreach my $f (@inFiles1) {
	my @tmp = split (/\//, $f);
	my $sample = $tmp[$#tmp];
	$sample =~ s/_R\d+.fastq.*.gz//;
	push (@{$sample_pairs{$sample}}, $f);
}
#print Dumper %sample_pairs;

my @sample_pairs = sort (keys %sample_pairs);
#print Dumper @sample_pairs;

#=pod
my $count1 = 0;
foreach my $s (@sample_pairs) {
	$count1++;
	#print "Processing $s... \n";
	my $cmd = "python /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/subsample.py $fraction " . @{$sample_pairs{$s}}[0] . " " . @{$sample_pairs{$s}}[1] . " ".  $outDir . "/". $s . "_R1_subsample$fraction.fastq " . $outDir . "/". $s . "_R2_subsample$fraction.fastq &";
	print CMDFILE "$cmd\n";
	$n_bg_processes++;
	print CMDFILE "p" . $n_bg_processes . "=\$!\n";
	if ($count1 == $nProc) {
		#print CMDFILE "\nwait;\n\n";
		print CMDFILE "wait";
		for my $bg (1..$n_bg_processes) { print CMDFILE " \$p" . $bg};
		print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
		$n_bg_processes = 0;
		$count1 = 0;
	}
	#system ("$cmd");
}
#print CMDFILE "\nwait;\n"; 
print CMDFILE "wait";
for my $bg (1..$n_bg_processes) { print CMDFILE " \$p" . $bg};
print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
$n_bg_processes = 0;

#=cut### Aligning the subsample files...
my $alnDir = "$outDir/bowtie2";
unless (-e $alnDir || mkdir ($alnDir, 0775)) {
	die "Unable to create $alnDir \n";
};
#=pod

print CMDFILE "\necho \"Aligning the subsample files...\";\n";

my $count2 = 1;
my $proc_count = $nProc;
my $sample_count = scalar(@sample_pairs);

my @inFiles3 = ();

foreach my $s (@sample_pairs ) {
	my $p = int($proc_count / $sample_count);
	#my $cmd = "bwa mem -t 1 -R \"\@RG	ID:ID	LB:LIB	SM:$s	PL:ILLUMINA\" -M " . $refGenomeSeq . " " . $outDir . "/". $s . "_R1_subsample" . $fraction . ".fastq.gz " . $outDir . "/". $s . "_R2_subsample" . $fraction . ".fastq.gz  > " . $alnDir . "/" . $s . "_subsample" . $fraction . ".sam &";
	if ($p >= 1){ # 1 or more processor/s per sample
		my $cmd = "bowtie2 -p $p -x $refGenomeBTindex -1 $outDir/${s}_R1_subsample$fraction.fastq -2 $outDir/${s}_R2_subsample$fraction.fastq -S $alnDir/${s}_subsample$fraction.sam > $alnDir/${s}_subsample_alignment_summary.txt 2>&1 &";
		print CMDFILE "$cmd\n";
		$n_bg_processes++;
		print CMDFILE "p" . $n_bg_processes . "=\$!\n";
		push (@inFiles3,  $alnDir . "/". $s . "_subsample" . $fraction . ".sam");
		$count2++;
		$proc_count = $proc_count-$p;
	}
	else { ## fewer processors than samples
		if($count2 <= $nProc) { ## not all processors used up yet
			my $cmd = "bowtie2 -p 1 -x $refGenomeBTindex -1 $outDir/${s}_R1_subsample$fraction.fastq -2 $outDir/${s}_R2_subsample$fraction.fastq -S $alnDir/${s}_subsample$fraction.sam > $alnDir/${s}_subsample_alignment_summary.txt 2>&1 &";
			print CMDFILE "$cmd\n";
			$n_bg_processes++;
			print CMDFILE "p" . $n_bg_processes . "=\$!\n";
			push (@inFiles3,  $alnDir . "/". $s . "_subsample" . $fraction . ".sam");
			$count2++;
			$proc_count = $proc_count-1;
		}
		if ($proc_count == 0) {
			#print CMDFILE "\nwait;\n\n";
			print CMDFILE "wait";
			for my $bg (1..$n_bg_processes) { print CMDFILE " \$p" . $bg};
			print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
			$n_bg_processes = 0;
			#print "wait!\n";
			$proc_count = $nProc;
			$count2 = 1;
		}
	}
	$sample_count--;
}
#print CMDFILE "\nwait;\n";
print CMDFILE "wait";
for my $bg (1..$n_bg_processes) { print CMDFILE " \$p" . $bg};
print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
$n_bg_processes = 0;
#=cut
#=pod
### Sorting aligned files...
#my @inFiles3 = glob( $alnDir . "/*_subsample*.sam");
#print Dumper @inFiles3;
print CMDFILE "\necho \"Sorting aligned files...\";\n";
#print "Sorting aligned files...\n";
my $count3 = 0;
my @inFiles4 = ();
foreach my $f (@inFiles3) {
	$count3 += 3;
	$f =~ s/.sam$//;
	my $cmd = "java -Xmx4g -Djava.io.tmpdir=$alnDir/tmp -jar /opt/bioinf/picard/1.77/SortSam.jar SO=coordinate INPUT=$f.sam OUTPUT=${f}_sorted.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &";
	print CMDFILE "$cmd\n";
	$n_bg_processes++;
	print CMDFILE "p" . $n_bg_processes . "=\$!\n";
	push (@inFiles4, "${f}_sorted.bam");
	if ($count3 >= $nProc) {
		#print CMDFILE "\nwait;\n\n";
		print CMDFILE "wait";
		for my $bg (1..$n_bg_processes) { print CMDFILE " \$p" . $bg};
		print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
		$n_bg_processes = 0;
		$count3 = 0;
	}
	#system ("$cmd")
}
#print CMDFILE "\nwait;\n";
print CMDFILE "wait";
for my $bg (1..$n_bg_processes) { print CMDFILE " \$p" . $bg};
print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
$n_bg_processes = 0;
#=cut
#print Dumper @inFiles4;

### Getting insert size metrics
##my @inFiles4 = glob( $alnDir . "/*_subsample*_sorted.bam");
#print Dumper @inFiles4;
print CMDFILE "\necho \"Getting insert size metrics...\";\n";
#print "Getting insert size metrics...\n";
my $count4 = 0;
foreach my $f (@inFiles4) {
	$count4 += 3;
	my @tmp = split (/\//, $f);
	my $sample = $tmp[$#tmp];
	$sample =~ m/(Sample_[A-Za-z]*\d+.*)_subsample/;
	my $s = $1;
	#print "Sample name is $s\n";
	my $cmd = "java -Xmx4g -jar -Djava.io.tmpdir=$alnDir/tmp -jar /opt/bioinf/picard/1.77/CollectInsertSizeMetrics.jar INPUT=$f HISTOGRAM_FILE=$alnDir/${s}_insert_histo.pdf OUTPUT=$alnDir/${s}_insert_dist_file &";
	print CMDFILE "$cmd\n";
	$n_bg_processes++;
	print CMDFILE "p" . $n_bg_processes . "=\$!\n";
	if ($count4 >= $nProc) {
		#print CMDFILE "\nwait;\n\n";
		print CMDFILE "wait";
		for my $bg (1..$n_bg_processes) { print CMDFILE " \$p" . $bg};
		print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
		$n_bg_processes = 0;
		$count4 = 0;
	}	
	#system ("$cmd")
}
print CMDFILE "wait";
for my $bg (1..$n_bg_processes) { print CMDFILE " \$p" . $bg};
print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
$n_bg_processes = 0;
close CMDFILE;
#print "nohup sh $GMD.sh 2>&1 > ${GMD}.log &\n"
#system("nohup sh $GMD.sh 2>&1 > ${GMD}.log &");

#print "I am back in the Perl script.\n";
