#! /usr/bin/perl 

### Manjusha Pande, mpande@umich.edu, 04/29/14
### Sort and index bam files
### Usage: $ perl sort_index_bams.pl <inDir> <bedoption> <nProc> [runInfoDir]
### Example: $ perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/sort_index_bams.pl /ccmb/BioinfCore/Projects/Lieberman_922_liebermn_mceachin_RS2/tophat 1 10 /ccmb/BioinfCore/Projects/Lieberman_922_liebermn_mceachin_RS2/runInfo
#########################################################################################################################
use strict;
use warnings;
use Data::Dumper;
my $start = time();

my $alnDir = $ARGV[0];
my $createBed = $ARGV[1];
my $nProc = $ARGV[2];
my $libDir = "/ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts";
my $runErrorLog;
my $cmdFile;

if ($ARGV[3]) {
	$cmdFile = $ARGV[3] . "/bam_sort";
} else { 
	$cmdFile = $alnDir . "/bam_sort";
}
if ($createBed) {
	$cmdFile = $cmdFile . "_bed";
}
$runErrorLog = $cmdFile . ".error";
open (CMDFILE, ">" . $cmdFile . ".sh");
print CMDFILE "#!/bin/sh\n";
print CMDFILE "\nPROGNAME=\$(basename \$0)
#echo \$PROGNAME
function error_exit
{
	echo \"\${PROGNAME}: \${1:-\"Unknown Error\"}\" >> $runErrorLog 2>&1
        ERROR=\"\${1:-\"Unknown Error\"} in \${PROGNAME}\" 
	DATE=\$(date +\"%a, %B %e, %Y\")
	TIME=\$(date +\"%I:%M %p\")
	perl ${libDir}/send_email_v1.pl \"Your RNA-seq run terminated (exit status = 1) on \${DATE} at \${TIME}. 
\${ERROR}.\" \"\" $runErrorLog
	exit 1
}
function error_log
{
	echo \"\${PROGNAME}: \${1:-\"Unknown Error\"}\" >> $runErrorLog 2>&1		
}\n";

if ($createBed) {
	print CMDFILE "\necho \"Sorting and indexing .bam files and creating .bed files...\";\n";	
} else { print CMDFILE "\necho \"Sorting and indexing .bam files ...\";\n";}
	
my $outDir1 = "$alnDir/SortedBams";
my $outDir2 = "$alnDir/BedFiles";
unless (-e $outDir1 || mkdir ($outDir1, 0775)) {
	die "Unable to create $outDir1 \n";
}
unless (-e $outDir2 || mkdir ($outDir2, 0775)) {
	die "Unable to create $outDir2 \n";
}

#print "$alnDir\n";
my @samples = glob($alnDir."/Sample*/Sample*accepted_hits.bam");
#print Dumper @samples;
my $count = 0;
foreach my $s (@samples) {
	my @tmp = split(/\//, $s);
	my $sortedBam = $tmp[$#tmp];
	$sortedBam =~ s/accepted_hits/sorted/;
	my $cmd1 = "java -Xmx4g -Djava.io.tmpdir=$outDir1/tmp -jar /opt/bioinf/picard/1.77/SortSam.jar SO=coordinate INPUT=$s OUTPUT=$outDir1/$sortedBam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &";
	print CMDFILE "$cmd1\n";
	$count++;
	print CMDFILE "p" . $count . "=\$!\n";	
	if (($count*3) > $nProc || $count == scalar (@samples)) {
		print CMDFILE "wait";
		for my $c (1..$count) { print CMDFILE " \$p" . $c};
		print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
		$count = 0;
	}
	print CMDFILE "\n";
}
if ($createBed) {
	foreach my $s (@samples) {
		my @tmp = split(/\//, $s);
		my $bed = $tmp[$#tmp];
		$bed =~ s/accepted_hits\.bam/bed/;
		my $cmd2 = "bamToBed -i $s > $outDir2/$bed &";
		print CMDFILE "$cmd2\n";
		$count++;
		print CMDFILE "p" . $count . "=\$!\n";	
		if ($count == $nProc || $count == scalar (@samples)) {
			print CMDFILE "wait";
			for my $c (1..$count) { print CMDFILE " \$p" . $c};
			print CMDFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
			$count = 0;
		}
	}
	print CMDFILE "\n";
}
my $cmd = "perl " . $libDir . "/send_email.pl" . " \"Sort_index_bams run completed on \" " . "NA" . " $runErrorLog $start";
print CMDFILE $cmd . "\n";
print CMDFILE "status=\$? 
if [ \$status != 0 ] ; then error_log \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n";
close CMDFILE;
#print "nohup sh $cmdFile.sh 2>&1 > ${cmdFile}.log &\n";
system("nohup sh ${cmdFile}.sh 2>&1 > ${cmdFile}.log &");

