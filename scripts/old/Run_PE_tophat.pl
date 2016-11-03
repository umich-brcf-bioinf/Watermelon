#! /usr/bin/perl

### Manjusha Pande, mpande@umich.edu, 05/23/15
### Read inner mate distance for parameter -r and write tophat command for PE samples
### Usage: $ perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/Run_PE_tophat.pl <projectDir> <gtfFile> <genome> <refGenomeBTindex> <novelTranscripts> <libraryType> <RABA> <nProc> <Date> <libDir> <max_intron_length> <option> <start> <error_log> <runInfoDir> <collaborator> <lite> [cmdFile]
### Example: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/Run_PE_tophat.pl /ccmb/BioinfCore/Projects/Manjusha/Rsq_test_v3/PEnovelTranscripts /ccmb/BioinfCore/Projects/Manjusha/Genomes/iGenomes/Caenorhabditis_elegans_UCSC_ce10/UCSC/ce10/Annotation/Genes/genes.gtf ce10 /ccmb/BioinfCore/Projects/Manjusha/Genomes/iGenomes/Caenorhabditis_elegans_UCSC_ce10/UCSC/ce10/Sequence/Bowtie2Index/genome 1 fr-unstranded 15 2014_05_15 /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts 3 1400182458 /ccmb/BioinfCore/Projects/Manjusha/Rsq_test_v3/PEnovelTranscripts/runInfo/Run_2014-05-15_15-34-18/Run_2014-05-15_15-34-18_error.log /ccmb/BioinfCore/Projects/Manjusha/Rsq_test_v3/PEnovelTranscripts/runInfo/Run_2014-05-15_15-34-18 0 [/ccmb/BioinfCore/Projects/Spence_959_spencejr_mceachin_RS4_all_compare/runInfo/Run_RSQ1_v3_part2_2014_05_14]
######################################################################################################################################################
use strict;
use warnings;
use Data::Dumper;

#=pod
my $projectDir = $ARGV[0];
my $gtfFile = $ARGV[1];
my $genome = $ARGV[2];
my $refGenomeBTindex = $ARGV[3];
my $novelTranscripts = $ARGV[4];
my $libraryType = $ARGV[5];
my $RABA = $ARGV[6];
my $nProc = $ARGV[7];
my $date = $ARGV[8];
my $libDir = $ARGV[9];
my $max_intron_length = $ARGV[10];
my $option = $ARGV[11];
my $start = $ARGV[12];
my $runErrorLog = $ARGV[13];
my $runInfoDir = $ARGV[14];	
my $collaborator = $ARGV[15];
my $lite = $ARGV[16];
#=cut

#my $projectDir = $ARGV[0];
#my $nProc = $ARGV[1];
#my $gtfFile = "/ccmb/BioinfCore/Projects/Manjusha/Genomes/iGenomes/Caenorhabditis_elegans_UCSC_ce10/UCSC/ce10/Annotation/Genes/#genes.gtf";
#my $genome = "ce10";
#my $refGenomeBTindex = "/ccmb/BioinfCore/Projects/Manjusha/Genomes/iGenomes/Caenorhabditis_elegans_UCSC_ce10/UCSC/ce10/Sequence/Bowtie2Index/genome";
#my $novelTranscripts = 1;
#my $libraryType = "fr-unstranded";
#my $date = "2014_01_28";
#my $libDir = "/ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts"; ### scripts dir
 
my $readsDir = "$ARGV[0]/reads";
my $alnDir = "$ARGV[0]/tophat";
my $qcDir = "$ARGV[0]/fastqc";
my $countDir = "$ARGV[0]/HTSeq";
my @tmp1 = split("/", $projectDir);
my $projectName = $tmp1[$#tmp1];

my $ss = "no"; ## unstranded library prep
if ($libraryType eq "fr-firststrand") {
	#$ss = "yes"; ## ??? not 100% sure
	$ss = "reverse";
} 
elsif ($libraryType eq "fr-secondstrand"){
	#$ss = "reverse";
	$ss = "yes";
} 

my $PET;
if ($ARGV[17]) {
	$PET = $ARGV[17];
}
else {
	$PET = $runInfoDir . "/Run_PE_tophat_" . $date;
}

open (OUTFILE, ">" . $PET . ".sh");
my $PETlog = "$PET.log";
my $runInfoFile =  glob($runInfoDir . "/Run*runinfo.txt");

my @samples = ();
my %samples = ();

open (INFILE1, "<$runInfoDir/Sample_list_" . $date . ".txt") || die "Can not open $runInfoDir/Sample_list_" . $date . ".txt for reading. $!\n";

while (<INFILE1>) {
	chomp $_;
	push(@samples, $_);
}
close INFILE1;
#=pod
open (INFILE2, "<$runInfoDir/innerMateDist.txt") || die "Can not open $runInfoDir/innerMateDist.txt for reading $!\n";
my @header = ();
	my $meanIDCol;
	my $sdIDCol;
while (<INFILE2>) {
	chomp $_;
	if ($_ =~ /^#/) {
		@header = split (/\t/, $_);
		for my $i (0..$#header) {
			if ($header[$i] =~ /Mate Inner Distance/) {
				$meanIDCol = $i;
			}
			elsif ($header[$i] =~ /Mate SD/) {
				$sdIDCol = $i;
			}
		}
	}
	else {
		my @tmp = split (/\t/, $_);
		foreach my $s (@samples) {
			if ($s eq $tmp[0]) {
				push (@{$samples{$s}}, @tmp[$meanIDCol,$sdIDCol]);
			}
		}
	}
	#print "\$meanIDCol is $meanIDCol and \$sdIDCol is $sdIDCol\n";
}
close INFILE2;
#=cut
#print Dumper %samples;

my $cmd;
print OUTFILE "#!/bin/sh\n";
print OUTFILE "\nPROGNAME=\$(basename \$0)
#echo \$PROGNAME
function error_exit
{
	echo \"\${PROGNAME}: \${1:-\"Unknown Error\"}\" >> $runErrorLog 2>&1
        ERROR=\"\${1:-\"Unknown Error\"} in \${PROGNAME}\" 
	DATE=\$(date +\"%a, %B %e, %Y\")
	TIME=\$(date +\"%I:%M %p\")
	perl ${libDir}/send_email.pl \"Your RNA-seq run ($projectName) terminated (exit status = 1) on \${DATE} at \${TIME}. 
\${ERROR}.\" \"\" $runErrorLog
	exit 1
}
function error_log
{
	echo \"\${PROGNAME}: \${1:-\"Unknown Error\"}\" >> $runErrorLog 2>&1		
}\n";

print OUTFILE "\necho \"Performing tophat alignment ...\";\n";
print OUTFILE "printf \"\nBowtie and TopHat versions:\n\" >> $runInfoFile\n";
print OUTFILE "bowtie2 --version | grep bowtie2 >> $runInfoFile\ntophat -v >> $runInfoFile\n";
print OUTFILE "printf \"\n\" >> $runInfoFile\n\n";	
if ($max_intron_length != 500000) {
	print OUTFILE "printf \"\nMaximum intron length re-set to $max_intron_length\n\" >> $runInfoFile\n";
}

my @tmpSamples1 = @samples;	
my $count = 1;
my $proc_count = $nProc;
my $sample_number = 0;
my $n_bg_processes = 0;

while (@tmpSamples1){
	my $sample_count = scalar(@tmpSamples1);
	my $p = int($proc_count/$sample_count);
	#print "\$sample_count is $sample_count, \$proc_count is $proc_count and \$p is $p \n";

	my $sample_name = shift(@tmpSamples1);
	$sample_number++;
	chomp $sample_name;
	my $tmpInFile1;
	my $tmpInFile2;	
	$tmpInFile1 = $readsDir . "/". $sample_name . "_R1.fastq.gz";
	$tmpInFile2 = $readsDir . "/". $sample_name . "_R2.fastq.gz";
	## for reads files that exist in the $readsDir (e. g. previously concatenated files, trimmed files) 
	my @inFiles = glob($readsDir."/*_R*.fastq*.gz");
	#print Dumper @inFiles;
	foreach my $f (@inFiles) {
		my @tmp = split(/\//, $f);
		my $file = $tmp[$#tmp];
		my $fileSample = $file;
		$fileSample =~ s/_R\d.fastq.*.gz//;
		if ($fileSample eq $sample_name) {
			#print "\$file is $file\n";
			if ($file =~ /R1/) {
				$tmpInFile1 = $f;
			} elsif ($file =~ /R2/) {
				$tmpInFile2 = $f;
			}
		}
	}

	my $tmpOutdir = $alnDir . "/" . $sample_name;
	my $innerMateDist = int(${$samples{$sample_name}}[0]);
	my $innerMateDistSd = int(${$samples{$sample_name}}[1]);
	
	if (!$novelTranscripts) { ## no novelTranscripts
		my $transcriptomeDir = $alnDir . "/transcriptome_data";
		if ($sample_number == 1){ ## align first sample 
			if (!$RABA) {	
				$cmd = "tophat -p $nProc --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd -G $gtfFile --transcriptome-index=$transcriptomeDir/${genome}_transcriptome -T --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &";
			} else { $cmd = "tophat -p $nProc --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd -G $gtfFile --transcriptome-index=$transcriptomeDir/${genome}_transcriptome --no-novel-juncs --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &"; }
			print OUTFILE $cmd . "\n";
			$n_bg_processes++;
			print OUTFILE "p" . $n_bg_processes . "=\$!\n";
			#print OUTFILE "\nwait;\n\n";
			print OUTFILE "wait";
			for my $bg (1..$n_bg_processes) { print OUTFILE " \$p" . $bg};
			print OUTFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
			$n_bg_processes = 0;
		} else { ## samples #2 onwards	
			if ($p >= 1){ ## more than one processor per sample
				if (!$RABA) {
					$cmd = "tophat -p $p --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd -G $gtfFile --transcriptome-index=$transcriptomeDir/${genome}_transcriptome -T --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &";
				} else { $cmd = "tophat -p $p --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd -G $gtfFile --transcriptome-index=$transcriptomeDir/${genome}_transcriptome --no-novel-juncs --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &"; }
				print OUTFILE $cmd . "\n";
				$n_bg_processes++;
				print OUTFILE "p" . $n_bg_processes . "=\$!\n";
				$count++;
				$proc_count = $proc_count-$p;	
			} else { 
				if($count <= $nProc) {
					#print "Assigning 1 processor\n";
					if (!$RABA) {
						$cmd = "tophat -p 1 --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd -G $gtfFile --transcriptome-index=$transcriptomeDir/${genome}_transcriptome -T --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &";
					} else { $cmd = "tophat -p 1 --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd -G $gtfFile --transcriptome-index=$transcriptomeDir/${genome}_transcriptome --no-novel-juncs --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &"; }					
					print OUTFILE $cmd . "\n";
					$n_bg_processes++;
					print OUTFILE "p" . $n_bg_processes . "=\$!\n";
					$count++;
					$proc_count = $proc_count-1;
				}
				if ($proc_count == 0) {
					#print OUTFILE "\nwait;\n";
					print OUTFILE "wait";
					for my $bg (1..$n_bg_processes) { print OUTFILE " \$p" . $bg};
					print OUTFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
					$n_bg_processes = 0;
					#print "wait!\n";
					$proc_count = $nProc;
					$count = 1;
				}
			}
				print OUTFILE "\n";
		}
	} else { ## novelTranscripts
		if ($RABA) {
			my $transcriptomeDir = $alnDir . "/transcriptome_data";
			if ($sample_number == 1){ ## align first sample 	
				$cmd = "tophat -p $nProc --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd -G $gtfFile --transcriptome-index=$transcriptomeDir/${genome}_transcriptome --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &";
				print OUTFILE $cmd . "\n";
				$n_bg_processes++;
				print OUTFILE "p" . $n_bg_processes . "=\$!\n";
				#print OUTFILE "\nwait;\n\n";
				print OUTFILE "wait";
				for my $bg (1..$n_bg_processes) { print OUTFILE " \$p" . $bg};
				print OUTFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
				$n_bg_processes = 0;
			} else { ## samples #2 onwards	
				if ($p >= 1){ ## more than one processor per sample
					$cmd = "tophat -p $p --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd -G $gtfFile --transcriptome-index=$transcriptomeDir/${genome}_transcriptome --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &";
					print OUTFILE $cmd . "\n";
					$n_bg_processes++;
					print OUTFILE "p" . $n_bg_processes . "=\$!\n";
					$count++;
					$proc_count = $proc_count-$p;	
				} else { 
					if($count <= $nProc) {
						#print "Assigning 1 processor\n";
						$cmd = "-I $max_intron_length  -p 1 --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd -G $gtfFile --transcriptome-index=$transcriptomeDir/${genome}_transcriptome --library-type $libraryType -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &";					
						print OUTFILE $cmd . "\n";
						$n_bg_processes++;
						print OUTFILE "p" . $n_bg_processes . "=\$!\n";
						$count++;
						$proc_count = $proc_count-1;
					}
					if ($proc_count == 0) {
						#print OUTFILE "\nwait;\n";
						print OUTFILE "wait";
						for my $bg (1..$n_bg_processes) { print OUTFILE " \$p" . $bg};
						print OUTFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
						$n_bg_processes = 0;
						#print "wait!\n";
						$proc_count = $nProc;
						$count = 1;
					}
				}
					print OUTFILE "\n";
			}
		} ## end if $RABA
		else { ## no $RABA
			if ($p >= 1){
				$cmd = "tophat -p $p --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &";
					print OUTFILE $cmd . "\n";
					$n_bg_processes++;
					print OUTFILE "p" . $n_bg_processes . "=\$!\n";
					$count++;
					$proc_count = $proc_count-$p;	
			} else {
				if($count <= $nProc) {
					#print "Assigning 1 processor\n";
					$cmd = "tophat -p 1 --b2-very-sensitive --no-coverage-search -r $innerMateDist --mate-std-dev $innerMateDistSd --library-type $libraryType -I $max_intron_length -o $tmpOutdir $refGenomeBTindex $tmpInFile1 $tmpInFile2 &";
					print OUTFILE $cmd . "\n";
					$n_bg_processes++;
					print OUTFILE "p" . $n_bg_processes . "=\$!\n";
					$count++;
					$proc_count = $proc_count-1;
				}
				if ($proc_count == 0) {
					#print OUTFILE "\nwait;\n\n";
					print OUTFILE "wait";
					for my $bg (1..$n_bg_processes) { print OUTFILE " \$p" . $bg};
					print OUTFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
					$n_bg_processes = 0;
					#print "wait!\n";
					$proc_count = $nProc;
					$count = 1;
				}
			}
			#print OUTFILE "\n";
		} ## end $RABA
	} ## end novelTranscripts
} ## end while
## wait for the tophat alignments to finish
#print OUTFILE "\nwait;\n";
print OUTFILE "wait";
for my $bg (1..$n_bg_processes) { print OUTFILE " \$p" . $bg};
print OUTFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
$n_bg_processes = 0;

## rename the accepted.hits.bam files output by tophat
my @tmpSamples2 = @samples;
print OUTFILE "\necho \"Re-naming alignment files ...\";\n";
while (@tmpSamples2){
	my $sample_name = shift(@tmpSamples2); 
	my $newName = $sample_name . "_accepted_hits.bam";
	print OUTFILE "mv $alnDir/$sample_name/accepted_hits.bam $alnDir/$sample_name/$newName;\n";
	print OUTFILE "status=\$? 
if [ \$status != 0 ] ; then error_exit \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n";
}

my $runningQC = 0;
#print "Running fastqc on aligned reads ...\n";
print OUTFILE "\necho \"Running fastqc on aligned reads ...\";\n";
foreach my $s (@samples) {
	my $alnFile = "$alnDir/$s/${s}_accepted_hits.bam";
	my $qcFile = "$qcDir/${s}_accepted_hits_fastqc.zip";
	#print "alnfile is $alnFile\nqcFile is $qcFile\n";
	#unless (-e $qcFile) {
		$cmd = "fastqc $alnFile -o $qcDir -t $nProc -q";
		print OUTFILE $cmd . ";\n";
		print OUTFILE "status=\$? 
if [ \$status != 0 ] ; then error_log \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n";
		if (!$runningQC) {
			$runningQC = 1; }
	#}
	#else { #print "Skipping $alnFile, fastqc report for this file exist/s in fastqc dir.\n";
	#	print OUTFILE "#Skipping $alnFile, fastqc report for this file exist/s in fastqc dir\n"; }
}		
#$cmd = "fastqc $alnDir/Sample*/*accepted_hits.bam -o $qcDir -t $nProc -q";
#print OUTFILE $cmd . ";\n";
open (RUNINF, "<$runInfoFile") || die "Can not open $runInfoFile for reading $!\n";
while (<RUNINF>) {
	if ($_ =~ /FastQC version/) {
		$runningQC = 0;
	}
}
if ($runningQC) {
	print OUTFILE "printf \"\nFastQC version:\n\" >> $runInfoFile\n";
	print OUTFILE "fastqc -v >> $runInfoFile\n\n";
}
else { print "\n*** Skipping fastqc on aligned reads, QC reports for samples exist in fastqc dir.***\n"; }

### Get QC and alignment metrics
print OUTFILE "\necho \"Get QC and alignment metrics ...\";\n";
$cmd = "perl $libDir/getQCMetrics.pl $alnDir $runInfoDir/Sample_list_" . $date . ".txt $runInfoDir $runInfoFile";
print OUTFILE $cmd . ";\n";
print OUTFILE "status=\$? 
if [ \$status != 0 ] ; then error_log \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n";

if(!$lite) {
	## Getting HTSeq counts from aligned reads ...

	print OUTFILE "\necho \"Sorting tophat output bam files by read names ...\";\n";	
	$count = 1;
	$proc_count = $nProc;
	$n_bg_processes = 0;
	
	foreach my $s (@samples) {
		my $bamFile = "$alnDir/$s/${s}_accepted_hits.bam";
		my $sortFile = "$countDir/${s}_namesorted";
		$cmd = "samtools sort -n $bamFile $sortFile &";
		print OUTFILE $cmd . "\n";
		$n_bg_processes++;
		print OUTFILE "p" . $n_bg_processes . "=\$!\n";
		$count++;
		$proc_count--;
		if ($proc_count == 0 || $count > scalar @samples) {
			print OUTFILE "wait";
			for my $bg (1..$n_bg_processes) { print OUTFILE " \$p" . $bg};
			print OUTFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
			$n_bg_processes = 0;
			$proc_count = $nProc;
		}	
	}
	print OUTFILE "\necho \"Getting HTSeq counts from aligned reads...\";\n";
	$count = 1;
	$proc_count = $nProc;
	$n_bg_processes = 0;
	foreach my $s (@samples) {
		my $sortedFile = "$countDir/${s}_namesorted.bam";
		my $countFile = "$countDir/${s}_counts.txt";
		$cmd = "python -m HTSeq.scripts.count -f bam -s $ss -m intersection-nonempty -q $sortedFile $gtfFile > $countFile &";
		print OUTFILE $cmd . "\n";
		$n_bg_processes++;
		print OUTFILE "p" . $n_bg_processes . "=\$!\n";
		$count++;
		$proc_count--;
		if ($proc_count == 0 || $count > scalar @samples) {
			print OUTFILE "wait";
			for my $bg (1..$n_bg_processes) { print OUTFILE " \$p" . $bg};
			print OUTFILE "; status=\$?\nif [ \$status != 0 ] ; then error_exit \"An error occurred near line \$((\${LINENO}-1))\" ; fi\n\n";
			$n_bg_processes = 0;
			$proc_count = $nProc;
		}	
	}	
	$cmd = "$libDir/mergeHTSeqCountFiles.pl $countDir";
			print OUTFILE $cmd . "\n";
			print OUTFILE "status=\$? 
	if [ \$status != 0 ] ; then error_exit \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n\n";
}
if ($option == 1) {
	print OUTFILE "\nprintf \"\nEnd time:\n\" >> $runInfoFile\n";
	print OUTFILE "date >> $runInfoFile\n";
	$cmd = "perl " . $libDir . "/tab2xls.pl " . $runInfoDir . " " . $runInfoDir . " " . $runInfoFile;
	print OUTFILE $cmd . "\n";
	$cmd = "perl " . $libDir . "/send_email.pl" . " \"Your RNA-seq run ($projectName) completed on \" " . "$collaborator $runErrorLog $start";
	print OUTFILE $cmd . "\n";
	print OUTFILE "status=\$? 
if [ \$status != 0 ] ; then error_log \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n";
}
close OUTFILE;
system("chmod +x $PET.sh");
#system("nohup sh $PET.sh 2>&1 > $PETlog &");





 
