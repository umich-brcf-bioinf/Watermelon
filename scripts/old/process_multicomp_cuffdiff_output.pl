#! /usr/bin/perl 

### Manjusha Pande, mpande@umich.edu, 05/02/14
### Separate multi-comparison cuffdiff output into separate files
### Usage: $ perl process_multicomp_cuffdiff_output.pl <inDir> <libDir> <genome> <gtfFile> <annotDir> <runInfoDir> <reportDir> <runInfoFile> <runErrorLog> [cmdFile];

### Example: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/process_multicomp_cuffdiff_output.pl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/Parse_multicomp_output/T0_v_T1_v_T2_v_T3_v_T4_v_T5_v_T6_v_T7 /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts hg19 /ccmb/BioinfCore/Projects/Manjusha/Genomes/iGenomes/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/Parse_multicomp_output/T0_v_T1_v_T2_v_T3_v_T4_v_T5_v_T6_v_T7 /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/Parse_multicomp_output/T0_v_T1_v_T2_v_T3_v_T4_v_T5_v_T6_v_T7 /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/Parse_multicomp_output/T0_v_T1_v_T2_v_T3_v_T4_v_T5_v_T6_v_T7/Run_2014-05-02_21-51-15_runinfo.txt
#########################################################################################################################
use strict;
use warnings;
use Data::Dumper;

my $inDir = $ARGV[0];
my $libDir = $ARGV[1];
my $genome = $ARGV[2];
my $gtfFile = $ARGV[3];
my $annotDir = $ARGV[4];
my $runInfoDir = $ARGV[5];
my $reportDir = $ARGV[6];
my $runInfoFile = $ARGV[7];
my $runErrorLog =$ARGV[8];

my $cmdFile;
if ($ARGV[9]) {
	$cmdFile = $ARGV[9];
	#print "Cmd file is $cmdFile\n";
	open (PMC, ">>$cmdFile") || die "Cannot open $cmdFile for writing. $!\n";
}
else {
	my $libDir =  "/ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts";
	my $runErrorLog = "$inDir/PMC_error.log";
	$cmdFile = "$inDir/process_multi_comps.sh";
	open (PMC, ">$cmdFile") || die "Cannot open $cmdFile for writing. $!\n";
	print PMC "#!/bin/sh\n";
	print PMC "\nPROGNAME=\$(basename \$0)
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
#print "Cmd file is $cmdFile\n";
my %comps;

my @inDirPath = split(/\//, $inDir);
my $inDirName = $inDirPath[$#inDirPath];
my @inFiles = glob($inDir . "/" . $inDirName. "*{gene,isoform}*.txt");
#print Dumper @inFiles;
#=pod
foreach my $inFile (@inFiles) {
	open (INFILE, "<$inFile") || die "Cannot open $inFile for reading. $!\n";
	my @inFilePath = split(/\//, $inFile);
	my $inFileName = $inFilePath[$#inFilePath];
	$inFileName =~ m/(_[gene|isoform].*)$/;
	my $ext = $1;
	my $foundHeader = 0;
	my $header;
	my $notes;
	my $sample1_col;
	my $sample2_col;
	my $group1;
	my $group2;

	#print "extension is $ext\n";
	while (<INFILE>) {
		my $isHeader = 0;
		if ($_ =~ /^test_id/) {
			$header = $_;
			chomp $_;
			$foundHeader = 1;
			$isHeader = 1;
			my @header = split (/\t/, $_);
			for my $h (0..$#header){
				if ($header[$h] =~ /sample_1/) {
					$sample1_col = $h;
				}
				elsif ($header[$h] =~ /sample_2/) {
					$sample2_col = $h;
				}	
			}
		}
		if (!$foundHeader) {
			$notes = $notes . $_;
		}
		else {
			if (!$isHeader) {
				chomp $_;
				my @tmp = split (/\t/, $_);
				my $group1 = $tmp[$sample1_col];
				$group1 =~ s/(\s)+//;
				my $group2 = $tmp[$sample2_col];
				$group2 =~ s/(\s)+//;
				my $comp = $group1 . "_v_" . $group2;
				push (@{$comps{$comp}}, $_);
			}
		}
	}
	close INFILE;
	#print "Notes is $notes\nHeader is $header\nSample1 col is $sample1_col and Sample2 col is $sample2_col\n";
	#print Dumper keys(%comps);
	foreach my $k (keys (%comps)) {
		my $compDir = "$inDir/$k";
		unless (-e $compDir || mkdir ($compDir, 0775)) {
			die "Unable to create $compDir \n";
		}
		open (OUTFILE, ">$compDir/${k}$ext") || die "Cannot create $compDir/${k}$ext. $!.\n";
		print OUTFILE $notes, $header;
		my @v = @{$comps{$k}};
		foreach my $v (@v) {
			print OUTFILE "$v\n";
		}
		close OUTFILE;
	}
}
my @comps = keys(%comps);
@comps = sort {$a cmp $b} (@comps);
#print Dumper @comps;
foreach my $c (@comps) {
	print PMC "\necho \"Processing comparison $c ...\";\n";
	my $dir = "$inDir/$c";
	#print "dir is $dir\n";
	my $cmd;
	$cmd = "perl $libDir/get_NCBI_gene_annotation.pl $dir $runInfoDir $genome $gtfFile $annotDir $runInfoFile";
			print PMC $cmd . "\n";
			print PMC "status=\$? 
if [ \$status != 0 ] ; then error_log \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n";

			$cmd = "perl " . $libDir . "/get_DAVIDchartReport.pl -dir " . $dir . " -type ENTREZ_GENE_ID -thd 0.05";
			print PMC $cmd . "\n";
			#print PMC "status=\$? 
#if [ \$status != 0 ] ; then error_log \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n";

			$cmd = "perl " . $libDir . "/tab2xls.pl " . $dir . " " . $reportDir;
			print PMC $cmd . "\n";
			print PMC "status=\$? 
if [ \$status != 0 ] ; then error_exit \"An error occurred at line \$((\${LINENO}-2))\" ; fi\n";
}

close PMC;
#=cut	
