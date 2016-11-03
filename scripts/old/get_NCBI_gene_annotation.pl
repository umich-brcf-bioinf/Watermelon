#!/usr/bin/perl

use strict;
use warnings;
use  Data::Dumper;
use POSIX;
use Cwd;
use File::Path;
use Archive::Extract;
use constant DATE => strftime("%Y_%m_%d", localtime);
umask 002; ## set default directory permissions set to 775 and default file permissions to 664

### Manjusha Pande, mpande@umich.edu, 03/18/13, modified 11/14/13
### 1. Mapping gene names from Cuffdiff output to refseq ids using the gtf file
### 2. Mapping refseq ids to NCBI Entrez Gene ids using NCBI gene2ref
### 3. Getting description for each gene/transcript and GO annotation for DE genes using NCBI gene_info and gene2GO
###### Usage:
### $ perl get_NCBI_gene_annotation.pl <path/to/inDir> <path/to/outDir> <genome> <path/to/gtfFile> <annotDir> [runInfoFile]
###### Example
### perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/get_NCBI_gene_annotation.pl /ccmb/home/mpande/RNA-seq_pipeline/Pipeline_v2/Sample_cuffdiff_output /ccmb/home/mpande/RNA-seq_pipeline/Pipeline_v2/Sample_cuffdiff_output/report mm10 /ccmb/BioinfCore/Reference_Genomes/programs/tophat/mouse/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf
### perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/get_NCBI_gene_annotation.pl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/AnnotateCuffdiffOutput/ENSEMBL_test /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/AnnotateCuffdiffOutput/ENSEMBL_test sacCer3 /ccmb/BioinfCore/Projects/Manjusha/Genomes/iGenomes/Saccharomyces_cerevisiae_Ensembl_R64-1-1/Ensembl/R64-1-1/Annotation/Archives/archive-2014-05-27-15-52-46/Genes/genes.gtf
### perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/get_NCBI_gene_annotation.pl /ccmb/BioinfCore/Projects/Tang_xtang2_RS2_mceachin_Run_1298/cuffdiff/G166_v_WT/NCBI_annotation_test /ccmb/BioinfCore/Projects/Tang_xtang2_RS2_mceachin_Run_1298/cuffdiff/G166_v_WT/NCBI_annotation_test Oryza_sativa /ccmb/BioinfCore/Projects/McEachin/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Annotation/Archives/archive-2014-05-28-13-55-57/Genes/genes.gtf /ccmb/BioinfCore/Projects/Tang_xtang2_RS2_mceachin_Run_1298/cuffdiff/G166_v_WT/NCBI_annotation_test/NCBI_annotation_2015_06_12 /ccmb/BioinfCore/Projects/Tang_xtang2_RS2_mceachin_Run_1298/cuffdiff/G166_v_WT/NCBI_annotation_test/test_runinfo.txt

my $inDir = $ARGV[0];
my $runInfoDir = $ARGV[1];
my $genome = $ARGV[2];
my $gtfFile = $ARGV[3];
my $annotDir = $ARGV[4];

#unless (-e $runInfoDir || mkdir ($runInfoDir, 0775)) {
#	die "Unable to create $runInfoDir \n";
#	};

#my $gtfFile = "/ccmb/BioinfCore/Reference_Genomes/programs/tophat/mouse/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf";
my $outFile;
my %geneId2geneInfo = (); ### map gene id to gene name and descritpion using NCBI gene info
my %refseq2geneInfo = ();  ### map refseq id to gene id, gen name and descritpion using NCBI gene2ref
my %geneId2GO = (); ### map gene id to GO annotation using NCBI gene2GO
my %geneName2refseq = (); ### map gene name in the cuffdiff output to transcript id using the isoform file
my @header;
my $taxId; # NCBI taxonomy Id, Human = 9606, Mouse = 10090, Worm = 6239, zebrafish = 7955
my $gene_count = 0;
my $mapped_gene_count = 0;
my $unmapped_DEGcount = 0;
my $isoform_count = 0;
#my $unannotatedGenome = 0;
#my $geneCountFile = $runInfoDir . "/DE_counts.txt";
my $geneCountFile;
if ($ARGV[5]) {
	$geneCountFile =  $ARGV[5];
} else {
 	$geneCountFile = glob($runInfoDir . "/*_runinfo.txt");
}
if ($genome =~ /^hg/ || $genome =~ /^GRCh/) { ## human
	$taxId = 9606; 
} elsif ($genome =~ /^mm/ || $genome =~ /^GRCm/) { ## mouse
	$taxId = 10090; 
} elsif ($genome =~ /^rn/) { ## rat
	$taxId = 10116; 
} elsif ($genome =~ /^OA/) { ## sheep
	$taxId = 9940;
} elsif ($genome =~ /^ce/ || $genome =~ /^WS/ || $genome =~ /^WBcel/) { ## c. elegans
	$taxId = 6239; 
} elsif ($genome =~ /^Zv/ || $genome =~ /^danRer/ || $genome =~ /^GRCz/) { ## zebrafish
	$taxId = 7955;
} elsif ($genome =~ /^R/ || $genome =~ /^sacCer/) { ## saccharomyces cerevisiae S288c
	$taxId = 559292;
} elsif ($genome =~ /^IRGSP/ || $genome =~ /^Oryza/) { ## Oryza sativa japonica
	$taxId = 39947;
} elsif ($genome =~ /^TAIR/) { ## Arabidopsis thaliana
	$taxId = 3702;
} elsif ($genome =~ /^Ecoli_K12/) { ## Escherichia coli O17 str. K12 substr. MG1655
	$taxId = 511145;
} elsif ($genome =~ /^Ecoli_LF82/) { ## Escherichia coli LF82
	$taxId = 591946;
#} elsif ($genome =~ /^PG_ATCC_33277/) { ## Porphyromonas gingivalis ATCC 33277 
	#$taxId = 431947;
#} elsif ($genome =~ /^AB_AB0057/) { ## Acinetobacter baumannii AB0057
	#$taxId = 480119;
} else {  
	#$unannotatedGenome = 1;
	die "Unsupported genome!\n";
}

######################################################################################################################################
my $lastUpdate = "2015_06_11";
my $gene2RefseqFile = "/nfs/med-bfx-common/pipelines/RNASeq/references/Annotation/${lastUpdate}/NCBI_gene2refseq_" . $lastUpdate;
my $geneInfoFile = "/nfs/med-bfx-common/pipelines/RNASeq/references/Annotation/${lastUpdate}/NCBI_gene_info_" . $lastUpdate;
my $gene2GOFile = "/nfs/med-bfx-common/pipelines/RNASeq/references/Annotation/${lastUpdate}/NCBI_gene2go_" . $lastUpdate;
#print "\$gene2RefseqFile is $gene2RefseqFile \n";

## update NCBI annotation
#my $annotDir = "/ccmb/BioinfCore/Projects/NCBI_annotation/" . DATE;
#unless (-e $annotDir || mkpath $annotDir) {
#	die "Unable to create $annotDir \n";
#}
my $cwd = cwd(); ## current working dir
chdir($annotDir) or die "Cannot change dir to $annotDir $!";
my $cmd;
my @NCBI_resources = ("gene2refseq", "gene_info", "gene2go");
my $updatedAnnot = 0;

foreach my $file (@NCBI_resources) {
	my $resource = "$annotDir/NCBI_${file}_" . DATE;
	unless (-e $resource) {
		print "Downloading NCBI $file...\n"; 
		$cmd = "wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/${file}.gz";
		#print "\$cmd is $cmd\n";
		system ($cmd);
		my $zipFile = "${file}.gz";
		print "Extracting $zipFile.....doesn't it?\n";
		my $archive = Archive::Extract->new(archive => $zipFile);
		my $ok = $archive->extract or die $archive->error;
		rename ("$annotDir/$file", "$annotDir/NCBI_${file}_" . DATE);
		$updatedAnnot = 1;
	} #else { print "$resource exists!\n"; }	 
}
if ($updatedAnnot) {print "Downloaded latest annotation from NCBI ftp site on " . DATE . "\n";}
foreach my $file (@NCBI_resources) {
	my $resource = "$annotDir/NCBI_${file}_" . DATE;
	if ($file =~ /gene2refseq/ ) {
		if (-e $resource && $gene2RefseqFile ne $resource) { 
			$gene2RefseqFile = $resource;
			print "NCBI $file updated on " . DATE . "\n";
		} else { print "Using NCBI $file updated on " . $lastUpdate . "\n"; }
	}
	elsif ($file =~ /gene_info/) {
		if ( -e $resource && $geneInfoFile ne $resource) { 
			$geneInfoFile = $resource;
			print "NCBI $file updated on " . DATE . "\n";
		} else { print "Using NCBI $file updated on " . $lastUpdate . "\n"; }
	}
	elsif ($file =~ /gene2go/ ) {
		if ( -e $resource && $gene2GOFile ne $resource) { 
			$gene2GOFile = $resource;
			print "NCBI $file updated on " . DATE . "\n";
		} else { print "Using NCBI $file updated on " . $lastUpdate . "\n"; }
	}
}
unlink glob "$annotDir/*.gz";	
chdir ($cwd) or die "Cannot change dir to $cwd $!"; ## back in $cwd

#######################################################################################################################################
#=pod
print "Getting gene info...\n";
### Step1: read NCBI_gene_info file, map gene Id to gene name and description for the particular tax and store in %geneId2geneInfo
open (INFILE1, "<".$geneInfoFile) or die "Cannot open $geneInfoFile $!\n";
while (<INFILE1>) {
	my $line = $_;
	if ($line=~/^$taxId/) {
		chomp($line);
		my @tmp = split(/\t/,$line);
		push (@{$geneId2geneInfo{$tmp[1]}}, @tmp[2,8]);  ### gene Id to gene name and description
	}
}
close INFILE1;
#print Dumper %geneId2geneInfo;
print "Total genes found for $genome in $geneInfoFile: " . scalar (keys (%geneId2geneInfo)) . "\n";
#=pod

print "Mapping gene info to refSeq ids...\n";
### Step2: read NCBI_gene2refseq file, map refSeq id to gene Id, gene name and description for the particular tax and store in %refseq2geneInfo
open (INFILE2, "<".$gene2RefseqFile) or die "Cannot open $gene2RefseqFile $!\n";
while (<INFILE2>) {
	my $line = $_;
	if ($line=~/^$taxId/) {
		chomp($line);
		my @tmp = split(/\t/,$line);
		$tmp[3] =~ s/\.(\d+)$//; ### remove the version number in the end
		#print "$tmp[3]\n";
		if (!exists ($refseq2geneInfo{$tmp[3]})) {
			push (@{$refseq2geneInfo{$tmp[3]}}, $tmp[1]);  ### RefSeq Id to gene Id
			if (exists ($geneId2geneInfo{$tmp[1]})) {
				push (@{$refseq2geneInfo{$tmp[3]}}, @{$geneId2geneInfo{$tmp[1]}}); ### RefSeq Id to gene Id, gene name and descritption
			}
		}
	}
}
close INFILE2;
#print Dumper %refseq2geneInfo;
print "Total refSeq ids mapped to genes for $genome in $gene2RefseqFile: " . scalar (keys (%refseq2geneInfo)) . "\n";

my $GO_category;
my $GO_entry;
print "Getting GO annotation...\n";
### Step3: read NCBI_gene2go file, map gene Ids to GO terms for the particular tax and store in %geneId2GO
open (INFILE3, "<".$gene2GOFile) or die "Cannot open $gene2GOFile $!\n";
while (<INFILE3>) {
	my $line = $_;
	if ($line=~/^$taxId/) {
		chomp($line);
		my @tmp = split(/\t/,$line);
		if ($tmp[$#tmp] eq "Component") {
			$GO_category = "GO_CC"; }
		elsif ($tmp[$#tmp] eq "Process") {
			$GO_category = "GO_BP"; }
		elsif ($tmp[$#tmp] eq "Function") {
			$GO_category = "GO_MF"; }
		$GO_entry = $GO_category . " (" . $tmp[2] . ")~ " . $tmp[5];
		my $found_entry = 0;
		if (exists ($geneId2GO{$tmp[1]})) {
			foreach my $t (@{$geneId2GO{$tmp[1]}}) {
				if ($t eq $GO_entry) {
					$found_entry = 1;
				}
			}
		}
		if (!$found_entry) {
			push (@{$geneId2GO{$tmp[1]}}, $GO_entry); ### gene Id to GO
		}
	}
}
close INFILE3;
print "Total gene ids mapped to GO terms for $genome in $gene2GOFile: " . scalar (keys (%geneId2GO)) . "\n";
#print Dumper %geneId2GO;
#=cut

print "Mapping gene names to transcript ids...\n";
### Step4: read gtf file, map gene names to transcript ids and store in %geneName2refseq, for ENSEMBL gtf, transcript ids are ENSEMBL transcript ids, not refSeq ids 
open (INFILE4, "<".$gtfFile) or die "Cannot open $gtfFile $!\n";
while (<INFILE4>) {
	my $line = $_;
	chomp($line);
	my @tmp1 = split(/\t/,$line);
	my @tmp2 = split(/;/,$tmp1[$#tmp1]);
	#print Dumper @tmp2;
	my $geneName = "";
	my $refseqId = ""; ## included in UCSC gtf, not ENSEMBL gtf
	for my $i (0..$#tmp2) {
		if ($tmp2[$i] =~ /gene_name/ )
		{ 
			$tmp2[$i] =~ s/gene_name|"|\s//g;
			$geneName = $tmp2[$i];
		}
		elsif ($tmp2[$i] =~ /transcript_id/ ) 
		{ 
			$tmp2[$i] =~ s/transcript_id|"|\s//g;
			$refseqId = $tmp2[$i]; 
		}
	}
	#print "gene name is $geneName\trefseq id is $refseqId\n";
	push (@{$geneName2refseq{$geneName}}, $refseqId);
}
close INFILE4;
print "Total gene names mapped to transcript ids in the gtf file: " . scalar (keys (%geneName2refseq)) . "\n";
#print Dumper %geneName2refseq;
	
#####################################################################################################################################
#=pod
my @tmp = split("/", $inDir);
open (COUNTINFILE, "<".$geneCountFile) or die "Can not read $geneCountFile $!";
my $DECountHeader = 0;
while (<COUNTINFILE>) 
{
	if ($_ =~ /^Number of DE genes and isoforms/)
	{
		$DECountHeader = 1;
	}
}
close COUNTINFILE;
open (COUNTFILE, ">>".$geneCountFile) or die "Can not write to $geneCountFile $!";
if (!$DECountHeader)
{
	print COUNTFILE "\nNumber of DE genes and isoforms:\nDataset_name\t#DE_genes\t#DE_isoforms\t#DE genes mapped to NCBI Ids\n";
}
my $comp = $tmp[$#tmp];

my $DEGIdFile = $inDir . "/" . $comp . "_DEG_ids.txt";
open (DEGFILE1, ">".$DEGIdFile) or die "Can not write to $DEGIdFile $!";
print DEGFILE1 "#Gene_id\n";

my $DEGNameFile = $inDir . "/" . $comp . "_DEG_names.txt";
open (DEGFILE2, ">".$DEGNameFile) or die "Can not write to $DEGNameFile $!";
print DEGFILE2 "#Gene_name\tFC\n";

### Step5: get annotation for genes and transcripts in the cuffdiff output files
### read formatted cuffdiff output files (for genes and isoforms), 
### map gene names in the output file to gene names in %geneName2refseq (gene name mapped to transcript id in the gtf file), 
### map corresponding refseq id (transcript id) to refseq ids in %refseq2geneInfo for more precise mapping, 
### if no match for transcript id in %refseq2geneInfo (ENSEMBL gtfs), map gene name from the cuffdiff output to gene names in %geneId2geneInfo

my @DEFiles = glob($inDir. "/*{gene,isoform}.foldchange*.txt");
#print Dumper @DEFiles;
print "Annotating CuffDiff output...\n";
for my $f (0..$#DEFiles) 
{ 
	#print "$f\n";
	unless ($DEFiles[$f] =~ /_annotated.txt$/) {
		#print "$f\n";
		$outFile = $DEFiles[$f];
		$outFile =~ s/.txt$/_annotated.txt/;
		open (OUTFILE, ">$outFile") || die "Can not open $outFile for writing $!\n";
		open (INFILE, "<".$DEFiles[$f]) or die "Cannot open $DEFiles[$f] $!\n";
		print "Annotating $DEFiles[$f]...\n";
		@header = ();
		my $header = 0;
		my $geneNameCol;
		my $DECol;
		my $FCCol;
		my $count = 0; 
		my @geneMatrix = ();
		my $genes;
		my $isoforms;
		my @unmapped = ();
		
		if ($DEFiles[$f] =~ /gene/) { 
			$genes = 1;
		} elsif ($DEFiles[$f] =~ /isoform/) { 
			$isoforms = 1;
		}

		while (<INFILE>) {
			my $refseqId = "";
			my $geneName = "";
			my $GO_annot = "";
			
			if (!$header) {
				if ($_ =~ /^test_id/ ) {
					chomp $_;
					my @tmp = split(/\t/, $_);
					push (@header, @tmp[0,2], "Entrez GeneId", "Description", @tmp[3..12,14..$#tmp]);	
					for my $i (0..$#tmp) {
						if ($tmp[$i] =~ /gene$/) {
							$geneNameCol = $i;
						}
						if ($tmp[$i] =~ /Diff Exp/) {
							$DECol = $i;
						}
						if ($tmp[$i] =~  /Fold Change/) {
							$FCCol = $i;
						}
					}
					$header = 1;
				}
				else { print OUTFILE $_; }
			}
			#print OUTFILE "\n";
			else {
				chomp $_;
				#print "Line is $_\n";
				my @tmp = split(/\t/, $_);
				$geneName = $tmp[$geneNameCol];
				$geneName =~ s/\s//g;
				my @tmp2 =  split(',', $geneName);
				$geneName = $tmp2[0];
				#print "Gene name is $geneName\n";
				my $FC = $tmp[$FCCol];
				my $foundId = 0;
				if (exists ($geneName2refseq {$geneName})) { 
					#print "Record exists for $geneName in geneName2refseq \n";
					foreach my $i (@{$geneName2refseq{$geneName}}) {
						#$refseqId = ${$geneName2refseq{$geneName}}[0];	
						$refseqId = $i;
						#print "refseqId before is $refseqId\n";
						$refseqId =~ s/\.(\d+)$//; ### remove the version number in the end (needed for NCBI .gtf)
						#print "refseqId after is $refseqId\n";
						#print "foundId is $foundId\n";	
						if (!$foundId && exists ($refseq2geneInfo{$refseqId})) {
							push (@{$geneMatrix[$count]}, @tmp[0,2], @{$refseq2geneInfo{$refseqId}}[0,2],  @tmp[3..12,14..$#tmp] );
							$foundId = 1;
						}
						#if (exists ($refseq2geneInfo{$refseqId})) {print "refseq2geneInfo is @{$refseq2geneInfo{$refseqId}}\n";}
						#else {"Print $refseqId does not exist in refseq2geneInfo\n";}
					}
				}
				if (!$foundId) { 
					#print "Could not map transcript id to gene Id..try mapping gene name.\n"; ## for ENSEMBL gtfs
					my $geneId = "";
					my $geneDesc = "";
					#my $geneName_ucFirstlc = ucfirst lc $geneName; ## this was needed for Oryza_sativa
					#my $geneName_ucLastlc = reverse(ucfirst(reverse(ucfirst(lc $geneName)))); ## this was needed for Oryza_sativa
					foreach my $g (keys %geneId2geneInfo) {
						#if (${$geneId2geneInfo{$g}}[0] && ($geneName eq ${$geneId2geneInfo{$g}}[0] || $geneName_ucFirstlc eq ${$geneId2geneInfo{$g}}[0] || $geneName_lc eq ${$geneId2geneInfo{$g}}[0])) {
						if (${$geneId2geneInfo{$g}}[0] && lc $geneName eq (lc ${$geneId2geneInfo{$g}}[0])) { ## case insensitive match
							#print "\$geneId is $g, \$geneName is ${$geneId2geneInfo{$g}}[0] and \$geneDesc is ${$geneId2geneInfo{$g}}[1]\n";
							$geneId = $g;
							$geneDesc = ${$geneId2geneInfo{$g}}[1];	
							$foundId = 1;
						}		
					}
 					#print "\$geneId is $geneId, \$geneName is $geneName and \$geneDesc is $geneDesc\n";
					#print "Adding $tmp[2]to the matrix..\n";
					push (@{$geneMatrix[$count]}, @tmp[0,2], $geneId, $geneDesc, @tmp[3..12,14..$#tmp]); 
				}
				if (! $foundId) { push (@unmapped, $geneName); }
				#print "geneMatrix[$count][2] is $geneMatrix[$count][2]\n";
				#print "DECol is $DECol and value is $tmp[$DECol]\n";

				if ( defined $tmp[$DECol] && $tmp[$DECol] eq "YES ") {
					if($isoforms) {				
						$isoform_count++;
					}
					#$tmp[0] =~ s/\s//g;
					#$tmp[2] =~ s/\s//g;
					#push (@{$geneName2refseq{$tmp[2]}}, $tmp[0]); ### map gene name to refseq id

					#$refseqId = $tmp[0];
					#$refseqId =~ s/\s//g;
					#if (exists ($refseq2geneInfo{$refseqId})) {
						#push (@{$geneMatrix[$count]}, @tmp[0,2], @{$refseq2geneInfo{$refseqId}}[0,2],  @tmp[3..12,14..$#tmp] );
					#}
					#else { push (@{$geneMatrix[$count]}, @tmp[0,2], "", "", @tmp[3..12,14..$#tmp]) }
					#}						
					elsif ($genes) {
						$gene_count++;			
						print DEGFILE2 "$geneName\t$FC\n";
						if ($geneMatrix[$count][1]) {
							#print "This is a DE gene: $geneMatrix[$count][2]\n";
							if ($geneMatrix[$count][2]) {
								print DEGFILE1 "$geneMatrix[$count][2]\n";
							} else { $unmapped_DEGcount++; }
						}
						if ( exists ($geneId2GO{$geneMatrix[$count][2]}) ) {
							my $term_count = scalar (@{$geneId2GO{$geneMatrix[$count][2]}});
							my $cnt = 1;			
							foreach my $g (@{$geneId2GO{$geneMatrix[$count][2]}}) {
								$GO_annot = $GO_annot . $g ;
								if ($cnt != $term_count) {
									$GO_annot = $GO_annot . "; " ;
								}
								$cnt++;
							}
						}
					}
					push (@{$geneMatrix[$count]}, $GO_annot);	
				}
				$count++;
			}
		}
		#print "$count\n";
		close INFILE;
		if ($genes && $gene_count) {
			push (@header, "GO Annotation");
		}
	#=pod
		foreach my $h (@header) {
			print OUTFILE "$h\t";
		}
		print OUTFILE "\n";
		#print Dumper @geneMatrix;
		#print scalar (@geneMatrix) . "\n";
		#my $r_count = 1;
		foreach my $r (@geneMatrix) {
			#my $c_count = 1;
			foreach my $c (@{$r}) {		
				#if (!defined($c)) {print "$r_count, $c_count\n";}
				print OUTFILE "$c\t";
				#$r_count++;
			}
			print OUTFILE "\n";
			#$r_count++;
		}
		close OUTFILE;
		 
		if ($genes) { 
			$mapped_gene_count = $gene_count - scalar (@unmapped);
			print "Out of total $count, " . scalar (@unmapped) . " genes did not map to NCBI gene Ids for $comp\n";
			print "Out of $gene_count DEGs, $unmapped_DEGcount genes did not map to NCBI gene Ids for $comp\n"; 
		}
		#print Dumper @unmapped;
		if (scalar @unmapped && $genes) {
			open (UNMAPPED, ">$inDir/${comp}_unmapped_ids.txt") or die "Can not write to $inDir/${comp}_unmapped_ids.txt $!";
			print UNMAPPED "Following genes did not map to NCBI gene ids:\n";
			foreach my $id (@unmapped) { print UNMAPPED "$id\n"; }
			close UNMAPPED;
		}			
	} ## end unless ($DEFiles[$f] =~ /_annotated.txt$/)
#=cut
} ## end for my $f (0..$#DEFiles) 
#print Dumper %geneName2refseq;	
close DEGFILE1;
close DEGFILE2;
my $mapped_ids = $gene_count - $unmapped_DEGcount;
print COUNTFILE "$comp\t$gene_count\t$isoform_count\t$mapped_ids\n";		
close COUNTFILE;
#=cut
print "All done!\n";
#<>;
