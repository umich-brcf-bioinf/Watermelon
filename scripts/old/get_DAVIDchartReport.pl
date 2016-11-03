#!/usr/bin/perl

##################################################################################################################################
### The original file is from DAVID http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_WebService.html
### Various operation (get and set functions) and web services are provided by DAVID for the Soap::Lite object
############################ Limitations #########################################################################################
### A job with more than 3000 genes to generate gene or term cluster report will not be handled by DAVID due to resource limit.
### No more than 200 jobs in a day from one user or computer.
##################################################################################################################################
### Modified by Manjusha Pande on 8/22/12
### Modified the script to take input gene list from .txt file(s) and generate a DAVID chart record file for each input file 
################################### Usage ######################################################################################## 
### > perl get_DAVIDchartReport.pl -dir [inDir] -type [idType] -thd [EASEScoreThreshold]
### inDir = path of the directory that contains input gene list (.txt file), header line (if any) should begin with #
### idType = DAIVID recognized id types of input ids; default: ENTREZ_GENE_ID
### EASEScoreThreshold = EASE score cut off for the enriched categories; default = 0.05
### Example: > perl get_DAVIDchartReport.pl -dir Input_geneLists -type ENTREZ_GENE_ID -thd 0.05
### Example: > perl scripts/get_DAVIDchartReport.pl -dir /ccmb/BioinfCore/Projects/Manjusha/Rsq_test_v2/noNovelTranscripts/cuffdiff/E7.5_WT_v_E7.5_mutant -type ENTREZ_GENE_ID -thd 0.05
##################################################################################################################################
 
use strict;
use warnings;
use SOAP::Lite; ## Perl's Web Services Toolkit Simple Object Access Protocol 
use HTTP::Cookies;
use Getopt::Long; ## Perl module to get command line arguments

# Get options  
my $inDir      = '';
my $idType        = 'ENTREZ_GENE_ID';
my $thd     = 0.05;
my $pipeline_email = "bfx-rnaseq-pipeline\@umich.edu";

GetOptions ( "dir=s"            => \$inDir,
             "type=s"           => \$idType,
			 "thd=s"            => \$thd );

#my $tmp = split("/", $inDir);
#my $dirName = $#tmp;

#my $outDir   = $inDir."_DAVIDchartReports/";
#if (-e "$outDir") 
#{	print "Warning: can not create output directory, directory $outDir already exists. \n"; # Check for existing Dir
#}
#else {
#	mkdir ($outDir);
#}

#=pod
my $soap = SOAP::Lite                             
     -> uri('http://service.session.sample')                
     -> proxy('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService',
                cookie_jar => HTTP::Cookies->new(ignore_discard=>1));

 #user authentication by email address
 #For new user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm
my $check = $soap->authenticate($pipeline_email)->result; # authenticate user by email address; return true if user has registered email with DAVID knowledge base
print "\nUser authentication: $check\n";

if (lc($check) eq "true") { 
	#list conversion types
	my $conversionTypes = $soap ->getConversionTypes()->result; # returns all acceptable idTypes
	# print  "\nConversion Types: \n$conversionTypes\n";  

	#list all annotation category names
	my $allCategoryNames= $soap ->getAllAnnotationCategoryNames()->result;	# returns all available annotation category names 	  	
	# print  "\nAll available annotation category names: \n$allCategoryNames\n"; 

	#list all background population names
        my $allBackgroundPopulations = $soap ->getAllPopulationNames()->result; 
	# print  "\nAll available background populations names: \n$allBackgroundPopulations\n"; 
	
	# Get input
	my @inFiles = glob ($inDir."/*_DEG_ids.txt");
	#print "input Dir is $inDir \n";
	#print "\@inFiles is @inFiles \n";
	for my $f (0..$#inFiles)
	{   	
		my @tmp1 = split (/\//, $inFiles[$f]);
		my @tmp2 = split (/\.txt/, $tmp1[$#tmp1]);
		my $inFile = $tmp2[0];
		
		print "Processing [".($f+1).'/'.($#inFiles+1)."]:$inFile... \n";
		open (INFILE, $inFiles[$f]);
		my @geneIDs = ();
		while(<INFILE>)
		{   my $line = $_;
			$line =~ s/\r|\n|\s//g;
			if ($line eq "" || $line =~ /^#/) {
				next;
			}
			if ($line =~ /\,/)
			{   my @tmp = split (/\,/, $line);
				push @geneIDs, @tmp;
			}else
			{   push @geneIDs, $line;
			}
		}   
		close INFILE;

		##### input list
		my $inputIds   = join (',', @geneIDs);
       		#print "Input ids are $inputIds \n";
		print "There are ". scalar (@geneIDs) . " genes in the input list.\n";
		#addList
		# my $inputIds = '31741_at,31734_at,32696_at,37559_at,41400_at,35985_at,39304_g_at,41438_at,35067_at,32919_at,35429_at,36674_at,967_g_at,36669_at,39242_at,39573_at,39407_at,33346_r_at,40319_at,2043_s_at,1788_s_at,36651_at,41788_i_at,35595_at,36285_at,39586_at,35160_at,39424_at,36865_at,2004_at,36728_at,37218_at,40347_at,36226_r_at,33012_at,37906_at,32872_at,989_at,32718_at,36957_at,32645_at,37628_at,33825_at,35687_at,32779_s_at,34493_at,31564_at,887_at,34712_at,32897_at,34294_at,41365_at,41446_f_at,34375_at,875_g_at,41099_at,919_at,38970_s_at,39159_at,34184_at,1018_at,38032_at,35956_s_at,35536_at,34562_at,1867_at,35957_at,39519_at,41657_at,38491_at,652_g_at,35776_at,34989_at,33455_at,39950_at,37723_at,31977_at,38629_at,34581_s_at,36210_g_at,35120_at,41532_at,37889_at,1332_f_at,40540_at,41105_s_at,1919_at,37542_at,39698_at,36711_at,36809_at,1167_s_at,31648_at,32364_at,40792_s_at,38685_at,41358_at,32931_at,35294_at,39870_at,38654_at,257_at,39071_at,35606_at,41726_at,33094_s_at,32405_at,1432_s_at,33698_at,408_at,39748_at,1953_at,36100_at,36101_s_at,1372_at,35314_at,40790_at,2030_at,179_at,1852_at,259_s_at,38024_at,35376_f_at,41779_at,39232_at,41159_at,40365_at,31626_i_at,40385_at,35613_at,37506_at,38207_at,887_at,600_at,1461_at,38691_s_at,1267_at,1177_at,1125_s_at,2036_s_at,31615_i_at,37283_at,40954_at,31758_at,36960_at,33143_s_at,37048_at,38538_at,1005_at,34963_at,39408_at,32464_at,706_at,1276_g_at,164_at,41445_at,40735_at,1891_at,1258_s_at,40856_at,1911_s_at,31562_at,32359_at,274_at,1804_at,41387_r_at,848_at,41499_at,39448_r_at,34537_at,36459_at,35500_at,37139_at,612_s_at,32133_at,39757_at,37629_at,38463_s_at,568_at,749_at,1939_at,38018_g_at,1857_at,32699_s_at,40661_at,1994_at,38373_g_at,33893_r_at,1388_g_at,35345_at,1385_at,36615_at,1263_at,37385_at,1774_at,37233_at,39753_at,32626_at,35915_at,35714_at,31669_s_at,36519_at,40473_at,1750_at,33751_at,37831_at,35472_at,41825_at,34666_at,35471_g_at,31888_s_at,37722_s_at,35414_s_at,39750_at,35726_at,37662_at,33802_at,352_at,31737_at,37938_at,36161_at,31558_at,34475_at,37223_at,38953_at,37857_at,189_s_at,41169_at,33092_at,38660_at,40895_g_at,37146_at,1936_s_at,38860_at,40210_at,41180_i_at,31586_f_at,33366_at,31521_f_at,762_f_at,1124_at,36009_at,41111_at,36749_at,37310_at,31522_f_at,35768_at,39421_at,39967_at,35992_at,38356_at,39331_at,34145_at,35378_at,199_s_at,35966_at,1866_g_at,37377_i_at,37378_r_at,833_at,31586_f_at,38062_at,34981_at,1569_r_at,1548_s_at,41446_f_at,36999_at,34226_at,33385_g_at,36173_r_at,1007_s_at,35149_at,38671_at,1973_s_at,37724_at,37317_at,33829_at,36532_at,39372_at,41717_at,38221_at,37418_at,33120_at,136_at,33492_at,1602_at,41505_r_at,41736_g_at,37862_at,31859_at,40913_at,35956_s_at,32193_at,1148_s_at,1244_at,38684_at,37440_at,32186_at,1242_at,39503_s_at,224_at,38374_at,36018_at,36603_at,33288_i_at,33662_at,33555_at,33539_at,430_at,471_f_at,1369_s_at,35372_r_at,38089_at,40310_at,41106_at,41216_r_at,32815_at,37463_r_at,33470_at,40522_at,1463_at,1743_s_at,1895_at,32583_at,35440_g_at,1091_at,1649_at,287_at,32119_at,131_at,38642_at,33922_at,35886_at,38326_at,38823_s_at,41088_at,41371_at,39841_at,32486_at,41234_at,41598_at,40478_at,37606_at,37170_at,34857_at,32062_at,37762_at,36052_at,40442_f_at,41550_at,36621_at,36929_at,38645_at,34438_at,39587_at,36562_at,37155_at,36055_at,36754_at,33545_at,1520_s_at,39402_at,32265_at,32679_at,1829_at,40669_at,31694_at,41382_at,41446_f_at,38391_at,34560_at,40098_at,32522_f_at,988_at,789_at,1270_at,1139_at,33665_s_at,1237_at,412_s_at,34688_at,31353_f_at,41856_at,32928_at,37584_at,32379_f_at,936_s_at,2082_s_at,36479_at,39175_at,32007_at,36103_at,37270_at,40840_at,37206_at,37365_at,37820_at,35848_at,37111_g_at,39522_at,36760_at,35018_at,31745_at,37424_at,36507_at,719_g_at,34165_at,41850_s_at';
		# my $idType = 'AFFYMETRIX_3PRIME_IVT_ID';
		# my $inputIds = '17381, 67717, 216343, 93726, 18430, 232680, 68800, 277898, 27222, 13004';		
		# my $idType = 'ENTREZ_GENE_ID';
		my $listName = $inFile;
		my $listType=0; # input list
		#to add background list, set listType=1
		my $list = $soap ->addList($inputIds, $idType, $listName, $listType)->result; # add a gene list or background list to current session
		$list *= 100;	# % of ids mapped												
		print "${list}% of list was mapped\n"; 
  	
		##### species
		#list all species  names
		#my $allSpecies= $soap ->getSpecies()->result;	# return species of the current list 	  	
		#print  "\nAll species: \n$allSpecies\n"; 
		#list current species names
		my $currentSpecies= $soap ->getCurrentSpecies()->result;	# return current species of the current list ??	  	
		#print  "\nCurrent species: \n$currentSpecies\n"; 
		#set user defined species 
		#my $species = $soap ->setCurrentSpecies("1")->result; # select specie(s) to use; argument is a string of integers delimited by comma
		#print "\nCurrent species: \n$species\n"; 
 
		##### annotation categories
		#set user defined categories 
		#my $categories = $soap ->setCategories("BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,UP_SEQ_FEATURE")->result;
		#to user DAVID default categories, send empty string to setCategories():
		my $categories = $soap ->setCategories("")->result; # select DAVID default categories; 
		# argument is a string with category names delimited by commas; return a list of validated category names
		#print "\nValid categories: \n$categories\n\n";  
        	
		##### EASE score threshold
		#my $thd=0.1;  ## threshold of EASE Score, a modified Fisher Exact P-Value, for gene-enrichment analysis.
		my $ct = 2; ## Minimum number of genes for the corresponding term

		my $chartReport = $soap->getChartReport($thd,$ct); # generate chart report; arguments include threshold, count
		my @chartRecords = $chartReport->paramsout; # @chartRecords is an array of n-1 (except the first one) $chartReport->result(hashes)
		
		my $n_records = @chartRecords;	
		print "Total chart records: $n_records\n\n";

		my $outFile = $inFile;
		$outFile =~ s/ids/DAVIDchartReport/;
		if (@geneIDs) { ### at least one gene id in the list
			open (chartReport, ">" . $inDir . "/" . $outFile . ".txt");
			print chartReport "Functional enrichment of ". scalar (@geneIDs) . " differentially expressed genes using DAVID (http://david.abcc.ncifcrf.gov)\n\n";
			if ($n_records) { ## at least one record
				print chartReport "#Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";
				my %chartRecord = %{$chartReport->result};  # this is the first record
				my $categoryName = $chartRecord{"categoryName"};
				my $termName = $chartRecord{"termName"};
				my $listHits = $chartRecord{"listHits"};
				my $percent = $chartRecord{"percent"};
				my $ease = $chartRecord{"ease"};
				my $Genes = $chartRecord{"geneIds"};
				my $listTotals = $chartRecord{"listTotals"};
				my $popHits = $chartRecord{"popHits"};
				my $popTotals = $chartRecord{"popTotals"};
				my $foldEnrichment = $chartRecord{"foldEnrichment"};
				my $bonferroni = $chartRecord{"bonferroni"};
				my $benjamini = $chartRecord{"benjamini"};
				my $FDR = $chartRecord{"afdr"};

					print chartReport "$categoryName\t$termName\t$listHits\t$percent\t$ease\t$Genes\t$listTotals\t$popHits\t$popTotals\t$foldEnrichment\t$bonferroni\t$benjamini\t$FDR\n";
				for my $j (0 .. (@chartRecords-1)) # these are remaining n-1 records
				{			
					%chartRecord = %{$chartRecords[$j]};
					$categoryName = $chartRecord{"categoryName"};
					$termName = $chartRecord{"termName"};
					$listHits = $chartRecord{"listHits"};
					$percent = $chartRecord{"percent"};
					$ease = $chartRecord{"ease"};
					$Genes = $chartRecord{"geneIds"};
					$listTotals = $chartRecord{"listTotals"};
					$popHits = $chartRecord{"popHits"};
					$popTotals = $chartRecord{"popTotals"};
					$foldEnrichment = $chartRecord{"foldEnrichment"};
					$bonferroni = $chartRecord{"bonferroni"};
					$benjamini = $chartRecord{"benjamini"};
					$FDR = $chartRecord{"afdr"};			
					print chartReport "$categoryName\t$termName\t$listHits\t$percent\t$ease\t$Genes\t$listTotals\t$popHits\t$popTotals\t$foldEnrichment\t$bonferroni\t$benjamini\t$FDR\n";				 
				}	
			}	
			else {print chartReport "No enrichment record found!!\n";}	
			#getChartReport 
			#print "\nchartReport.txt generated\n";
			close chartReport;
		}
	}
} 
print "All done!\n";
__END__
#=cut		
