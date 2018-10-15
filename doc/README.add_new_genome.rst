Add a new genome
================

This doc shows the steps used to add hg38 genome to Watermelon.
These same steps can be adapted to any new genome.

You will need:

- About 3 hours to download and unpack the genome (depending on size).
- About 2 hours and 150Gb RAM to build the HISAT2 indices (again, depending on
  size; FWIW, GCRh38 took ).

Establish directory and download sequence/annotation data
---------------------------------------------------------

Download/unpack genome::

  cd /nfs/med-bfx-common/pipelines/Watermelon/references
  mkdir hg38
  cd hg38
  mkdir igenome
  cd igenome
  # (This will take several hours.)
  wget ftp://igenome:G3nom3s4u@ussd-
  ftp.illumina.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
  # (This will also take awhile.)
  tar -xvzf Homo_sapiens_UCSC_hg38.tar.gz

Prep Watermelon dirs::

  cd ..
  mkdir Sequence
  cd Sequence/
  cp -lr ../igenome/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/ .
  mv Bowtie2Index/genome.fa .
  chmod 664 genome.fa
  module load samtools
  samtools faidx genome.fa
  cd Bowtie2Index/
  ln -s ../genome.fa
  ln -s ../genome.fa.fai
  cd ../..
  ln igenome/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf
  chmod 664 genes.gtf

Capture the steps from above as a README.txt file and confirm your top-level results look correct::

  ls -l
  -rw-rw-r-- 2 cgates Isilon-bfxcore 151467494 Aug 18  2015 genes.gtf
  drwxrwsr-x 3 cgates Isilon-bfxcore       105 Oct 14 14:10 igenome
  -rw-rw-r-- 1 cgates Isilon-bfxcore       682 Oct 14 14:20 README.txt
  drwxrwsr-x 3 cgates Isilon-bfxcore        88 Oct 14 14:05 Sequence

Build HISAT2 indices
--------------------

Add a new block to config/genome_references.yaml::

  hg38:
      genome: hg38
      fastq_screen:
          species: human
      hisat2_index_spliced: true
      references:
          fasta: /nfs/med-bfx-common/pipelines/Watermelon/references/hg38/Sequence/genome.fa
          gtf: /nfs/med-bfx-common/pipelines/Watermelon/references/hg38/genes.gtf
          bowtie2_index: /nfs/med-bfx-common/pipelines/Watermelon/references/hg38/Sequence/Bowtie2Index
          entrez_gene_info: /nfs/med-bfx-common/pipelines/Watermelon/references/entrez_gene_info/2016_09_02/gene_info
          hisat2_index: /nfs/med-bfx-common/pipelines/Watermelon/references/hg38/HISAT2_index/hg38_genes

Run rule to build HISAT2 indices::

  conda activate watermelon
  snakemake --use-conda --snakefile /nfs/med-bfx-common/pipelines/Watermelon/hisat2_index.smk --cores 8 -p

Note that you can build the indices for multiple genomes at once by setting the
cores to 16 for 2 genomes, 24 for three, etc. N.B. a reasonably large
genome on 8 cores can take 1.5 hours and 150Gb RAM to build, so be deliberate.

Benchmarks
----------

You can review benchmarks like so::

  cd /ccmb/BioinfCore/Common/pipelines/Watermelon/references
  awk 'BEGIN {OFS="\t"} NR==1 {print "#file",$0} FNR>1 {split(FILENAME, a, "/"); print a[6],$0}' `find . -maxdepth 5 -path '*/benchmarks/*' -name '*benchmark.txt' | sort` | column -t

  #file                                              s           h:m:s    max_rss    max_vms    max_uss    max_pss    io_in     io_out     mean_load
  align_create_hisat2_index.ce10.benchmark.txt       242.2707    0:04:02  13695.97   14586.54   13691.81   13692.32   97.55     2857.31    0.00
  align_create_hisat2_index.ce11.benchmark.txt       206.4044    0:03:26  5446.36    9531.71    5442.03    5442.57    97.58     2874.23    0.00
  align_create_hisat2_index.ecoMG1655.benchmark.txt  2.0948      0:00:02  48.01      1321.11    45.46      45.61      5.69      25.44      0.00
  align_create_hisat2_index.ecoUTI89.benchmark.txt   1.4071      0:00:01  16.43      569.97     10.38      12.17      0.00      0.00       0.00
  align_create_hisat2_index.GRCh37.benchmark.txt     7192.2192   1:59:52  156094.35  161058.81  156091.22  156091.38  3007.23   81670.04   0.00
  align_create_hisat2_index.GRCh38.benchmark.txt     25859.3212  7:10:59  326911.32  333640.07  326908.18  326908.34  68777.61  158933.12  0.00
  align_create_hisat2_index.GRCz10.benchmark.txt     2854.1559   0:47:34  75249.83   80204.65   75245.82   75246.42   1465.97   38305.20   0.00
  align_create_hisat2_index.hg19.benchmark.txt       6787.2697   1:53:07  156088.25  162661.54  156084.02  156084.61  13079.12  81587.84   0.00
  align_create_hisat2_index.hg38.benchmark.txt       6365.8881   1:46:05  158338.64  163503.16  158334.55  158335.28  7.66      83703.42   0.00
  align_create_hisat2_index.mm10.benchmark.txt       5584.6798   1:33:04  149324.60  154711.56  149320.48  149321.07  2731.88   75414.35   0.00
  align_create_hisat2_index.rn5.benchmark.txt        5566.8741   1:32:46  146936.95  152294.62  146932.88  146933.48  19608.36  73093.34   0.00
  align_create_hisat2_index.rn6.benchmark.txt        6058.8278   1:40:58  149589.25  155123.49  149585.15  149585.88  3951.16   75507.94   0.00
  align_create_hisat2_index.WBS235.benchmark.txt     204.6275    0:03:24  5434.53    8302.07    5430.21    5430.75    97.24     2878.59    0.00

- Note that hg38 took just under 2 hours (h:m:s) to build and just under 160Gb
  RAM (max_rss).
- Also, note GRCh38 took ~7 hours to build and ~326Gb RAM.
