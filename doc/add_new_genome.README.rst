Add a new genome
================

This doc shows the steps used to add hg38 genome to Watermelon.
These same steps can be adapted to any new genome.

You will need:

- About 3 hours to download and unpack the genome (depending on size).
- About 2 hours and 150Gb RAM to build the HISAT2 indices (again, depending on size).

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
