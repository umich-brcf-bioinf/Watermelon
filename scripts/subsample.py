#!/usr/bin/python

### From Simon Anders http://seqanswers.com/forums/showthread.php?t=12070
### Python tutorial at http://docs.python.org/2/tutorial/
### Usage: python subsample.py <fraction> <input file 1> <input file 2> <output file 1> <output file 2>
### where <fraction> is a number between 0 and 1, giving the sampling fraction
### Example: python /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/SubSampleFastq/subsample.py 0.01 /ccmb/BioinfCore/Projects/Hu_826_pathu_mceachin_mpande_RS2/reads/Sample_28339_R1.fastq.gz /ccmb/BioinfCore/Projects/Hu_826_pathu_mceachin_mpande_RS2/reads/Sample_28339_R2.fastq.gz /ccmb/BioinfCore/Projects/Hu_826_pathu_mceachin_mpande_RS2/reads/subsample/Sample_28339_R1_25.fastq.gz /ccmb/BioinfCore/Projects/Hu_826_pathu_mceachin_mpande_RS2/reads/subsample/Sample_28339_R2_25.fastq.gz


import sys, random, itertools
import HTSeq

fraction = float( sys.argv[1] )
in1 = iter( HTSeq.FastqReader( sys.argv[2] ) )
in2 = iter( HTSeq.FastqReader( sys.argv[3] ) )
out1 = open( sys.argv[4], "w" )
out2 = open( sys.argv[5], "w" )

for read1, read2 in itertools.izip( in1, in2 ):
   if random.random() < fraction:
      read1.write_to_fastq_file( out1 )
      read2.write_to_fastq_file( out2 )
      
out1.close()
out2.close()

### for SE data

## import sys, random, itertools
## import HTSeq

## fraction = float( sys.argv[1] )
## in1 = iter( HTSeq.FastqReader( sys.argv[2] ) )
## out1 = open( sys.argv[3], "w" )

## for read1 in in1:
  ##  if random.random() < fraction:
      ## read1.write_to_fastq_file( out1 )

## out1.close()
