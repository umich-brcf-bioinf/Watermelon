import binascii

'''This rule assures that all concatenated fastqs are gzipped.
This assures that downstream rules aren't confused'''
rule align_standardize_gz:
    input:
        ALIGNMENT_DIR + "01-raw_reads/{sample}_R{read}.fastq.gz", # TWS - Despite the suffix, these may or may not be gzipped
    output:
        ALIGNMENT_DIR + "02-gz_reads/{sample}_R{read}.fastq.gz",
    threads: 1
    log:
        JOB_LOG_DIR + "align_standardize_gz_{sample}_R{read}.log"
    resources: time_min=240
    run:
        #TWS - This assumption will not work if a sample has mixed gzipped / plaintext fastqs!
        #Addressed via InputFastqManager._get_sample_fastq_paths throwing an error in that case
        with open(str(input), 'rb') as test_f:
            if binascii.hexlify(test_f.read(2)) == b'1f8b':
                #Is compressed. Just symlink
                os.symlink(os.path.abspath(str(input)), os.path.abspath(str(output)))
            else:
                #Is uncompressed. Use Gzip to compress
                cmd = "gzip -c {} > {}".format(input, output)
                shell(cmd)
