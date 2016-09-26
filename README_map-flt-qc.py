
### Map reads with BWA
import psychiplib.mapping as mp
test_map = mp.BwaMapper("test_data/R1.fastq.gz", # fastq read 1
                        "test_data/R2.fastq.gz", # fastq read 2
                        "bwa_index/hg38_sponge", # bwa index
                        "test_out")              # output directory
test_map.run(16) # CPU for BWA
                 # Output: R1.raw.sam.gz
                 # Format: <prefix>.raw.sam.gz


### Post alignment filtering
import psychiplib.post_filter as pf
test_filter = pf.PostFilter("test_out/R1.raw.sam.gz", # raw sam file
                            "test_out")               # output directory

## Step 1: remove reads with bad CIGAR strings
test_filter.remove_badcigar(16) # Threads for SAMtools
                                # Output: R1.raw.bam
                                # Format: <prefix>.raw.bam

## Step 2: remove unmapped, low quality, orphan, and supplemental reads
test_filter.remove_aberrant_reads(16) # Threads for SAMtools
                                      # Output: R1.flt.bam
                                      # Format: <prefix>.flt.bam

## Step 3: remove PCR duplicates and sponge mapping reads
##         (currently, "java -Xmx4G". Please change it if needed)
test_filter.remove_artifacts("path/to/picard.jar",  # path to picard.jar
                             "path/to/sponges.txt") # sponge name list
                                                    # Output: R1.final.[bam|bai], R1.pcrDups.QC.txt

## Step 4: generate BED and BEDPE files
##         (those files are *unsorted*)
test_filter.generate_bed(16) # Threads for SAMtools
                             # Output: R1.bed.gz (for pseudorep)
                             # Output: R1.bedpe.gz (bedtool raw output)

### QC part
import psychiplib.qc as qc
test_qc = qc.QC("test_out") # output direcotry

## Step 1: Get map stats
test_qc.get_mapstats("test_out/R1.raw.bam",   # Output: R1.raw.flagstat.QC.txt
                     "test_out/R1.final.bam") # Output: R1.final.flagstat.QC.txt

## Step 2: Compute library complexity
test_qc.get_library_complexity("test_out/R1.flt.bam", 16) # Threads for SAMtools
                                                          # Output: R1.libComplexity.QC.txt

## Step 3: *!NOT TESTED!* Cross correlation coefficient
test_qc.cross_correlation("test_out/R1.bedpe.gz",                  # BEDPE file
                          "phantompeakqualtools/run_spp_nodups.R", # phantompeakqualtools
                          32)                                      # Number of CPU for phantompeakqualtools

## Clean up intermediate files (if needed)
## It removes R1.raw.sam.gz, R1.raw.bam, R1.flt.bam
test_filter.clean()



