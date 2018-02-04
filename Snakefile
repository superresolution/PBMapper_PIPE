
"""
Author: Hao Chang
Affiliation: Yale
Aim: A simple Snakemake workflow to process paired-end stranded PB-Mapper.
Date: Mon Jan 30th 2018
Run: snakemake   -s Snakefile   
Latest modification: 
  - todo
"""
# package

import os


# functions

def message(ip):
	print (ip)

def hamdist(str1, str2):
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs





# This should be placed in the Snakefile.

##-----------------------------------------------##
## Working directory                             ##
## Adapt to your needs                           ##
##-----------------------------------------------##

BASE_DIR = "/media/changhao/MYSTORE/Project_4_B16_screen"
WDIR = BASE_DIR + "/B16_screen"
INDEX_MOUSE = WDIR + "/index/mm10"

# workdir: WDIR
#message("The current working directory is " + WDIR)



##--------------------------------------------------------------------------------------##
## Variables declaration                          
## Declaring some variables used by topHat and other tools... 
## (GTF file, INDEX, chromosome length)
##--------------------------------------------------------------------------------------##
# Adapt the path to your needs



##--------------------------------------------------------------------------------------##
## The list of samples to be processed
##--------------------------------------------------------------------------------------##

R1, = glob_wildcards("1_or_seq/R1_{smp}.fastq.gz")
R2, = glob_wildcards("1_or_seq/R2_{smp}.fastq.gz")

NUM = 100

#message(R1)
#message(R2)

##

from Bio import SeqIO
from Bio import Seq
import gzip
import operator

# parameter

Q5 = ["TAATGTGG",
              "GCACTCAG", 
              "AACAGCGG", 
              "CCATATGA", 
              "TGGAAAGC", 
              "AGCAACGC", 
              "CCCTTGCA", 
              "CCCTCTTG", 
              "TTCGAGCC", 
              "AGTAGTTA"]

Q7R = ["ATGCGGAT",
        "ACTTCAAT",
        "ACGCAGAG",
        "TGCGTATC",
        "CGAACGTA",
        "AGGATTCA",
        "TTATAGCC",
        "AGGTTGTG",
        "GAGGCTGA",
        "AACAATCC"]


Q5_name = []
Q7_name = []
filenames = []



for i in range(1,11):
    Q5_name.append("Q5_" + str(i))
    Q7_name.append("Q7_" + str(i))
for i in range(0,100):
    filenames.append(Q5_name[i//10] + "_" + Q7_name[i%10] + ".fastq")


##

rule all:
    input:
        # '2_merge/R1.fastq.gz',
        # '2_merge/R2_R.fastq.gz',
        expand('4_joint_spliter_trimmer/{fn}', fn = filenames),
        # expand('6_bowtie2/{fn}.sam', fn = filenames),       
        # expand('7_bowtie2_dedup/{fn}', fn = bamdedup),
        expand('8_bowtie2_dedup_sorted/{fn}_dedup_sorted.bam', fn = filenames),
        expand('9_bigwigfile/{fn}_dedup_sorted.bw', fn = filenames),
        expand('10_bedfile/{fn}_dedup_sorted.bed', fn = filenames)

## merge fastq.gz files to R1 and R2

rule merge_R1:
    input:
        r1 = expand('1_or_seq/R1_{smp}.fastq.gz', smp=R1)
    output:
        R1_merged = '2_merge/R1.fastq.gz'

    shell:
        'cat {input.r1} > {output.R1_merged}'

rule merge_R2:
    input:
        r2 = expand('1_or_seq/R2_{smp}.fastq.gz', smp=R2)
    output:
        R2_merged = '2_merge/R2.fastq.gz'

    shell:
        'cat {input.r2} > {output.R2_merged}'

##

rule reverse_complement:
    input:
        R2 = '2_merge/R2.fastq.gz'
    output:
        R2_R = '2_merge/R2_R.fastq.gz'

    shell:
        'gunzip -c {input.R2} \
        | fastx_reverse_complement -z -o {output.R2_R}'


## joint_split_trimmer

rule joint_spliter_trimmer:
    input:
        R1 = '2_merge/R1.fastq.gz',
        R2 = '2_merge/R2_R.fastq.gz'
    output:
        expand('4_joint_spliter_trimmer/{fn}', fn = filenames)

    run:
        # package

        from Bio import SeqIO
        from Bio import Seq
        import gzip
        import operator
        import os


        # parameter

        directory = "4_joint_spliter_trimmer/"

        if not os.path.exists(directory):
            os.makedirs(directory)

        filedata = {filename: open("4_joint_spliter_trimmer/" + filename, 'w') for filename in filenames}

        mismatch = 3
        plus = 37
        minus = -28


        i = 1
        rec2 = SeqIO.parse(gzip.open(input[1], 'rt'),"fastq")
        # print (input[:])

        for record in SeqIO.parse(gzip.open(input[0], 'rt'),"fastq"): 
            if i < 100:
                record2 = rec2.__next__() # for python 2 is next()
                # if record2.id == record.id:
                # print "yes"
                sequence = str(record.seq)
                letter_annotations = record.letter_annotations
                sequence2 = str(record2.seq)
                letter_annotations2 = record2.letter_annotations

                # 
                Q5_diff = []
                Q7_diff = []

                for index, value in enumerate(Q5):
                    Q5_diffs = hamdist(sequence[:8], value)
                    Q5_diff.append(Q5_diffs) 

                for index, value in enumerate(Q7R):
                    Q7_diffs = hamdist(sequence2[-8:], value)
                    Q7_diff.append(Q7_diffs)

                Q5_min_index, Q5_min_value = min(enumerate(Q5_diff), key=operator.itemgetter(1))
                Q7_min_index, Q7_min_value = min(enumerate(Q7_diff), key=operator.itemgetter(1))
                
                # print (Q5_min_index, Q5_min_value, Q7_min_index, Q7_min_value)

                if (Q5_min_value < mismatch) & (Q7_min_value < mismatch):
                    # print ("yes")
                    record.letter_annotations = {}
                    new_sequence = sequence[plus:] + sequence2[:minus]
                    record.seq = Seq.Seq(new_sequence)
                    new_letter_annotations = {'phred_quality': (letter_annotations['phred_quality'][plus:] + letter_annotations2['phred_quality'][:minus])}
                    record.letter_annotations = new_letter_annotations

                    f_name = filenames[Q5_min_index*10 + Q7_min_index]
                    SeqIO.write(record, filedata[f_name], "fastq")
                    # print (record)
            else:
                break

            # i = i + 1 # just for program test
        
        print ('Out of loop')

        for file in filedata.values():
            file.close()

# joint_split_trimmer_TTAA_detection

rule joint_spliter_trimmer_withttaa:
    input:
        R1 = '2_merge/R1.fastq.gz',
        R2 = '2_merge/R2_R.fastq.gz'
    output:
        expand('5_joint_spliter_trimmer_withttaa/{fn}', fn = filenames)

    run:
        # package

        from Bio import SeqIO
        from Bio import Seq
        import gzip
        import operator
        import os


        # parameter

        directory = "5_joint_spliter_trimmer_withttaa/"

        if not os.path.exists(directory):
            os.makedirs(directory)

        filedata = {filename: open("5_joint_spliter_trimmer_withttaa/" + filename, 'w') for filename in filenames}

        mismatch = 3
        plus = 37
        minus = -28
        left = 2
        right = 12


        i = 1
        rec2 = SeqIO.parse(gzip.open(input[1], 'rt'),"fastq")
        # print (input[:])

        for record in SeqIO.parse(gzip.open(input[0], 'rt'),"fastq"): 
            if i < 100:
                record2 = rec2.__next__() # for python 2 is next()
                # if record2.id == record.id:
                # print "yes"
                sequence = str(record.seq)
                letter_annotations = record.letter_annotations
                sequence2 = str(record2.seq)
                letter_annotations2 = record2.letter_annotations

                # 
                Q5_diff = []
                Q7_diff = []

                for index, value in enumerate(Q5):
                    Q5_diffs = hamdist(sequence[:8], value)
                    Q5_diff.append(Q5_diffs) 

                for index, value in enumerate(Q7R):
                    Q7_diffs = hamdist(sequence2[-8:], value)
                    Q7_diff.append(Q7_diffs)

                Q5_min_index, Q5_min_value = min(enumerate(Q5_diff), key=operator.itemgetter(1))
                Q7_min_index, Q7_min_value = min(enumerate(Q7_diff), key=operator.itemgetter(1))
                
                # print (Q5_min_index, Q5_min_value, Q7_min_index, Q7_min_value)

                str1 = sequence[(plus-left):(plus+right)]

                if (Q5_min_value < mismatch) & (Q7_min_value < mismatch) & (str1.find("TTAA") > -1):
                    # print ("yes")
                    record.letter_annotations = {}
                    new_sequence = sequence[plus:] + sequence2[:minus]
                    record.seq = Seq.Seq(new_sequence)
                    new_letter_annotations = {'phred_quality': (letter_annotations['phred_quality'][plus:] + letter_annotations2['phred_quality'][:minus])}
                    record.letter_annotations = new_letter_annotations

                    f_name = filenames[Q5_min_index*10 + Q7_min_index]
                    SeqIO.write(record, filedata[f_name], "fastq")
                    # print (record)
            else:
                break

            # i = i + 1 # just for program test
        
        print ('Out of loop')

        for file in filedata.values():
            file.close()

## 

rule bowtie2:
    input:
        # idx = '/media/changhao/MYSTORE/Project_4_B16_screen/B16_Screen/index/mm10',
        expand('5_joint_spliter_trimmer_withttaa/{fn}', fn = filenames)
    output:
        expand('6_bowtie2/{fn}.sam', fn = filenames)
    run:
        for i in range(0,num):
            r1 = input[i]
            o1 = output[i]
            shell('bowtie2 -x /media/changhao/MYSTORE/Project_4_B16_screen/B16_Screen/index/mm10 -U {r1} -S {o1}')

##

rule sam_to_bam:
    input:
        expand('6_bowtie2/{fn}.sam', fn = filenames)
    output:
        expand('6_bowtie2_bam/{fn}.bam', fn = filenames)
    run:
        num = 100
        for i in range(0,num):
            r1 = input[i]
            o1 = output[i]
            shell('samtools view -S -b {r1} > {o1}')
##

rule dedup:
    input:
        R1 = expand('6_bowtie2_bam/{fn}.bam', fn = filenames)
    output:
        expand('7_bowtie2_dedup/{fn}_dedup.bam', fn = filenames)
    run:
        num = 100
        for i in range(0,num):
            r1 = input.R1[i]
            o1 = output[i]
            shell('samtools rmdup {r1} {o1}')

##

rule sorted_bam:
    input:
        R1 = expand('7_bowtie2_dedup/{fn}_dedup.bam', fn = filenames)
    output:
        expand('8_bowtie2_dedup_sorted/{fn}_dedup_sorted.bam', fn = filenames)
    run:
        num = 100
        for i in range(0,num):
            r1 = input.R1[i]
            o1 = output[i]
            shell('samtools sort {r1} > {o1}')

##

rule index_bam:
    input:
        R1 = expand('8_bowtie2_dedup_sorted/{fn}_dedup_sorted.bam', fn = filenames)
    output:
        expand('8_bowtie2_dedup_sorted/{fn}_dedup_sorted.bam.bai', fn = filenames)
    run:
        num = 100
        for i in range(0,num):
            r1 = input.R1[i]
            o1 = output[i]
            shell('samtools index {r1}')


##

rule sorted_bam_to_bw: 
    input:
        R = expand('8_bowtie2_dedup_sorted/{fn}_dedup_sorted.bam', fn = filenames),
        bai = expand('8_bowtie2_dedup_sorted/{fn}_dedup_sorted.bam.bai', fn = filenames)
    output:
        expand('9_bigwigfile/{fn}_dedup_sorted.bw', fn = filenames)
    run:
        num = 100
        for i in range(0,num):
            r1 = input.R[i]
            o1 = output[i]
            shell('bamCoverage -b {r1} -o {o1}')

##

rule sorted_bam_to_bed:
    input:
        R = expand('8_bowtie2_dedup_sorted/{fn}_dedup_sorted.bam', fn = filenames),
        bai = expand('8_bowtie2_dedup_sorted/{fn}_dedup_sorted.bam.bai', fn = filenames)
    output:
        expand('10_bedfile/{fn}_dedup_sorted.bed', fn = filenames)
    run:
        num = 100
        for i in range(0,num):
            r1 = input.R[i]
            o1 = output[i]
            shell('bedtools bamtobed -i {r1} > {o1}')
