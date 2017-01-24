#!/usr/bin/env python

import subprocess

#import pysam
#samfile = pysam.AlignmentFile("/t1-data1/WTSA_Dev/jkerry/HiC/HiCUP/hicup_v0.5.9/Mifsud2015/CD34_HiC_1_2.hicup.sorted.bam", "rb")
#
#variable = samfile.fetch("chr1", 10000, 200000)
#for x in variable:
#    print (str(x))
    
samfile = "/t1-data1/WTSA_Dev/jkerry/HiC/HiCUP/hicup_v0.5.9/Mifsud2015/CD34_HiC_1_2.hicup.bam"

ViewFile = subprocess.Popen("samtools view "+samfile, shell=True, stdout=subprocess.PIPE)
#FileToVar = ViewFile.stdout.read().rstrip('\n')
FileToVar = ViewFile.stdout.read()