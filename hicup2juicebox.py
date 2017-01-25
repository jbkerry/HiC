#!/usr/bin/env python

#import subprocess

import pysam
samfile = pysam.AlignmentFile("/t1-data1/WTSA_Dev/jkerry/HiC/HiCUP/hicup_v0.5.9/Mifsud2015/CD34_HiC_1_2.hicup.bam", "rb")

outFile = open("ForJB.txt","w")
Counter = 1    
for read in samfile.fetch(until_eof=True):
    if read.is_read1:
        readStr = ""
        mateStr = ""
        if read.is_reverse:
            readStr = "1"
        else:
            readStr = "0"
            
        if read.mate_is_reverse:
            mateStr = "1"
        else:
            mateStr = "0"
        if Counter<=20:
            outFile.write("{0} {1} {2} {3} resFrag1 {4} {5} {6} resFrag2 {7} ".format(read.query_name,readStr,read.reference_name,read.reference_start,mateStr,read.next_reference_name,read.next_reference_start,read.mapping_quality))
            Counter=Counter+1
    else:
        outFile.write("{0}\n".format(read.mapping_quality))
        if Counter>20:
            break
outFile.close()