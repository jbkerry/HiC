#!/usr/bin/env python

import pysam
samfile = pysam.AlignmentFile("/t1-data1/WTSA_Dev/jkerry/HiC/HiCUP/hicup_v0.5.9/Mifsud2015/CD34_HiC_1_2.hicup.bam", "rb")

variable = samfile.fetch("chr1", 10, 20)
for x in variable:
    print (str(x))