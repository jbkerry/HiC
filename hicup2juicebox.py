#!/usr/bin/env python

import pysam, pandas as pd
import subprocess
from collections import OrderedDict

#samfile = pysam.AlignmentFile("/t1-data1/WTSA_Dev/jkerry/HiC/HiCUP/hicup_v0.5.9/Mifsud2015/CD34_HiC_1_2.hicup.trunc.bam", "rb")
#samfile = pysam.AlignmentFile("/t1-data1/WTSA_Dev/jkerry/HiC/Rao2014/HiC_firstRun/bowtie_results/bwt2/K562_batch1/SRR1658693_hg19.bwt2pairs.bam", "rb")
samfile = pysam.AlignmentFile("/t1-data1/WTSA_Dev/jkerry/HiC/HiCUP/hicup_v0.5.9/Rao2014/K562_1/SRR1658693_R1_2.hicup.bam", "rb")

#resfile = "/t1-data1/WTSA_Dev/jkerry/HiC/HiC-Pro_2.7.8/annotation/HindIII_resfrag_hg19.bed"
#resfile = "/t1-data1/WTSA_Dev/jkerry/HiC/HiC-Pro_2.7.8/annotation/DpnII_resfrag_hg19.bed"
#resdf = pd.read_table(resfile,sep="\t",header=None)
#
#def GetResFrag(readChr,readStart,mateChr,mateStart):
#    Exp,Chr,RFNum = resdf[(resdf[0]==readChr) & (resdf[1]<=int(readStart)) & (resdf[2]>int(readStart))][3].item().split('_')
#    MateExp,MateChr,MateRFNum = resdf[(resdf[0]==mateChr) & (resdf[1]<=int(mateStart)) & (resdf[2]>int(mateStart))][3].item().split('_')
#    return RFNum,MateRFNum

topDict = {}
#Counter = 1

for read in samfile.fetch(until_eof=True):
    #if Counter>100:
        #break
    if read.query_name not in topDict.keys():
        topDict[read.query_name] = {}
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
        
        #readResSite,mateResSite = GetResFrag(read.reference_name,read.reference_start,read.next_reference_name,read.next_reference_start)
        readResSite = mateResSite = 3
        ChrNum = read.reference_name[3:]
        MateChrNum = read.next_reference_name[3:]
        if ChrNum=="X":
            ChrNum = 24
        elif ChrNum=="Y":
            ChrNum = 25
        elif ChrNum=="M":
            ChrNum = 23
            
        if MateChrNum=="X":
            MateChrNum = 24
        elif MateChrNum=="Y":
            MateChrNum = 25
        elif MateChrNum=="M":
            MateChrNum = 23
        MAPQF = read.mapping_quality
        ReadList = [read.reference_name, read.reference_start, readResSite, readStr, read.mapping_quality]
        if int(ChrNum)<=int(MateChrNum):
            topDict[read.query_name].update({"read1": ReadList})
        elif int(MateChrNum)<int(ChrNum):
            topDict[read.query_name].update({"read2": ReadList})
    else:
        ReadList = [read.reference_name, read.reference_start, mateResSite, mateStr, read.mapping_quality]
        if int(ChrNum)<=int(MateChrNum):
            topDict[read.query_name].update({"read2": ReadList})
        elif int(MateChrNum)<int(ChrNum):
            topDict[read.query_name].update({"read1": ReadList})
    #Counter=Counter+1

sortedDict = OrderedDict(sorted(topDict.iteritems(), key=lambda x: (x[1]["read1"][0],x[1]["read2"][0])))
outFile = open("K562_batch1_HiCUP.jb.txt","w")
for i in sortedDict:
    #outFile.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}".format(i,sortedDict[i]["read1"][3],sortedDict[i]["read1"][0],sortedDict[i]["read1"][1],sortedDict[i]["read1"][2],sortedDict[i]["read2"][3],sortedDict[i]["read2"][0],sortedDict[i]["read2"][1],sortedDict[i]["read2"][2],sortedDict[i]["read1"][4],sortedDict[i]["read2"][4]))
    outFile.write("{0} {1} {2} {3} 0 {4} {5} {6} 1 {7} {8}\n".format(i,sortedDict[i]["read1"][3],sortedDict[i]["read1"][0],sortedDict[i]["read1"][1],sortedDict[i]["read2"][3],sortedDict[i]["read2"][0],sortedDict[i]["read2"][1],sortedDict[i]["read1"][4],sortedDict[i]["read2"][4]))
outFile.close()