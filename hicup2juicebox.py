#!/usr/bin/env python

import pysam, pandas as pd

samfile = pysam.AlignmentFile("/t1-data1/WTSA_Dev/jkerry/HiC/HiCUP/hicup_v0.5.9/Mifsud2015/CD34_HiC_1_2.hicup.bam", "rb")
resfile = "/t1-data1/WTSA_Dev/jkerry/HiC/HiC-Pro_2.7.8/annotation/HindIII_resfrag_hg19.bed"
resdf = pd.read_table(resfile,sep="\t",header=None)

def GetResFrag(readChr,readStart,mateChr,mateStart):
    Exp,Chr,RFNum = resdf[(resdf[0]==readChr) & (resdf[1]<=int(readStart)) & (resdf[2]>int(readStart))][3].item().split('_')
    MateExp,MateChr,MateRFNum = resdf[(resdf[0]==mateChr) & (resdf[1]<=int(mateStart)) & (resdf[2]>int(mateStart))][3].item().split('_')
    return RFNum,MateRFNum
    
outFile = open("CD34_HiC_1_2.hicup.JB.txt","w")
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
        
        readResSite,mateResSite = GetResFrag(read.reference_name,read.reference_start,read.next_reference_name,read.next_reference_start)    
        
        #if Counter<=20:
        outFile.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} ".format(read.query_name,readStr,read.reference_name,read.reference_start,readResSite,mateStr,read.next_reference_name,read.next_reference_start,mateResSite,read.mapping_quality))
            #Counter=Counter+1
    else:
        outFile.write("{0}\n".format(read.mapping_quality))
        #if Counter>20:
            #break
outFile.close()

#print(resdf[resdf[0]=="chr1"][resdf[1]<20000])
#if Counter==1:
#Exp,Chr,RFNum = resdf[(resdf[0]==read.reference_name) & (resdf[1]<=int(read.reference_start)) & (resdf[2]>int(read.reference_start))][3].item().split('_')
#readFragID = resdf[(resdf[0]==read.reference_name) & (resdf[1]<=int(read.reference_start)) & (resdf[2]>int(read.reference_start))][3].item()
#[(resdf[1]>=int(read.reference_start)) & (resdf[2]<int(read.reference_start))]
#print(Exp+"_"+Chr+"_"+RFNum+" : "+RFNum)
#print(readFragID)