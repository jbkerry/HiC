#!/usr/bin/env python

from __future__ import division
import numpy as np
import pandas as pd
import math
import numpy.random
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata

### Calculate observed contact matrix

df = pd.DataFrame()
Matrix = [MatrixLine.rstrip('\n') for MatrixLine in open('/t1-data1/WTSA_Dev/jkerry/HiC/Rao2014/HiC_GM12878/hic_results/matrix/GM12878_batch1/raw/500000/GM12878_batch1_500000_chr1.matrix')]
#Matrix = [MatrixLine.rstrip('\n') for MatrixLine in open('datHalf.txt')]
PreviousY = 0
PreviousX = 0
for ThisMatLine in Matrix:
    x, y, z = ThisMatLine.split('\t')
    x = int(x)
    y = int(y)
    z = int(z)
    
    if (PreviousY<y) & (PreviousY!=(y-1)):
        difference = y-PreviousY
        counter=1
        while counter<=difference:
            df.loc[(PreviousY+counter),x] = 0
            if (PreviousY+counter)!=x:
                df.loc[x,(PreviousY+counter)] = 0
            counter=counter+1
    elif (PreviousY>y) & (PreviousY!=499):
        difference = 499-PreviousY
        counter=1
        while counter<=difference:
            df.loc[(PreviousY+counter),PreviousX] = 0
            if (PreviousY+counter)!=PreviousX:
                df.loc[PreviousX,(PreviousY+counter)] = 0
            counter=counter+1
    df.loc[y,x]=z
    if x!=y:
        df.loc[x,y]=z
    PreviousX=x
    PreviousY=y

#### Calculate Vanilla Coverage
#SumList = []
#RowCounter = 1
#while RowCounter<=499:
#    SumList.append(df[RowCounter].sum())
#    RowCounter+=1
#    
#print(max(SumList))
#SumList = SumList/max(SumList)
#print(max(SumList))




df_VC = df
#df_VC["sum"] = df_VC.sum(axis=1)
#df_VC.loc[500] = df_VC.sum(axis=0)
#df_VC_new=df_VC.div(df_VC["sum"], axis=0)
#df_VC_new=df_VC_new.div(df_VC.loc[500], axis=1)
i_pos = 0
while i_pos<=498:
    j_pos = 0
    while j_pos<=498:
        df_VC.iloc[i_pos,j_pos] = (df.loc[i_pos+1,j_pos+1]*10000000000)/(df[i_pos+1].sum()*df[j_pos+1].sum())
        j_pos+=1
    i_pos+=1
#print(df_VC[:10])
#print(df_VC[495:])
#print(df_VC_new[:10])
#new_df = (df.T/df.T.sum()).T
#new_df["sum"] = new_df.sum(axis=1)
#print(new_df[:10])
#print(df[:10])
#print(new_df[2].sum())
#print(new_df[5].sum())
#corrected_df = df.div(new_df["sum"],axis=0)
#print(corrected_df[:10])
#print(df.loc[2,2]/(sum(df[2])*sum(df[2])))

        

#df_sub = df.loc[:200,:200]
df.to_csv("chr1MatrixFull_GM12878.csv",sep="\t")

#outfile = open("chr1Matrix.txt","w")
#matCounter = 1
#while matCounter<=100:
#    outfile.write(str(df.loc[matCounter][1:100]))
#    outfile.write("\n")
#    matCounter=matCounter+1
#outfile.close()

### Calculate expected contact matrix

distance = 500000
chrlength = df.shape[0]*distance

#exp_df = pd.DataFrame()
colNames = df.columns.values
rowNames = df.index.values
###for j in rowNames:
###    for i in colNames:
###        idist = i*distance
###        jdist = j*distance
###        #interaction_pos = math.ceil((chrlength-distance/2)/((abs(idist-jdist)/distance)+1))
###        #interaction_pos = math.ceil(((chrlength-distance/2)/distance)/((abs(idist-jdist)/distance)+1))
###        interaction_pos = math.ceil((chrlength-distance/2)/(abs(idist-jdist)+distance))
###        #pos_square = interaction_pos*interaction_pos
###        pos_square = interaction_pos*(interaction_pos+1)/2
###        exp_df.loc[j,i]=pos_square
###        #print("i = {0}, j = {1}, idist = {2}, jdist = {3}, int_pos = {4}".format(i,j,idist,jdist,interaction_pos))
###        if i!=j:
###            exp_df.loc[i,j]=pos_square
###        #if i>10:
###            #break
###    #if j>10:
        #break
        
ExpectedFile = [ExpLine.rstrip('\n') for ExpLine in open("./chr1_expected.txt")]
Expected = []
for ThisExp in ExpectedFile:
    x = float(ThisExp)
    x = math.floor(x)
    Expected.append(x)
#Expected = range(1,498)
#Expected.reverse()
#[float(i) for i in Expected]
#print("length = {0}".format(len(Expected)))
#print(Expected)

exp_df = pd.DataFrame(columns=range(1,499),index=range(1,499))
YCounter=1
XCounter=1
while YCounter<499:
    FirstList = Expected[1:YCounter]
    SecondList = Expected[0:(499-YCounter)]
    FirstList.reverse()
    AppendList = FirstList + SecondList
    exp_df.loc[YCounter]=AppendList
    YCounter=YCounter+1

#oe_df = np.log2(df/exp_df)
#oe_df.to_csv("chr1MatrixFull_OE.csv",sep="\t")

fig, ax = plt.subplots()

ax.set_yticks((df.index.values[0::50]-1), minor=False)
Ylabels = (df.index.values[0::50]-1)/2
addString = " MB"
YlabelsText = [str(yitem)+addString for yitem in Ylabels]
ax.set_yticklabels(YlabelsText)

ax.set_xticks((df.columns.values[0::50])-1, minor=False)
Xlabels = (df.columns.values[0::50]-1)/2
XlabelsText = [str(xitem)+addString for xitem in Xlabels]
ax.set_xticklabels(XlabelsText)

#ax.set_title('In situ Hi-C, K562 chr1',fontdict={'verticalAlignment':5})
plt.title('In situ Hi-C, GEO12878 chr1, 500kb (Obs)', y=1.07)
#plt.title('In situ Hi-C, K562 chr1, 500kb (Exp)', y=1.07)
#plt.title('In situ Hi-C, K562 chr1, 500kb (O/E)', y=1.07)
#print(str(df.shape[0])+", y="+str(df.shape[1]))
ax.set_xlim(-0.5, df.shape[1]-0.5)
#ax.set_ylim(df.shape[0]-0.5, -0.5)
ax.set_ylim(-0.5, df.shape[0]-0.5)
ax.invert_yaxis()
#ax.xaxis.set_label_position('top')
ax.xaxis.set_ticks_position('top')
ax.tick_params(axis='both', direction='out', length=6)

heatmap = plt.imshow(df, cmap=plt.cm.Reds, interpolation='nearest', vmax=247)
#heatmap = plt.imshow(oe_df, cmap='seismic', interpolation='nearest', vmax=4, vmin=-4)
cbar = plt.colorbar(heatmap)
#cbar.set_label('Observed/Expected', rotation=270, labelpad=16)
cbar.set_label('Observed value', rotation=270, labelpad=16)
#cbar.set_label('Expected value', rotation=270, labelpad=16)
#plt.plot([2,50], [2,2], 'k-')

#ax.minorticks_off()
#plt.imshow(df, cmap=plt.cm.Reds, interpolation='nearest', norm=colors.LogNorm())

#print(a)
#plt.colorbar()
plt.show()
