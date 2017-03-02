#!/usr/bin/env python
import numpy as np
import pandas as pd
import math
import numpy.random
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import getopt,sys
from scipy.interpolate import griddata
import scipy.stats as stats
import matplotlib.mlab as mlab

def usage():
    print("usage: Proximity_SD.py -m <mode (individual or grouped)> -l <lower bound> -u <upper bound>")

mode = ""
lowPass = 0
upPass = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], 'm:l:u:h',)
except getopt.GetoptError:
    usage()
    sys.exit(2)
    
if not opts:
    usage()
    sys.exit(2)
else:
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt == '-m':
            mode = arg
        elif opt == '-l':
            lowPass = arg
        elif opt == '-u':
            upPass = arg
        else:
            usage()
            sys.exit(2)

## Import chr1 matrix to datafram (Rao et al. K562)

df = pd.read_table("chr1MatrixFull.csv",sep="\t",index_col=0)

def Individual(lower,upper):

    ## Set up dataframe for storing graphical numbers
    
    rnames = np.arange(-5,5.5,0.5)
    StartCoors = np.arange(lower,upper+1,1)
    cnames = [str(i) for i in StartCoors]
    graph_df = pd.DataFrame(columns=cnames,index=rnames)
    
    ## Calculate reads varying with distance from intersection
    
    Counter = lower
    for i in StartCoors:
        xstart = i
        ystart = xstart-10
        yloop = ystart
        ValueList = []
        while yloop<=xstart+10:
            ValueList.append(df.iloc[xstart,yloop])
            yloop+=1
        
        ## Normalise graph to maximum value (typcially centre value)
        ValueList = ValueList/max(ValueList)
        
        graph_df[str(Counter)] = ValueList
        Counter+=1
        
    print(graph_df)
    
    ## Define graph/axis parameters
    
    fig, ax = plt.subplots()
    #xticks = np.arange(-5,6,1)
    xticks = np.arange(0,21,2)
    addString = "MB"
    #xticknames = [str(xitem)+addString for xitem in xticks]
    #ax.set_xlim(-5, 5)
    ax.set_ylim(0, graph_df.values.max())
    ax.set_xticks(xticks, minor=False)
    xarray = np.arange(-5,6,1)
    xticknames = [str(xitem)+addString for xitem in xarray]
    ax.set_xticklabels(xticknames)
    
    ## Plot graph
    
    colours = ['b','r','g','y']
    ColCounter = 0
    for i in graph_df:
        plt.plot(graph_df[i], colours[ColCounter], lw=2, label="Bin"+i)
        ColCounter = 0 if ColCounter==3 else ColCounter+1
    #plt.legend()
    plt.show()
    
def Grouped():
    y_list = []
    error_list = []
    cnames = np.arange(-15,15.5,0.5)
    cnames_str = [str(i) for i in cnames]
    graph_df = pd.DataFrame(columns=cnames_str, index=range(40,242))
    
    Resolution = 500000
    
    for j in graph_df:
        ValueList = []
        DistanceAway = int(abs(float(j))*1000000)
        for i in graph_df.index.values:
            xstart = i
            xcheck = DistanceAway/Resolution
            if j[0]=="-":
                Pix = (df.iloc[xstart,xstart]-df.iloc[xstart,xstart+xcheck])/df.iloc[xstart,xstart]
            else:
                Pix = (df.iloc[xstart,xstart]-df.iloc[xstart,xstart-xcheck])/df.iloc[xstart,xstart]   
            ValueList.append(Pix)
        graph_df[j] = ValueList
        #n, x, _ = plt.hist(graph_df[j], bins=25, normed=1)
        #density = stats.gaussian_kde(ValueList)
        #plt.plot(x, density(x))
        
        #plt.hist(graph_df[j], 25, label=j+"MB")
        #plt.legend()
        y_list.append(graph_df[j].mean())
        error_list.append(graph_df[j].std())
        #print("dist = {0}, sd = {1}, mean = {2}".format(j,graph_df[j].std(),graph_df[j].mean()))
    rev_list = [1-i for i in y_list]
    add_list = [x+y for x,y in zip(rev_list,error_list)]
    min_list = [x-y for x,y in zip(rev_list,error_list)]
    fig, ax = plt.subplots()
    #ax.set_xlim(graph_df.columns.values.min(), graph_df.columns.values.max())
    #ax.set_xlim(-2, 2)
    ax.set_ylim(0, 1)
    #plt.errorbar(graph_df.columns.values, rev_list, error_list, fmt='o', marker='^')
    plt.plot(graph_df.columns.values, rev_list, 'b')
    float_values = [float(i) for i in graph_df.columns.values]
    ax.set_xlim(min(float_values), max(float_values))
    ax.fill_between(float_values, add_list, min_list, facecolor='blue', alpha=0.5)
    ax.set_title("Read drop-off varying with distance from viewpoint,\nmeasured at positions 40-242 on chr1 (500kb resolution)")
    ax.set_xlabel("Distance from viewpoint (MB)")
    ax.set_ylabel("Fraction of reads at viewpoint")
    plt.show()
    
if mode=="individual":    
    Individual(int(lowPass),int(upPass))
elif mode=="grouped":
    Grouped()
else:
    print("Error: mode option not recognised, please choose between 'individual' or 'grouped'")
    sys.exit(2)

#fig, ax = plt.subplots()
#
#ax.set_yticks((df.index.values[0::50]-1), minor=False)
#Ylabels = (df.index.values[0::50]-1)/2
#addString = " MB"
#YlabelsText = [str(yitem)+addString for yitem in Ylabels]
#ax.set_yticklabels(YlabelsText)
#
#plt.title('In situ Hi-C, K562 chr1, 500kb (Obs)', y=1.07)
#ax.set_xlim(-0.5, df.shape[1]-0.5)
#ax.set_ylim(-0.5, df.shape[0]-0.5)
#ax.invert_yaxis()
#ax.xaxis.set_ticks_position('top')
#ax.tick_params(axis='both', direction='out', length=6)
#
#heatmap = plt.imshow(df, cmap=plt.cm.Reds, interpolation='nearest', vmax=247)
#cbar = plt.colorbar(heatmap)
#cbar.set_label('Observed value', rotation=270, labelpad=16)
#print(df[:10])
#print(df.iloc[1,1])
#plt.show()

## Test code



