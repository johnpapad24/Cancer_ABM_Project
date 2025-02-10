import os 
import pandas as pd 
import math 
from bayes_opt import BayesianOptimization 
from statsmodels.nonparametric.smoothers_lowess import lowess as  sm_lowess
import matplotlib as mplb
import matplotlib.pyplot as plt 
import subprocess 
import warnings 
import seaborn as sb 
import sys
import signal
import numpy as np
import scipy as sp
import glob
from pathlib import Path
import re
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerBase

class HandlerColormap(HandlerBase):
    def __init__(self, cmap, num_stripes=8, **kw):
        HandlerBase.__init__(self, **kw)
        self.cmap = cmap
        self.num_stripes = num_stripes
    def create_artists(self, legend, orig_handle, 
                       xdescent, ydescent, width, height, fontsize, trans):
        stripes = []
        for i in range(self.num_stripes):
            s = Rectangle([xdescent + i * width / self.num_stripes, ydescent], 
                          width / self.num_stripes, 
                          height, 
                          fc=self.cmap((2 * i + 1) / (2 * self.num_stripes)), 
                          transform=trans)
            stripes.append(s)
        return stripes

def translate(string):
    string=string.strip()
    dictt={"N_cells_pheno_1":"Number of cells",
          "N_cells_pheno_1_Ap":"Number of apoptosed cells",
          "N_cells_pheno_1_G1":"Number of cells in G1 phase",
          "N_cells_pheno_1_Sy":"Number of cells in S phase",
          "N_cells_pheno_1_G2":"Number of cells in G2 phase",
          }
    return dictt[string]
    
def movingaverage(data,windowsize):
    i=0
    smoothed=[]
    while i < len(data) - windowsize + 1:
        window_average = round(np.sum(data[i:i+windowsize]) / windowsize, 2)
        smoothed.append(window_average)
        i=i+1
    return smoothed

def smoothTriangle(data, degree):
    triangle=np.concatenate((np.arange(degree + 1), np.arange(degree)[::-1])) # up then down
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        smoothed.append(np.sum(point)/np.sum(triangle))
    # Handle boundaries
    smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed

def unique(lst):
    unique_list=[]
    for x in lst:
        if x not in unique_list:
            unique_list.append(x)
    return(unique_list)     
   
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def targetbestgraphpoints3(graphdir,points,wow=False,semilog=False):
    df4444=pd.read_csv(graphdir+'target_mean.csv',comment='#')
    df3333=pd.read_csv(graphdir+'bestiteration.csv',comment='#')
    df5555=pd.read_csv(graphdir+"output.csv")
    df5555=df5555.sort_values(by="Score",ascending=True)
    if(wow==False):
        worstiteration=df5555.iloc[0,0]
        print(worstiteration)
    medianiteradion=math.ceil(np.median(df5555.iloc[:,0]))
    print(medianiteradion)
    df7777=pd.read_csv(graphdir+"iteration "+str(medianiteradion)+".csv",comment='#')
    legendok=False
    if(wow==False):
        df6666=pd.read_csv(graphdir+"iteration "+str(worstiteration)+".csv",comment='#')
    filenames=glob.glob(graphdir+"target*.csv")
    plt.figure()
    for point in df4444.iloc[:,0]:
            if point in points:
               plt.scatter(x=point,y=df4444.iloc[point,2], label="Target mean",color="orangered")
               plt.scatter(x=point,y=df3333.iloc[point,2], label="Best Iteration",color="green", marker="d")
               if(wow==False):
                plt.scatter(x=point,y=df6666.iloc[point,2], label="Worst Iteration",color="red", marker="X")
               plt.scatter(x=point,y=df7777.iloc[point,2], label="Median Iteration",color="cyan", marker="o")
               if(semilog==True):
                 print("yscale")



    for filename in filenames:
        if "mean" in filename:
            continue
        df1111=pd.read_csv(filename,comment='#')
        for point in df1111.iloc[:,0]:
            if point in points:
                plt.scatter(x=point,y=df1111.iloc[point,2],color="blue",alpha=1/(len(filenames)-1))
                if(semilog==True):
                    print("yscale")


                if legendok==False:
                     plt.legend()
                legendok=True


    plt.ylabel("Number of cells")
    plt.xlabel("Time")
    if(wow==False):
        if(semilog==False):
            plt.savefig(graphdir+'targetbestwithpoints3.pdf',format='pdf', bbox_inches='tight')
        else:
            plt.savefig(graphdir+'targetbestwithpoints3semilog.pdf',format='pdf', bbox_inches='tight')
    else:
        if(semilog==False):
            plt.savefig(graphdir+'targetbestwithpoints3wow.pdf',format='pdf', bbox_inches='tight')
        else:
            plt.savefig(graphdir+'targetbestwithpoints3wowsemilog.pdf',format='pdf', bbox_inches='tight')
        





graphdir="/home/john/res/"
#graphs=" filtergenetic variablescore3genetic scorehistgenetic scorehistgenetictr variablescore3geneticnosmooth"
#graphs=" filterbaysian variablescore3new scorehist variablescore3newnosmooth"
graphs="targetbestgraphpoints3allwownew"


points=[0,3,7,11]
steps=11
step=1
graphstok=graphs.split(" ")

if "targetmeangraph" in graphstok:
    df2=pd.read_csv(graphdir+'target_mean.csv',comment='#')
    plt.figure()
    plt.plot(df2.iloc[:,2])
    plt.ylabel("Number of cells")
    plt.xlabel("Time")
    plt.title('tarßget')
    plt.savefig(graphdir+'target_mean.pdf',format='pdf', bbox_inches='tight')

if "scorehist" in graphstok:
    df8022=pd.read_csv(graphdir+"output.csv")
    plt.figure()
    plt.title("score histogram")
    #bins=
    a=np.array(df8022.iloc[:,1])
    a=a[~np.isnan(a)]
    upper_lim=np.percentile(a,90)
    lower_lim=np.percentile(a,10)
    iqr=(upper_lim-lower_lim)*1.5
    qset=(lower_lim-iqr,upper_lim+iqr)
    plist=[]
    for y in a.tolist():
        if y>=qset[0] and y<=qset[1]:
            plist.append(y)
    plt.hist(plist, bins=100)
    plt.axvline(np.median(plist),color="red", linestyle="dashed",linewidth="1", label="median score")
    plt.legend()
    plt.savefig(graphdir+"scorehist"+".pdf",format='pdf', bbox_inches='tight')

if "filterbaysian" in graphstok:
    df1001u=pd.read_csv(graphdir+"output.csv")
    df1001umax=df1001u.sort_values("Score",ascending=False)
    df1001umax=df1001umax.head(1)    
    df1001umax.to_csv(graphdir+'oldbest.csv',index=False, mode="w",header=True)

if "filtergenetic" in graphstok:
    df1001u=pd.read_csv(graphdir+"output.csv")
    df1001f=df1001u[df1001u["can_apoptose_prob"]>=0.1]
    df1001f=df1001f[df1001f["can_apoptose_prob"]<=0.95]
    df1001f=df1001f[df1001f["can_apoptose_prob_inc"]>=0.0001] 
    df1001f=df1001f[df1001f["can_apoptose_prob_inc"]<=0.1]
    df1001f=df1001f[df1001f["can_apoptose_time_window"]>=0]
    df1001f=df1001f[df1001f["can_apoptose_time_window"]<=11]
    df1001f=df1001f[df1001f["can_apoptose_time_window_to_delete"]>=0]
    df1001f=df1001f[df1001f["can_apoptose_time_window_to_delete"]<=11]
    df1001f=df1001f[df1001f["can_divide_prob"]>=0.1]
    df1001f=df1001f[df1001f["can_divide_prob"]<=0.95]
    df1001f=df1001f[df1001f["can_divide_prob_inc"]>=0.0001] 
    df1001f=df1001f[df1001f["can_divide_prob_inc"]<=0.1]
    df1001f=df1001f[df1001f["can_divide_time_window"]>=0]
    df1001f=df1001f[df1001f["can_divide_time_window"]<=11]
    df1001f=df1001f[df1001f["can_divide_max"]>=0]
    df1001f=df1001f[df1001f["can_divide_max"]<=5]
    df1001f=df1001f[df1001f["can_grow_prob"]>=0.1]
    df1001f=df1001f[df1001f["can_grow_prob"]<=0.95]
    df1001f=df1001f[df1001f["diameter_rate"]>=0.5] 
    df1001f=df1001f[df1001f["diameter_rate"]<=19.5]
    df1001f=df1001f[df1001f["can_migrate_half_range"]>=0]
    df1001f=df1001f[df1001f["can_migrate_half_range"]<=19.5]
    df1001f.to_csv(graphdir+'newoutput.csv',index=False, mode="w",header=True)
    df1001fmax=df1001f.sort_values("Score",ascending=False)
    df1001fmax=df1001fmax.head(1)    
    df1001fmax.to_csv(graphdir+'newbest.csv',index=False, mode="w",header=True)
    df1001umax=df1001u.sort_values("Score",ascending=False)
    df1001umax=df1001umax.head(1)    
    df1001umax.to_csv(graphdir+'oldbest.csv',index=False, mode="w",header=True)


if "scorehistgenetic" in graphstok:
    df8022=pd.read_csv(graphdir+"output.csv")
    plt.figure()
    plt.title("score histogram")
    #bins=
    a=np.array(df8022.iloc[:,2])
    a=a[~np.isnan(a)]
    upper_lim=np.percentile(a,90)
    lower_lim=np.percentile(a,10)
    iqr=(upper_lim-lower_lim)*1.5
    qset=(lower_lim-iqr,upper_lim+iqr)
    plist=[]
    for y in a.tolist():
        if y>=qset[0] and y<=qset[1]:
            plist.append(y)
    plt.hist(plist, bins=100)
    plt.axvline(np.median(plist),color="red", linestyle="dashed",linewidth="1", label="median score")
    plt.legend()
    plt.savefig(graphdir+"scorehist"+".pdf",format='pdf', bbox_inches='tight')

if "scorehistgenetictr" in graphstok:
    df8022=pd.read_csv(graphdir+"newoutput.csv")
    plt.figure()
    plt.title("score histogram")
    #bins=
    a=np.array(df8022.iloc[:,2])
    a=a[~np.isnan(a)]
    upper_lim=np.percentile(a,90)
    lower_lim=np.percentile(a,10)
    iqr=(upper_lim-lower_lim)*1.5
    qset=(lower_lim-iqr,upper_lim+iqr)
    plist=[]
    for y in a.tolist():
        if y>=qset[0] and y<=qset[1]:
            plist.append(y)
    plt.hist(plist, bins=100)
    plt.axvline(np.median(plist),color="red", linestyle="dashed",linewidth="1", label="median score")
    plt.legend()
    plt.savefig(graphdir+"scorehisttr"+".pdf",format='pdf', bbox_inches='tight')


    
if "scorehistsimple" in graphstok:
    df8022=pd.read_csv(graphdir+"output.csv")
    plt.figure()
    plt.title("score histogram")
    #bins=
    a=np.array(df8022.iloc[:,1])
    a=a[~np.isnan(a)]
    plist=[]
    plt.hist(a, bins=1000)
    plt.axvline(np.median(a),color="red", linestyle="dashed",linewidth="1", label="median score")
    plt.legend()
    plt.savefig(graphdir+"scorehistsimple"+".pdf",format='pdf', bbox_inches='tight')


if "targetgraph" in graphstok:
    filenames=glob.glob(graphdir+"target*.csv")
    plt.figure()
    iii=1.25
    legendok=False
    n_lines=len(filenames)-1
    c = np.arange(1, n_lines + 1)
    norm = mplb.colors.Normalize(vmin=c.min(), vmax=c.max())
    cmap = mplb.cm.ScalarMappable(norm=norm, cmap=mplb.cm.Blues)     
    cmap.set_array([])
    for filename in filenames:
        print(filename)
        df111=pd.read_csv(filename,comment='#')
        if "mean" in filename:
            for point in df111.iloc[:,0]:
                if point in points:     
                 plt.scatter(x=point,y=df111.iloc[point,2],color="orangered", label="Target mean")
                 if legendok==False:
                     plt.legend()
                     legendok=True
        else:
            for point in df111.iloc[:,0]:
                if point in points:
                 plt.scatter(x=point,y=df111.iloc[point,2],color="blue",alpha=1/len(filenames))
                 
            iii=iii+1/n_lines
    plt.ylabel("Number of cells")
    plt.xlabel("Time")
    plt.savefig(graphdir+"targetsgraph"+'.pdf',format='pdf', bbox_inches='tight')

if "targetgraphall" in graphstok:
    filenames=glob.glob(graphdir+"target*.csv")
    plt.figure()
    iii=1.25
    legendok=[False,False,False,False,False,False]
    #n_lines=len(filenames)-1
    #c = np.arange(1, n_lines + 1)
   ##
    # norm = mplb.colors.Normalize(vmin=c.min(), vmax=c.max())
    #cmap = mplb.cm.ScalarMappable(norm=norm, cmap=mplb.cm.Blues)     
    #cmap.set_array([])
    fig,axs=plt.subplots(6,sharex=True,sharey=True)
    plt.subplots_adjust(hspace=0.5,top=2)
    fig.supylabel("Number of cells")
    fig.supxlabel("Time")
    a=0
    for filename in filenames:
        df111=pd.read_csv(filename,comment='#')
        if "mean" in filename:
            for point in df111.iloc[:,0]:
                a=0
                if point in points:
                 for j in range(2,9):
                  if j==3:
                   continue
                  axs[a].set_title(str(df111.columns[j]))
                  axs[a].scatter(x=point,y=df111.iloc[point,j],color="orangered", label="Target mean")
                  if legendok[a]==False:
                    axs[a].legend()
                    legendok[a]=True
                  a=a+1
                            
        else:
            for point in df111.iloc[:,0]:
                a=0
                if point in points:
                 for j in range(2,9):
                  if j==3:
                   continue
                  axs[a].set_title(str(df111.columns[j]))
                  axs[a].scatter(x=point,y=df111.iloc[point,j],color="blue",alpha=1/len(filenames))
                  a=a+1
                 
            #iii=iii+1/n_lines
    plt.savefig(graphdir+"targetsgraphall"+'.pdf',format='pdf', bbox_inches='tight')


if "iterationsgraph" in graphstok:
    for filename in glob.glob(graphdir+"iteration*.csv"):
        filn=Path(filename).stem
        df1=pd.read_csv(filename,comment='#')
        print(filn)
        plt.figure()
        plt.plot(df1.iloc[:,2])
        plt.title(filn.capitalize())
        plt.ylabel("Number of cells")
        plt.xlabel("Time")
        plt.savefig(graphdir+filn+'.pdf',format='pdf', bbox_inches='tight')

if "targetbestgraph" in graphstok:
    df4=pd.read_csv(graphdir+'target_mean.csv',comment='#')
    df3=pd.read_csv(graphdir+'bestiteration.csv',comment='#')

    plt.figure()
    plt.plot(df4.iloc[:,2])
    plt.plot(df3.iloc[:,2])
    plt.title("Comparison of mean target and best iteration")
    plt.ylabel("Number of cells")
    plt.xlabel("Time")
    plt.savefig(graphdir+'targetbest.pdf',format='pdf', bbox_inches='tight')

if "targetbestgraphpoints" in graphstok:
    df444=pd.read_csv(graphdir+'target_mean.csv',comment='#')
    df333=pd.read_csv(graphdir+'bestiteration.csv',comment='#')


    plt.figure()
    plt.plot(df444.iloc[:,2],label="Target mean" , c="blue")
    plt.plot(df333.iloc[:,2],'-*', markevery=points, mfc='r', mec="r", label="Best iteration with points", c="green")
    plt.title("Comparison of mean target and best iteration")
    plt.legend()
    plt.ylabel("Number of cells")
    plt.xlabel("Time")
    plt.savefig(graphdir+'targetbestwithpoints.pdf',format='pdf', bbox_inches='tight')

if "targetbestgraphpoints3" in graphstok:
    targetbestgraphpoints3(graphdir,points)

if "targetbestgraphpoints3wow" in graphstok:
    targetbestgraphpoints3(graphdir,points,wow=True)

if "targetbestgraphpoints3semilog" in graphstok:
    targetbestgraphpoints3(graphdir,points,semilog=True)

if "targetbestgraphpoints3wowsemilog" in graphstok:
    targetbestgraphpoints3(graphdir,points,wow=True,semilog=True)

if "targetbestgraphpoints3allwow" in graphstok:
    df44444=pd.read_csv(graphdir+'target_mean.csv',comment='#')
    df33333=pd.read_csv(graphdir+'bestiteration.csv',comment='#')
    df55555=pd.read_csv(graphdir+"output.csv")
    df55555=df55555.sort_values(by="Score",ascending=True)
    legendok=[False,False,False,False,False,False]
    plt.figure()
    figg,axxs=plt.subplots(6,sharex=True,sharey=True)
    plt.subplots_adjust(hspace=0.5,top=5)
    figg.supylabel("Number of cells")
    figg.supxlabel("Time")
    a=0
    ll=0
    medianiteradion=math.ceil(np.median(df55555.iloc[:,0]))
    print(medianiteradion)
    df77777=pd.read_csv(graphdir+"iteration "+str(medianiteradion)+".csv",comment='#')
    filenames=glob.glob(graphdir+"target*.csv")

    for filename in filenames:
        if "mean" in filename:
            continue
        df11111=pd.read_csv(filename,comment='#')
        for point in df11111.iloc[:,0]:
            a=0
            if point in points:
                for j in range(2,9):
                    maxll=len(range(2,9))
                    if j==3:
                     continue
                    axxs[a].set_title(str(df11111.columns[j]))
                    if ll<maxll-1: 
                        axxs[a].scatter(x=point,y=df11111.iloc[point,j],color="blue",label="Iteration",alpha=1/(len(filenames)-1))
                    else:
                        axxs[a].scatter(x=point,y=df11111.iloc[point,j],color="blue",alpha=1/(len(filenames)-1))
                    a=a+1
                    ll=ll+1


    for point in df44444.iloc[:,0]:
            a=0
            if point in points:
               for j in range(2,9):
                    if j==3:
                     continue
                    axxs[a].set_title(str(df44444.columns[j]))
                    axxs[a].scatter(x=point,y=df44444.iloc[point,j], label="Target mean",color="orangered")
                    axxs[a].scatter(x=point,y=df33333.iloc[point,j], label="Best Iteration",color="green", marker="d")
                    axxs[a].scatter(x=point,y=df77777.iloc[point,j], label="Median Iteration",color="cyan", marker="p")
                    if legendok[a]==False:
                        #axxs[a].legend()
                        print("test")
                    legendok[a]=True
                    a=a+1

    lines_labels=[ax.get_legend_handles_labels() for ax in axxs]
    lines,labels=[sum(lol,[]) for lol in zip(*lines_labels)]
    lines=unique(lines)
    labels=unique(labels)
    #plt.legend(lines,labels,loc='lower center', ncol=4,fancybox=True,shadow=True)
    plt.figlegend(lines,labels,loc='outside upper center', ncol=4,fancybox=True,shadow=True,bbox_to_anchor=(0.5,5.2,0,0))
    plt.savefig(graphdir+'targetbestwithpoints3allwow.pdf',format='pdf', bbox_inches='tight')


if "targetbestgraphpoints3all" in graphstok:
    df44444=pd.read_csv(graphdir+'target_mean.csv',comment='#')
    df33333=pd.read_csv(graphdir+'bestiteration.csv',comment='#')
    df55555=pd.read_csv(graphdir+"output.csv")
    df55555=df55555.sort_values(by="Score",ascending=True)
    legendok=[False,False,False,False,False,False]
    plt.figure()
    figg,axxs=plt.subplots(6,sharex=True,sharey=True)
    plt.subplots_adjust(hspace=0.5,top=5)
    figg.supylabel("Number of cells")
    figg.supxlabel("Time")
    a=0
    ll=0
    worstiteration=df55555.iloc[0,0]
    print(worstiteration)
    medianiteradion=math.ceil(np.median(df55555.iloc[:,0]))
    print(medianiteradion)
    df66666=pd.read_csv(graphdir+"iteration "+str(worstiteration)+".csv",comment='#')
    df77777=pd.read_csv(graphdir+"iteration "+str(medianiteradion)+".csv",comment='#')
    filenames=glob.glob(graphdir+"target*.csv")
    for filename in filenames:
        if "mean" in filename:
            continue
        df11111=pd.read_csv(filename,comment='#')
        for point in df11111.iloc[:,0]:
            a=0
            if point in points:
                for j in range(2,9):
                    maxll=len(range(2,9))
                    if j==3:
                     continue
                    axxs[a].set_title(str(df11111.columns[j]))
                    if ll<maxll-1:
                         axxs[a].scatter(x=point,y=df11111.iloc[point,j],color="blue",label="Target", alpha=1/(len(filenames)-1))
                    else:
                         axxs[a].scatter(x=point,y=df11111.iloc[point,j],color="blue", alpha=1/(len(filenames)-1))
                    a=a+1
                    ll=ll+1


    for point in df44444.iloc[:,0]:
            a=0
            if point in points:
               for j in range(2,9):
                    if j==3:
                     continue
                    axxs[a].set_title(str(df44444.columns[j]))
                    axxs[a].scatter(x=point,y=df44444.iloc[point,j], label="Target mean",color="orangered")
                    axxs[a].scatter(x=point,y=df33333.iloc[point,j], label="Best Iteration",color="green", marker="d")
                    axxs[a].scatter(x=point,y=df66666.iloc[point,j], label="Worst Iteration",color="red", marker="X")
                    axxs[a].scatter(x=point,y=df77777.iloc[point,j], label="Median Iteration",color="cyan", marker="o")
                    if legendok[a]==False:
                        axxs[a].legend()
                    legendok[a]=True
                    a=a+1


    plt.savefig(graphdir+'targetbestwithpoints3all.pdf',format='pdf', bbox_inches='tight')

if "targetbestgraphpoints3allsemilog" in graphstok:
    df4444444=pd.read_csv(graphdir+'target_mean.csv',comment='#')
    df3333333=pd.read_csv(graphdir+'bestiteration.csv',comment='#')
    df5555555=pd.read_csv(graphdir+"output.csv")
    df5555555=df5555555.sort_values(by="Score",ascending=True)
    legendok=[False,False,False,False,False,False]
    plt.figure()
    figg,axxs=plt.subplots(6,sharex=True,sharey=True)
    plt.subplots_adjust(hspace=0.5,top=5)
    figg.supylabel("Number of cells")
    figg.supxlabel("Time")
    a=0
    ll=0
    worstiteration=df5555555.iloc[0,0]
    print(worstiteration)
    medianiteradion=math.ceil(np.median(df5555555.iloc[:,0]))
    print(medianiteradion)
    df6666666=pd.read_csv(graphdir+"iteration "+str(worstiteration)+".csv",comment='#')
    df7777777=pd.read_csv(graphdir+"iteration "+str(medianiteradion)+".csv",comment='#')
    filenames=glob.glob(graphdir+"target*.csv")
    for filename in filenames:
        if "mean" in filename:
            continue
        df1111111=pd.read_csv(filename,comment='#')
        for point in df1111111.iloc[:,0]:
            a=0
            if point in points:
                for j in range(2,9):
                    maxll=len(range(2,9))
                    if j==3:
                     continue
                    axxs[a].set_title(str(df1111111.columns[j]))
                    if ll<maxll-1:
                         axxs[a].scatter(x=point,y=df1111111.iloc[point,j],color="blue",label="Target", alpha=1/(len(filenames)-1))
                    else:
                         axxs[a].scatter(x=point,y=df1111111.iloc[point,j],color="blue", alpha=1/(len(filenames)-1))
                    #axxs[a].set_yscale("log")
                    a=a+1
                    ll=ll+1


    for point in df4444444.iloc[:,0]:
            a=0
            if point in points:
               for j in range(2,9):
                    if j==3:
                     continue
                    axxs[a].set_title(str(df4444444.columns[j]))
                    axxs[a].scatter(x=point,y=df4444444.iloc[point,j], label="Target mean",color="orangered")
                    axxs[a].scatter(x=point,y=df3333333.iloc[point,j], label="Best Iteration",color="green", marker="d")
                    axxs[a].scatter(x=point,y=df6666666.iloc[point,j], label="Worst Iteration",color="red", marker="X")
                    axxs[a].scatter(x=point,y=df7777777.iloc[point,j], label="Median Iteration",color="cyan", marker="o")
                    if legendok[a]==False:
                        axxs[a].legend()
                    legendok[a]=True
                    axxs[a].set_yscale("log")
                    a=a+1



    plt.savefig(graphdir+'targetbestwithpoints3allsemilog.pdf',format='pdf', bbox_inches='tight')


if "iterationsgraph2" in graphstok:
    filenames=glob.glob(graphdir+"iteration*.csv")
    totalgraphs=len(filenames)+1
    print (totalgraphs)
    i=1
    z=1
    while i*i<totalgraphs:
        i=i+1
    filenames.sort(key=natural_keys)
    df5=pd.read_csv(graphdir+'target_mean.csv',comment='#')
    plt.figure()

    fig,ax=plt.subplots(i,i, figsize=(i+i+8,i+i+8))
    fig.tight_layout(pad=7.0)
    fig.supylabel("Number of cells",fontsize=20)
    fig.supxlabel("Time", fontsize=20)
    plt.subplot(i,i,z)
    plt.plot(df5.iloc[:,2])
    
    plt.title('Target')
    z=z+1
    for filename in filenames:
        filn=Path(filename).stem
        df6=pd.read_csv(filename,comment='#')
        print(filn)
        plt.subplot(i,i,z)
        plt.plot(df6.iloc[:,2])
        plt.gca().set_title(filn.capitalize())
        z=z+1
    plt.savefig(graphdir+'alliterations.pdf',format='pdf', bbox_inches='tight')

if "variableiteration" in graphstok:
    df7=pd.read_csv(graphdir+"output.csv")
    for col in df7.columns[2:]:
        print(col)
        plt.figure()
        plt.title(col)
        df7.plot(x='Iteration',y=col,ylabel="Value",xlabel="Iteration")
        plt.savefig(graphdir+'iter_'+col+'.pdf',format="pdf", bbox_inches="tight")


if "variableiteration2" in graphstok:
    df8=pd.read_csv(graphdir+"output.csv")
    totalgraphs=len(df8.columns)-2
    i=1
    z=1
    while i*i<totalgraphs:
        i=i+1
    plt.figure()
    fig,ax=plt.subplots(i,i, figsize=(i+i+8,i+i+8))
    fig.tight_layout(pad=7.0)
    fig.supylabel("Value",fontsize=20)
    fig.supxlabel("Iteration", fontsize=20)
    for col in df8.columns[2:]:
        print(col)
        plt.subplot(i,i,z)
        plt.title(col)
        plt.plot(df8["Iteration"],df8[col])
        z=z+1
    plt.savefig(graphdir+'allvars_iter'+'.pdf',format="pdf", bbox_inches="tight")

if "variablescore" in graphstok:
    df9=pd.read_csv(graphdir+"output.csv")
    for col in df9.columns[2:]:
        print(col)
        plt.figure()
        plt.title(col)
        df9.plot(x='Score',y=col,ylabel="Value",xlabel="Score",kind="box")
        plt.savefig(graphdir+'targ_'+col+'.pdf',format="pdf", bbox_inches="tight")


if "variablescore2" in graphstok:
    df10=pd.read_csv(graphdir+"output.csv")
    totalgraphs=len(df10.columns)-2
    i=1
    z=1
    while i*i<totalgraphs:
        i=i+1
    plt.figure()
    fig,ax=plt.subplots(i,i, figsize=(i+i+8,i+i+8))
    fig.tight_layout(pad=7.0)
    fig.supylabel("Value",fontsize=20)
    fig.supxlabel("Score", fontsize=20)
    for col in df10.columns[2:]:
        print(col)
        plt.subplot(i,i,z)
        plt.title(col)
        plt.scatter(df10["Score"],df10[col])
        z=z+1
    plt.savefig(graphdir+'allvars_targ'+'.pdf',format="pdf", bbox_inches="tight")

if "variablescorereg" in graphstok:
    df11=pd.read_csv(graphdir+"output.csv")
    for col in df11.columns[2:]:
        print(col)
        plt.figure()
        plt.title(col)
        ax=sb.regplot(x=df11[col],y=df11['Score'],ci=95)
        plt.savefig(graphdir+'targreg_'+col+'.pdf',format="pdf", bbox_inches="tight")


if "variablescorereg2" in graphstok:
    df12=pd.read_csv(graphdir+"output.csv")
    totalgraphs=len(df12.columns)-2
    i=1
    z=1
    while i*i<totalgraphs:
        i=i+1
    plt.figure()
    fig,ax=plt.subplots(i,i, figsize=(i+i+8,i+i+8))
    fig.tight_layout(pad=7.0)
    fig.supylabel("Score",fontsize=20)
    fig.supxlabel("Value", fontsize=20)
    for col in df12.columns[2:]:
        print(col)
        plt.subplot(i,i,z)
        plt.title(col)
        ax=sb.regplot(x=df12[col],y=df12['Score'],ci=95).set(xlabel=None,ylabel=None)
        z=z+1
    plt.savefig(graphdir+'allvars_targreg'+'.pdf',format="pdf", bbox_inches="tight")

if "variablescore3" in graphstok:
    df1001=pd.read_csv(graphdir+"output.csv")
    totalgraphs=len(df1001.columns)-2
    i=0
    plt.figure()
    ffig,aaxs=plt.subplots(totalgraphs)
    plt.subplots_adjust(hspace=0.8,top=3)
    df1001=df1001.dropna()
    df1001max=df1001.sort_values("Score",ascending=False)
    df1001max=df1001max.head(1)
    for col in df1001.columns[2:]:
        xx=df1001[col]
        yy=df1001["Score"]
        yymax=df1001max["Score"]
        xxmax=df1001max[col]
        y_s=yy
        x_s=xx
        y_sm = sm_lowess(y_s,x_s, frac=1./5., it=5,
                     return_sorted = False)
        
        xgrid = np.linspace(xx.min(),xx.max(),num=6400)
        y_grid = sp.interpolate.interp1d(x_s, y_sm, 
                                        fill_value='extrapolate')(xgrid)
        K = 1
        smooths = np.stack([y_grid for k in range(K)]).T
        
        print(len(smooths))
        mean = np.nanmean(smooths, axis=1)
        print(len(mean))
        print("mean")
        print(np.nanmean(smooths,axis=1))
        print("mean+std")
        print(np.nanstd(smooths))
        print("mean-std")
        print(np.nanstd(smooths))
        aaxs[i].fill_between(xgrid, np.nanmean(smooths,axis=1)+50*np.nanstd(smooths),
                     np.nanmean(smooths,axis=1)-50*np.nanstd(smooths), alpha=0.30)
        aaxs[i].set_title(col)
        aaxs[i].plot(xgrid, mean, color='tomato', label='')
        aaxs[i].plot(xxmax,yymax, label="Best Iteration",color="blue", marker="d")
        aaxs[i].scatter(x=xx, y=yy,color="green", alpha=0.25,s=2,label="Iteration")
        aaxs[i].legend()
        i=i+1
    
    plt.savefig(graphdir+'allvars_targ_vertical'+'.pdf',format="pdf", bbox_inches="tight")




if "variablescore3geneticnosmooth" in graphstok:
    df1001=pd.read_csv(graphdir+"newoutput.csv")
    totalgraphs=len(df1001.columns)-3
    i=0
    cmaps=[mplb.colormaps['viridis_r']]
    plt.figure()
    ffig,aaxs=plt.subplots(totalgraphs)
    plt.subplots_adjust(hspace=0.8,top=5)
    df1001=df1001.dropna()
    df1001max=df1001.sort_values("Score",ascending=False)
    df1001max=df1001max.head(1)
    for col in df1001.columns[3:]:
        xx=df1001[col]
        yy=df1001["Score"]
        zz=df1001["Generation"]
        yymax=df1001max["Score"]
        xxmax=df1001max[col]
        y_s=yy
        x_s=xx
        y_sm = sm_lowess(y_s,x_s, frac=1./5., it=5,
                     return_sorted = False)
        
        xgrid = np.linspace(xx.min(),xx.max(),num=len(df1001.index))
        y_grid = sp.interpolate.interp1d(x_s, y_sm, 
                                        fill_value='extrapolate')(xgrid)
        K = 1
        smooths = np.stack([y_grid for k in range(K)]).T
        mean2 = np.nanmean(smooths, axis=1)
        deviations=[]
        iii=0
        for mm in yy:
            devit=math.sqrt((mm -mean2[iii])**2)
            iii=iii+1
            deviations.append(devit)
        
        print(deviations)
        #deviations=smoothTriangle(deviations,55)
        #print(deviations)
        #deviations=sp.interpolate.interp1d(deviations,xgrid)
        #print(len(mean))
        print("mean")
        #print(np.nanmean(smooth
        print("mean+std")
        #print(np.nanstd(newarr,axis=1))
        print("mean-std")
        #print(np.nanstd(np.array(xx,yy),axis=1))
        aaxs[i].fill_between(xgrid, mean2+deviations,
                     mean2-deviations, alpha=0.30)
        aaxs[i].set_title(col)
        aaxs[i].plot(xgrid, mean2, color='tomato', label='Mean Iteration')
        aaxs[i].plot(xxmax,yymax, label="Best Iteration",color="blue", marker="d")
        aaxs[i].scatter(x=xx, y=yy, alpha=0.1,s=2,c=zz,cmap='viridis_r')
        aaxs[i].set_ylim(-15000.0,0)
        i=i+1
    

    lines_labels=[ax.get_legend_handles_labels() for ax in aaxs]
    lines,labels=[sum(lol,[]) for lol in zip(*lines_labels)]
    lines=unique(lines)
    labels=unique(labels)
    cmap_handles = [Rectangle((0, 0), 1, 1) for _ in cmaps]
    cmap_labels=["Iteration"]
    handler_map = dict(zip(cmap_handles, 
                       [HandlerColormap(cm, num_stripes=1000) for cm in cmaps]))
    #plt.legend(lines,labels,loc='lower center', ncol=4,fancybox=True,shadow=True)
    plt.figlegend(handles=cmap_handles,labels=cmap_labels, handler_map=handler_map,loc='outside upper left', ncol=4,fancybox=True,shadow=True,bbox_to_anchor=(0,5.2,0,0))
    plt.figlegend(handles=lines,labels=labels, loc='outside upper center', ncol=4,fancybox=True,shadow=True,bbox_to_anchor=(0.5,5.2,0,0))

    
    plt.savefig(graphdir+'allvars_targ_vertical_no_smooth'+'.pdf',format="pdf", bbox_inches="tight")

if "variablescore3genetic" in graphstok:
    df1001=pd.read_csv(graphdir+"newoutput.csv")
    totalgraphs=len(df1001.columns)-3
    i=0
    cmaps=[mplb.colormaps['viridis_r']]
    plt.figure()
    ffig,aaxs=plt.subplots(totalgraphs)
    plt.subplots_adjust(hspace=0.8,top=5)
    df1001=df1001.dropna()
    df1001max=df1001.sort_values("Score",ascending=False)
    df1001max=df1001max.head(1)
    for col in df1001.columns[3:]:
        xx=df1001[col]
        yy=df1001["Score"]
        zz=df1001["Generation"]
        yymax=df1001max["Score"]
        xxmax=df1001max[col]
        y_s=yy
        x_s=xx
        y_sm = sm_lowess(y_s,x_s, frac=1./5., it=5,
                     return_sorted = False)
        
        xgrid = np.linspace(xx.min(),xx.max(),num=len(df1001.index))
        y_grid = sp.interpolate.interp1d(x_s, y_sm, 
                                        fill_value='extrapolate')(xgrid)
        K = 1
        smooths = np.stack([y_grid for k in range(K)]).T
        mean2 = np.nanmean(smooths, axis=1)
        deviations=[]
        iii=0
        for mm in yy:
            devit=math.sqrt((mm -mean2[iii])**2)
            iii=iii+1
            deviations.append(devit)
        
        print(deviations)
        deviations=smoothTriangle(deviations,55)
        #print(deviations)
        #deviations=sp.interpolate.interp1d(deviations,xgrid)
        #print(len(mean))
        print("mean")
        #print(np.nanmean(smooth
        print("mean+std")
        #print(np.nanstd(newarr,axis=1))
        print("mean-std")
        #print(np.nanstd(np.array(xx,yy),axis=1))
        aaxs[i].fill_between(xgrid, mean2+deviations,
                     mean2-deviations, alpha=0.30)
        aaxs[i].set_title(col)
        aaxs[i].plot(xgrid, mean2, color='tomato', label='Mean Iteration')
        aaxs[i].plot(xxmax,yymax, label="Best Iteration",color="blue", marker="d")
        aaxs[i].scatter(x=xx, y=yy, alpha=0.1,s=2,c=zz,cmap='viridis_r')
        aaxs[i].set_ylim(-15000.0,0)
        i=i+1
    

    lines_labels=[ax.get_legend_handles_labels() for ax in aaxs]
    lines,labels=[sum(lol,[]) for lol in zip(*lines_labels)]
    lines=unique(lines)
    labels=unique(labels)
    cmap_handles = [Rectangle((0, 0), 1, 1) for _ in cmaps]
    cmap_labels=["Iteration"]
    handler_map = dict(zip(cmap_handles, 
                       [HandlerColormap(cm, num_stripes=1000) for cm in cmaps]))
    #plt.legend(lines,labels,loc='lower center', ncol=4,fancybox=True,shadow=True)
    plt.figlegend(handles=cmap_handles,labels=cmap_labels, handler_map=handler_map,loc='outside upper left', ncol=4,fancybox=True,shadow=True,bbox_to_anchor=(0,5.2,0,0))
    plt.figlegend(handles=lines,labels=labels, loc='outside upper center', ncol=4,fancybox=True,shadow=True,bbox_to_anchor=(0.5,5.2,0,0))

    
    plt.savefig(graphdir+'allvars_targ_vertical'+'.pdf',format="pdf", bbox_inches="tight")

if "variablescore3new" in graphstok:
    df1001=pd.read_csv(graphdir+"output.csv")
    totalgraphs=len(df1001.columns)-2
    i=0
    plt.figure()
    ffig,aaxs=plt.subplots(totalgraphs)
    plt.subplots_adjust(hspace=0.8,top=5)
    df1001=df1001.dropna()
    df1001max=df1001.sort_values("Score",ascending=False)
    df1001max=df1001max.head(1)
    for col in df1001.columns[2:]:
        xx=df1001[col]
        yy=df1001["Score"]
        yymax=df1001max["Score"]
        xxmax=df1001max[col]
        y_s=yy
        x_s=xx
        y_sm = sm_lowess(y_s,x_s, frac=1./5., it=5,
                     return_sorted = False)
        
        xgrid = np.linspace(xx.min(),xx.max(),num=len(df1001.index))
        y_grid = sp.interpolate.interp1d(x_s, y_sm, 
                                        fill_value='extrapolate')(xgrid)
        K = 1
        smooths = np.stack([y_grid for k in range(K)]).T
        
        print(len(smooths))
        mean2 = np.nanmean(smooths, axis=1)
        deviations=[]
        iii=0
        for mm in yy:
            devit=math.sqrt((mm -mean2[iii])**2)
            iii=iii+1
            deviations.append(devit)
        
        deviations=smoothTriangle(deviations,55)
        aaxs[i].fill_between(xgrid, mean2+deviations,
                     mean2-deviations, alpha=0.30)
        aaxs[i].set_title(col)
        aaxs[i].plot(xgrid, mean2, color='tomato', label='')
        aaxs[i].plot(xxmax,yymax, label="Best Iteration",color="blue", marker="d")
        aaxs[i].scatter(x=xx, y=yy,color="green", alpha=0.25,s=2,label="Iteration")
        aaxs[i].set_ylim(-15000.0,0)
        i=i+1
    
    lines_labels=[ax.get_legend_handles_labels() for ax in aaxs]
    lines,labels=[sum(lol,[]) for lol in zip(*lines_labels)]
    lines=unique(lines)
    labels=unique(labels)
    plt.figlegend(handles=lines,labels=labels, loc='outside upper center', ncol=4,fancybox=True,shadow=True,bbox_to_anchor=(0.5,5.2,0,0))
    plt.savefig(graphdir+'allvars_targ_vertical'+'.pdf',format="pdf", bbox_inches="tight")

if "variablescore3newnosmooth" in graphstok:
    df1001=pd.read_csv(graphdir+"output.csv")
    totalgraphs=len(df1001.columns)-2
    i=0
    plt.figure()
    ffig,aaxs=plt.subplots(totalgraphs)
    plt.subplots_adjust(hspace=0.8,top=5)
    df1001=df1001.dropna()
    df1001max=df1001.sort_values("Score",ascending=False)
    df1001max=df1001max.head(1)
    for col in df1001.columns[2:]:
        xx=df1001[col]
        yy=df1001["Score"]
        yymax=df1001max["Score"]
        xxmax=df1001max[col]
        y_s=yy
        x_s=xx
        y_sm = sm_lowess(y_s,x_s, frac=1./5., it=5,
                     return_sorted = False)
        
        xgrid = np.linspace(xx.min(),xx.max(),num=len(df1001.index))
        y_grid = sp.interpolate.interp1d(x_s, y_sm, 
                                        fill_value='extrapolate')(xgrid)
        K = 1
        smooths = np.stack([y_grid for k in range(K)]).T
        
        print(len(smooths))
        mean2 = np.nanmean(smooths, axis=1)
        deviations=[]
        iii=0
        for mm in yy:
            devit=math.sqrt((mm -mean2[iii])**2)
            iii=iii+1
            deviations.append(devit)
        
        #deviations=smoothTriangle(deviations,55)
        aaxs[i].fill_between(xgrid, mean2+deviations,
                     mean2-deviations, alpha=0.30)
        aaxs[i].set_title(col)
        aaxs[i].plot(xgrid, mean2, color='tomato', label='')
        aaxs[i].plot(xxmax,yymax, label="Best Iteration",color="blue", marker="d")
        aaxs[i].scatter(x=xx, y=yy,color="green", alpha=0.25,s=2,label="Iteration")
        aaxs[i].set_ylim(-15000.0,0)
        i=i+1
    
    lines_labels=[ax.get_legend_handles_labels() for ax in aaxs]
    lines,labels=[sum(lol,[]) for lol in zip(*lines_labels)]
    lines=unique(lines)
    labels=unique(labels)
    plt.figlegend(handles=lines,labels=labels, loc='outside upper center', ncol=4,fancybox=True,shadow=True,bbox_to_anchor=(0.5,5.2,0,0))
    plt.savefig(graphdir+'allvars_targ_vertical_no_smooth'+'.pdf',format="pdf", bbox_inches="tight")



if "targetbestgraphpoints3allwownew" in graphstok:
    df44444=pd.read_csv(graphdir+'target_mean.csv',comment='#')
    df33333=pd.read_csv(graphdir+'l1bo.csv',comment='#')
    df55555=pd.read_csv(graphdir+'l2bo.csv',comment='#')
    df66666=pd.read_csv(graphdir+'l1gen-cut.csv',comment='#')
    df77777=pd.read_csv(graphdir+'l2gen-cut.csv',comment='#')

    legendok=[False,False,False,False,False,False]
    plt.figure()
    figg,axxs=plt.subplots(3,sharex=True,sharey=True)
    plt.subplots_adjust(hspace=0.5,top=5)
    figg.supylabel("Number of cells")
    figg.supxlabel("Time")
    a=0
    ll=0
    filenames=glob.glob(graphdir+"target*.csv")

    for j in range(4,8):
        if j==5:
            continue
        axxs[a].plot(df33333.iloc[:,0],df33333.iloc[:,j], label="Best Iteration B.O. L1",color="green")
        axxs[a].plot(df55555.iloc[:,0],df55555.iloc[:,j], label="Best Iteration B.O. L2",color="cyan")
        axxs[a].plot(df66666.iloc[:,0],df66666.iloc[:,j], label="Best Iteration G.A. L1",color="maroon")
        axxs[a].plot(df77777.iloc[:,0],df77777.iloc[:,j], label="Best Iteration G.A. L2",color="violet")
        a=a+1
    a=0
    
    filenames=glob.glob(graphdir+"target*.csv")
    for filename in filenames:
        if "mean" in filename:
            continue
        df11111=pd.read_csv(filename,comment='#')
        for point in df11111.iloc[:,0]:
            a=0
            if point in points:
                for j in range(4,8):
                    maxll=len(range(4,8))
                    if j==5:
                     continue
                    axxs[a].set_title(translate(str(df11111.columns[j])))
                    if ll<maxll-1:
                         axxs[a].scatter(x=point,y=df11111.iloc[point,j],color="blue",label="Target", alpha=1/(len(filenames)-1))
                    else:
                         axxs[a].scatter(x=point,y=df11111.iloc[point,j],color="blue", alpha=1/(len(filenames)-1))
                    a=a+1
                    ll=ll+1
        for point in df44444.iloc[:,0]:
         a=0
         if point in points:
            for j in range(4,8):
                if j==5:
                    continue
                axxs[a].set_title(translate(str(df44444.columns[j])))
                axxs[a].scatter(x=point,y=df44444.iloc[point,j], label="Target mean",color="orangered")
                if legendok[a]==False:
                        #axxs[a].legend()
                    print("test")
                legendok[a]=True
                a=a+1



    lines_labels=[ax.get_legend_handles_labels() for ax in axxs]
    lines,labels=[sum(lol,[]) for lol in zip(*lines_labels)]
    lines=unique(lines)
    labels=unique(labels)
    #plt.legend(lines,labels,loc='lower center', ncol=4,fancybox=True,shadow=True)
    plt.figlegend(lines,labels,loc='outside upper center', ncol=4,fancybox=True,shadow=True,bbox_to_anchor=(0.5,5.2,0,0))
    plt.savefig(graphdir+'targetbestwithpoints3allwow.pdf',format='pdf', bbox_inches='tight')

