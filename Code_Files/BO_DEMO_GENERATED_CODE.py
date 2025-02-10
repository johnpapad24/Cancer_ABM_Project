

import os 
import pandas as pd 
import math 
from bayes_opt import BayesianOptimization 
import matplotlib.pyplot as plt 
import subprocess 
import warnings 
import seaborn as sb 
import sys
import signal
import numpy as np
import pickle


level=2
steps=11
sample=8
dirforresults='/home/john/vasvav-invitro_neuro-fd03fdfd4ce4/res/'
ignoref=False
onlymean=False
sims=3
samplelistt=[0, 3, 5, 7, 11]
withsamplelist=True
bestiteration=float(-10000000000000000)
ctiterations=0
os.chdir('/home/john/vasvav-invitro_neuro-fd03fdfd4ce4/examples/HeLa_cells')
pd.set_option('mode.chained_assignment',None)
warnings.simplefilter(action='ignore', category=FutureWarning)

df2=[]
dfm=pd.DataFrame()
targetcsvs=['/home/john/restreated/t11.csv', '/home/john/restreated/t12.csv', '/home/john/restreated/t13.csv', '/home/john/restreated/t14.csv', '/home/john/restreated/t15.csv', '/home/john/restreated/t16.csv', '/home/john/restreated/t17.csv', '/home/john/restreated/t18.csv']
ii=0
for targetcsv in targetcsvs:
  ii=ii+1
  df22=pd.read_csv(targetcsv,comment='#')
  print(df22.iloc[:,0])
  for jj in range(steps+1):
    if (not df22['current_time'].isin([jj]).any().any()):
      df22=pd.concat([df22,pd.DataFrame([{"current_time":jj}])],ignore_index=False)
  df22=df22.fillna(0).astype(int)
  df22=df22.sort_values('current_time')
  df22=df22.reset_index(drop=True)
  df22.to_csv(dirforresults+"target"+str(ii)+".csv",index=False)
  df2.append(df22)
  print(df22)
  #plt.figure()
  #plt.plot(df2.iloc[:,2])
  #plt.title('target')
  #plt.savefig(dirforresults+'target'+str(ii)+'.pdf',format='pdf', bbox_inches='tight')
print(ii)
dfm=pd.concat(df2)
dfm1=pd.DataFrame({"current_time":dfm["current_time"].unique().astype(int)})
dfm2=pd.DataFrame(dfm.groupby("current_time").mean().round().astype(int))
dfma=[dfm1,dfm2]
dfm=pd.DataFrame()
dfm=pd.concat(dfma,axis=1)
dfm.to_csv(dirforresults+"target_mean"+".csv",index=False)
print(dfm)



df4=pd.DataFrame({'Iteration':[0],'Score':None,'can_apoptose_prob':None,'can_apoptose_prob_inc':None,'can_apoptose_time_window':None,'can_apoptose_time_window_to_delete':None,'can_divide_prob':None,'can_divide_prob_inc':None,'can_divide_time_window':None,'can_divide_max':None,'can_grow_prob':None,'diameter_rate':None,'can_migrate_half_range':None})

df4.to_csv(dirforresults+'output.csv',index=False, mode="w")




def l1generator(steps,step,df2,df3,ignorefirst=False):
  sum=0
  i=0
  st=step
  if (st>steps):
    st=steps
  if(ignorefirst==True and i==0):
    i=i+st
  while (i<=steps):
    sum=sum+abs((df2.iloc[:,3][i])-(df3.iloc[:,3][i]))+ \
    abs((df2.iloc[:,4][i])-(df3.iloc[:,4][i]))+ \
    abs((df2.iloc[:,5][i])-(df3.iloc[:,5][i]))+ \
    abs((df2.iloc[:,6][i])-(df3.iloc[:,6][i]))+ \
    abs((df2.iloc[:,7][i])-(df3.iloc[:,7][i]))+ \
    abs(((df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,8][i]+df3.iloc[:,9][i])))

    i=i+st
  if ((i-st)<steps):
    sum=sum+abs((df2.iloc[:,3][steps])-(df3.iloc[:,3][steps]))+ \
    abs((df2.iloc[:,4][steps])-(df3.iloc[:,4][steps]))+ \
    abs((df2.iloc[:,5][steps])-(df3.iloc[:,5][steps]))+ \
    abs((df2.iloc[:,6][steps])-(df3.iloc[:,6][steps]))+ \
    abs((df2.iloc[:,7][steps])-(df3.iloc[:,7][steps]))+ \
    abs(((df2.iloc[:,8][steps]+df2.iloc[:,9][steps])-(df3.iloc[:,8][steps]+df3.iloc[:,9][steps])))

  return float(sum)

def l2generator(steps,step,df2,df3,ignorefirst=False):
  sum=0
  i=0
  st=step
  if (st>steps):
    st=steps
  if(ignorefirst==True and i==0):
    i=i+st
  while (i<=steps):
    sum=sum+pow(abs((df2.iloc[:,3][i])-(df3.iloc[:,3][i])),2)+ \
    pow(abs((df2.iloc[:,4][i])-(df3.iloc[:,4][i])),2)+ \
    pow(abs((df2.iloc[:,5][i])-(df3.iloc[:,5][i])),2)+ \
    pow(abs((df2.iloc[:,6][i])-(df3.iloc[:,6][i])),2)+ \
    pow(abs((df2.iloc[:,7][i])-(df3.iloc[:,7][i])),2)+ \
    pow(abs(((df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,8][i]+df3.iloc[:,9][i]))),2)
    i=i+st
  if ((i-st)<steps):
    sum=sum+pow(abs((df2.iloc[:,3][steps])-(df3.iloc[:,3][steps])),2)+ \
    pow(abs((df2.iloc[:,4][steps])-(df3.iloc[:,4][steps])),2)+ \
    pow(abs((df2.iloc[:,5][steps])-(df3.iloc[:,5][steps])),2)+ \
    pow(abs((df2.iloc[:,6][steps])-(df3.iloc[:,6][steps])),2)+ \
    pow(abs((df2.iloc[:,7][steps])-(df3.iloc[:,7][steps])),2)+ \
    pow(abs(((df2.iloc[:,8][steps]+df2.iloc[:,9][steps])-(df3.iloc[:,8][steps]+df3.iloc[:,9][steps]))),2)
  return math.sqrt(sum)



def l1generator2(stepslist,steps,df2,df3):
  sum=0
  i=0
  while (i<=steps):
    if(i in stepslist):
      sum=sum+abs((df2.iloc[:,3][i])-(df3.iloc[:,3][i]))+ \
      abs((df2.iloc[:,4][i])-(df3.iloc[:,4][i]))+ \
      abs((df2.iloc[:,5][i])-(df3.iloc[:,5][i]))+ \
      abs((df2.iloc[:,6][i])-(df3.iloc[:,6][i]))+ \
      abs((df2.iloc[:,7][i])-(df3.iloc[:,7][i]))+ \
      abs(((df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,8][i]+df3.iloc[:,9][i])))

    i=i+1
  return float(sum)

def l2generator2(stepslist,steps,df2,df3):
  sum=0
  i=0
  while (i<=steps):
    if(i in stepslist):
      sum=sum+pow(abs((df2.iloc[:,3][i])-(df3.iloc[:,3][i])),2)+ \
      pow(abs((df2.iloc[:,4][i])-(df3.iloc[:,4][i])),2)+ \
      pow(abs((df2.iloc[:,5][i])-(df3.iloc[:,5][i])),2)+ \
      pow(abs((df2.iloc[:,6][i])-(df3.iloc[:,6][i])),2)+ \
      pow(abs((df2.iloc[:,7][i])-(df3.iloc[:,7][i])),2)+ \
      pow(abs(((df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,8][i]+df3.iloc[:,9][i]))),2)
    i=i+1
  return float(sum)




def black_box_function(x_can_apoptose_prob,x_can_apoptose_prob_inc,x_can_apoptose_time_window,x_can_apoptose_time_window_to_delete,x_can_divide_prob,x_can_divide_prob_inc,x_can_divide_time_window,x_can_divide_max,x_can_grow_prob,x_diameter_rate,x_can_migrate_half_range):
    global ctiterations
    global df4
    df1=pd.read_csv('/home/john/vasvav-invitro_neuro-fd03fdfd4ce4/examples/HeLa_cells/input.csv',comment='#')
    i1=(df1[df1['parameter_name']=='number_of_time_steps'].index[0])
    df1.iloc[i1][2]=steps


    i1=(df1[df1['parameter_name']=='cancer_cell/can_apoptose/probability'].index[0])
    df1.iloc[i1][2]=x_can_apoptose_prob
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_apoptose/probability_increment_with_age'].index[0])
    df1.iloc[i1][2]=x_can_apoptose_prob_inc
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_apoptose/time_window'].index[0])
    df1.iloc[i1][2]=x_can_apoptose_time_window.round(0).astype(int)
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_apoptose/time_window/to_delete'].index[0])
    df1.iloc[i1][2]=x_can_apoptose_time_window_to_delete.round(0).astype(int)
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_divide/probability'].index[0])
    df1.iloc[i1][2]=x_can_divide_prob
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_divide/probability_increment_with_age'].index[0])
    df1.iloc[i1][2]=x_can_divide_prob_inc
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_divide/time_window'].index[0])
    df1.iloc[i1][2]=x_can_divide_time_window.round(0).astype(int)
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_divide/max'].index[0])
    df1.iloc[i1][2]=x_can_divide_max.round(0).astype(int)
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_grow/probability'].index[0])
    df1.iloc[i1][2]=x_can_grow_prob
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_grow/diameter_rate'].index[0])
    df1.iloc[i1][2]=x_diameter_rate
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_migrate/half_range'].index[0])
    df1.iloc[i1][2]=x_can_migrate_half_range
    print(df1.iloc[i1][2])

    

    li1=[]
    for s in range(sims):
      i1=(df1[df1['parameter_name']=='output_directory'].index[0])
      df1.iloc[i1][2]="results_"+str(s)
      df1.to_csv('/home/john/vasvav-invitro_neuro-fd03fdfd4ce4/examples/HeLa_cells/input.csv', index=False)
      try:
        p = subprocess.Popen("make run", start_new_session=True, shell=True, bufsize=1)
        p.wait(timeout=300)
      except:
        os.killpg(os.getpgid(p.pid), signal.SIGKILL)
        return float("-10000000000000000")
    for s in range(sims):
      dfffp=pd.read_csv('/home/john/vasvav-invitro_neuro-fd03fdfd4ce4/examples/HeLa_cells/results_'+str(s)+'/stats.csv',comment='#')
      li1.append(dfffp)
    dfpf=pd.concat(li1)
    byrowindex=dfpf.groupby(dfpf.index)
    df3=byrowindex.mean().round(0).astype(int)
    title="iteration "+str(ctiterations+1)
    df3.to_csv(dirforresults+title+".csv",index=False)

    #plt.figure()
    #plt.plot(df3.iloc[:,2])
    #plt.savefig(dirforresults+title+".pdf",format="pdf", bbox_inches="tight")
    
    global bestiteration
    if(level==1):
        score=0
        if (onlymean==False):
          ii=0
          for df22 in df2:
            if (withsamplelist==False):
              score=score+l1generator(steps,sample,df22,df3,ignoref)
            else:
              score=score+l1generator2(samplelistt,steps,df22,df3)
            ii=ii+1
          score=-(score/ii)
        else:
          if (withsamplelist==False):
            score=score+l1generator(steps,sample,dfm,df3,ignoref)
          else:
            score=score+l1generator2(samplelistt,steps,dfm,df3)
        new_row={'Iteration':(ctiterations+1),'Score':score,'can_apoptose_prob':x_can_apoptose_prob,'can_apoptose_prob_inc':x_can_apoptose_prob_inc,'can_apoptose_time_window':x_can_apoptose_time_window,'can_apoptose_time_window_to_delete':x_can_apoptose_time_window_to_delete,'can_divide_prob':x_can_divide_prob,'can_divide_prob_inc':x_can_divide_prob_inc,'can_divide_time_window':x_can_divide_time_window,'can_divide_max':x_can_divide_max,'can_grow_prob':x_can_grow_prob,'diameter_rate':x_diameter_rate,'can_migrate_half_range':x_can_migrate_half_range}
        print(new_row)
        if (new_row['Score']>=bestiteration):
          bestiteration=new_row['Score']
          df3.to_csv(dirforresults+"bestiteration.csv",index=False)
        lendf4=len(df4)
        df4.loc[lendf4]=new_row
        print(df4)
        if lendf4>=200:
          df4[1:].to_csv(dirforresults+'output.csv',index=False,mode="a",header=False)
          df4=pd.DataFrame({'Iteration':[0],'Score':None,'can_apoptose_prob':None,'can_apoptose_prob_inc':None,'can_apoptose_time_window':None,'can_apoptose_time_window_to_delete':None,'can_divide_prob':None,'can_divide_prob_inc':None,'can_divide_time_window':None,'can_divide_max':None,'can_grow_prob':None,'diameter_rate':None,'can_migrate_half_range':None})
          pickle.dump(optimizer, open(dirforresults+'model.pkl', 'wb'))
        ctiterations=ctiterations+1
        return score
    
    if(level==2):
        score=0
        if (onlymean==False):
          ii=0
          for df22 in df2:
            if (withsamplelist==False):
              score=score+l2generator(steps,sample,df22,df3,ignoref)
            else:
              score=score+l2generator2(samplelistt,steps,df22,df3)
            ii=ii+1
          score=-(score/ii)
        else:
          if (withsamplelist==False):
            score=score+l2generator(steps,sample,dfm,df3,ignoref)
          else:
            score=score+l2generator2(samplelistt,steps,dfm,df3)
        new_row={'Iteration':(ctiterations+1),'Score':score,'can_apoptose_prob':x_can_apoptose_prob,'can_apoptose_prob_inc':x_can_apoptose_prob_inc,'can_apoptose_time_window':x_can_apoptose_time_window,'can_apoptose_time_window_to_delete':x_can_apoptose_time_window_to_delete,'can_divide_prob':x_can_divide_prob,'can_divide_prob_inc':x_can_divide_prob_inc,'can_divide_time_window':x_can_divide_time_window,'can_divide_max':x_can_divide_max,'can_grow_prob':x_can_grow_prob,'diameter_rate':x_diameter_rate,'can_migrate_half_range':x_can_migrate_half_range}
        print(new_row)
        if (new_row['Score']>=bestiteration):
          bestiteration=new_row['Score']
          df3.to_csv(dirforresults+"bestiteration.csv",index=False)
        lendf4=len(df4)
        df4.loc[lendf4]=new_row
        print(df4)
        if lendf4>=200:
          df4[1:].to_csv(dirforresults+'output.csv',index=False,mode="a",header=False)
          df4=pd.DataFrame({'Iteration':[0],'Score':None,'can_apoptose_prob':None,'can_apoptose_prob_inc':None,'can_apoptose_time_window':None,'can_apoptose_time_window_to_delete':None,'can_divide_prob':None,'can_divide_prob_inc':None,'can_divide_time_window':None,'can_divide_max':None,'can_grow_prob':None,'diameter_rate':None,'can_migrate_half_range':None})
          pickle.dump(optimizer, open(dirforresults+'model.pkl', 'wb'))
        ctiterations=ctiterations+1
        return score

    

pbounds={'x_can_apoptose_prob': (0.1,0.95),'x_can_apoptose_prob_inc': (0.0001,0.1),'x_can_apoptose_time_window': (0,11),'x_can_apoptose_time_window_to_delete': (0,11),'x_can_divide_prob': (0.1,0.95),'x_can_divide_prob_inc': (0.0001,0.1),'x_can_divide_time_window': (0,11),'x_can_divide_max': (0,5),'x_can_grow_prob': (0.1,0.95),'x_diameter_rate': (0.5,19.5),'x_can_migrate_half_range': (0.0,20.0)}

optimizer = BayesianOptimization(
    f=black_box_function,
    pbounds=pbounds,
    random_state=1,
    #allow_duplicate_points=True,
)
optimizer.maximize(
    init_points=0,
    n_iter=6400,
)
print(df4)
df5=df4.copy(deep=True)
df5['Score']=df5['Score'].replace(float(-10000000000000000),np.nan)
df5['Score']=df5['Score'].replace('',np.nan)
df5=df5.dropna()
print(df5)
df4[1:].to_csv(dirforresults+'output.csv',index=False, mode="a",header=False)
#plt.figure()
#df4.plot(x="Iteration",y="Score")
#plt.savefig(dirforresults+"itertarg"+".pdf",format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_apoptose_prob')
#plt.savefig(dirforresults+'iter_can_apoptose_prob'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_apoptose_prob', kind="box")
#plt.savefig(dirforresults+'targ_can_apoptose_prob'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_apoptose_prob'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_apoptose_prob'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_apoptose_prob_inc')
#plt.savefig(dirforresults+'iter_can_apoptose_prob_inc'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_apoptose_prob_inc', kind="box")
#plt.savefig(dirforresults+'targ_can_apoptose_prob_inc'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_apoptose_prob_inc'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_apoptose_prob_inc'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_apoptose_time_window')
#plt.savefig(dirforresults+'iter_can_apoptose_time_window'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_apoptose_time_window', kind="box")
#plt.savefig(dirforresults+'targ_can_apoptose_time_window'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_apoptose_time_window'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_apoptose_time_window'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_apoptose_time_window_to_delete')
#plt.savefig(dirforresults+'iter_can_apoptose_time_window_to_delete'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_apoptose_time_window_to_delete', kind="box")
#plt.savefig(dirforresults+'targ_can_apoptose_time_window_to_delete'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_apoptose_time_window_to_delete'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_apoptose_time_window_to_delete'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_divide_prob')
#plt.savefig(dirforresults+'iter_can_divide_prob'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_divide_prob', kind="box")
#plt.savefig(dirforresults+'targ_can_divide_prob'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_divide_prob'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_divide_prob'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_divide_prob_inc')
#plt.savefig(dirforresults+'iter_can_divide_prob_inc'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_divide_prob_inc', kind="box")
#plt.savefig(dirforresults+'targ_can_divide_prob_inc'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_divide_prob_inc'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_divide_prob_inc'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_divide_time_window')
#plt.savefig(dirforresults+'iter_can_divide_time_window'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_divide_time_window', kind="box")
#plt.savefig(dirforresults+'targ_can_divide_time_window'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_divide_time_window'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_divide_time_window'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_divide_max')
#plt.savefig(dirforresults+'iter_can_divide_max'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_divide_max', kind="box")
#plt.savefig(dirforresults+'targ_can_divide_max'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_divide_max'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_divide_max'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_grow_prob')
#plt.savefig(dirforresults+'iter_can_grow_prob'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_grow_prob', kind="box")
#plt.savefig(dirforresults+'targ_can_grow_prob'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_grow_prob'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_grow_prob'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='diameter_rate')
#plt.savefig(dirforresults+'iter_diameter_rate'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='diameter_rate', kind="box")
#plt.savefig(dirforresults+'targ_diameter_rate'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['diameter_rate'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_diameter_rate'+'.pdf',format="pdf", bbox_inches="tight")



#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='can_migrate_half_range')
#plt.savefig(dirforresults+'iter_can_migrate_half_range'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='can_migrate_half_range', kind="box")
#plt.savefig(dirforresults+'targ_can_migrate_half_range'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['can_migrate_half_range'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_can_migrate_half_range'+'.pdf',format="pdf", bbox_inches="tight")



print(optimizer.max)
pickle.dump(optimizer, open(dirforresults+'model.pkl', 'wb'))

