import os
import pandas as pd
import math
from bayes_opt import BayesianOptimization
import matplotlib.pyplot as plt 
import subprocess
import warnings
import seaborn as sb
import json
from prodict import Prodict

f=open('/your/dir/Bayesian_Optimization_Example_Conf.json')

json=Prodict.from_dict(json.load(f))

strfigs=""


str1= """

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
"""

str2="""

level="""+str(json.normlevel)+"""
steps="""+str(json.steps)+"""
sample="""+str(json.sample)+"""
dirforresults='"""+str(json.resultsdir)+"""'
ignoref="""+str(json.ignorefirst)+"""
onlymean="""+str(json.onlymean)+"""
sims="""+str(json.simulations)+"""
samplelistt="""+str(json.samplelist)+"""
withsamplelist="""+str(json.withsamplelist)+"""
bestiteration=float(-10000000000000000)
ctiterations=0
os.chdir('"""+str(json.exampledir)+"""')
pd.set_option('mode.chained_assignment',None)
warnings.simplefilter(action='ignore', category=FutureWarning)

df2=[]
dfm=pd.DataFrame()
targetcsvs="""+str(json.targetcsvs)+"""
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

"""

strlgenerators2="""

def l1generator2(stepslist,steps,df2,df3):
  sum=0
  i=0
  while (i<=steps):
    if(i in stepslist):
      sum=sum+abs((df2.iloc[:,3][i])-(df3.iloc[:,3][i]))+ \\
      abs((df2.iloc[:,4][i])-(df3.iloc[:,4][i]))+ \\
      abs((df2.iloc[:,5][i])-(df3.iloc[:,5][i]))+ \\
      abs((df2.iloc[:,6][i])-(df3.iloc[:,6][i]))+ \\
      abs(((df2.iloc[:,7][steps]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][steps]+df3.iloc[:,8][i]+df3.iloc[:,9][i])))

    i=i+1
  return float(sum)

def l2generator2(stepslist,steps,df2,df3):
  sum=0
  i=0
  while (i<=steps):
    if(i in stepslist):
      sum=sum+pow(abs((df2.iloc[:,3][i])-(df3.iloc[:,3][i])),2)+ \\
      pow(abs((df2.iloc[:,4][i])-(df3.iloc[:,4][i])),2)+ \\
      pow(abs((df2.iloc[:,5][i])-(df3.iloc[:,5][i])),2)+ \\
      pow(abs((df2.iloc[:,6][i])-(df3.iloc[:,6][i])),2)+ \\
      pow(abs(((df2.iloc[:,7][steps]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][steps]+df3.iloc[:,8][i]+df3.iloc[:,9][i]))),2)
    i=i+1
  return float(sum)


"""

strlgenerators="""

def l1generator(steps,step,df2,df3,ignorefirst=False):
  sum=0
  i=0
  st=step
  if (st>steps):
    st=steps
  if(ignorefirst==True and i==0):
    i=i+st
  while (i<=steps):
    sum=sum+abs((df2.iloc[:,3][i])-(df3.iloc[:,3][i]))+ \\
    abs((df2.iloc[:,4][i])-(df3.iloc[:,4][i]))+ \\
    abs((df2.iloc[:,5][i])-(df3.iloc[:,5][i]))+ \\
    abs((df2.iloc[:,6][i])-(df3.iloc[:,6][i]))+ \\
    abs(((df2.iloc[:,7][steps]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][steps]+df3.iloc[:,8][i]+df3.iloc[:,9][i])))

    i=i+st
  if ((i-st)<steps):
    sum=sum+abs((df2.iloc[:,3][steps])-(df3.iloc[:,3][steps]))+ \\
    abs((df2.iloc[:,4][steps])-(df3.iloc[:,4][steps]))+ \\
    abs((df2.iloc[:,5][steps])-(df3.iloc[:,5][steps]))+ \\
    abs((df2.iloc[:,6][steps])-(df3.iloc[:,6][steps]))+ \\
    abs(((df2.iloc[:,7][steps]+df2.iloc[:,8][steps]+df2.iloc[:,9][steps])-(df3.iloc[:,7][steps]+df3.iloc[:,8][steps]+df3.iloc[:,9][steps])))

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
    sum=sum+pow(abs((df2.iloc[:,3][i])-(df3.iloc[:,3][i])),2)+ \\
    pow(abs((df2.iloc[:,4][i])-(df3.iloc[:,4][i])),2)+ \\
    pow(abs((df2.iloc[:,5][i])-(df3.iloc[:,5][i])),2)+ \\
    pow(abs((df2.iloc[:,6][i])-(df3.iloc[:,6][i])),2)+ \\
    pow(abs(((df2.iloc[:,7][steps]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][steps]+df3.iloc[:,8][i]+df3.iloc[:,9][i]))),2)
    i=i+st
  if ((i-st)<steps):
    sum=sum+pow(abs((df2.iloc[:,3][steps])-(df3.iloc[:,3][steps])),2)+ \\
    pow(abs((df2.iloc[:,4][steps])-(df3.iloc[:,4][steps])),2)+ \\
    pow(abs((df2.iloc[:,5][steps])-(df3.iloc[:,5][steps])),2)+ \\
    pow(abs((df2.iloc[:,6][steps])-(df3.iloc[:,6][steps])),2)+ \\
    pow(abs(((df2.iloc[:,7][steps]+df2.iloc[:,8][steps]+df2.iloc[:,9][steps])-(df3.iloc[:,7][steps]+df3.iloc[:,8][steps]+df3.iloc[:,9][steps]))),2)
  return math.sqrt(sum)

"""


strnewrow="""{'Iteration':(ctiterations+1),'Score':"""
strnewrowl2=strnewrow+"""score"""
strnewrowl1=strnewrow+"""score"""



strdf4="""

df4=pd.DataFrame({'Iteration':[0],'Score':None"""

strdf4init="""

df4.to_csv(dirforresults+'output.csv',index=False, mode="w")


"""
strfundef="""

def black_box_function("""


strbody="""
    global ctiterations
    global df4
    df1=pd.read_csv('"""+str(json.inputcsv)+"""',comment='#')
    i1=(df1[df1['parameter_name']=='number_of_time_steps'].index[0])
    df1.iloc[i1][2]=steps

"""
strbounds="""

pbounds={"""

for variable in json.variables:
    variable=Prodict.from_dict(variable)
    strdf4=strdf4+""",'"""+str(variable.name)+"""':None"""
    if(variable.static==False):
        strfundef=strfundef+"""x_"""+str(variable.name)+""","""
        strbounds=strbounds+"""'x_"""+str(variable.name)+"""': ("""+str(variable.lbound)+""","""+str(variable.ubound)+"""),"""
        if(variable.type=="float"):
            strbody=strbody+"""
    i1=(df1[df1['parameter_name']=='"""+str(variable.realname)+"""'].index[0])
    df1.iloc[i1][2]="""+"""x_"""+str(variable.name)+"""
    print(df1.iloc[i1][2])

    """
        if(variable.type=="int"):
            strbody=strbody+"""
    i1=(df1[df1['parameter_name']=='"""+str(variable.realname)+"""'].index[0])
    df1.iloc[i1][2]="""+"""x_"""+str(variable.name)+""".round(0).astype(int)
    print(df1.iloc[i1][2])

    """
        if(variable.type=="bool"):
            strbody=strbody+"""
    i1=(df1[df1['parameter_name']=='"""+str(variable.realname)+"""'].index[0])
    if(x_"""+str(variable.name)+""">=("""+str(variable.ubound)+"""/2):
      df1.iloc[i1][2]=True
    else:
      df1.iloc[i1][2]=False
    print(df1.iloc[i1][2])

    """
    if(variable.static==True):
        strbody=strbody+"""
    x_"""+str(variable.name)+"""="""+str(variable.value)+"""
    i1=(df1[df1['parameter_name']=='"""+str(variable.realname)+"""'].index[0])
    df1.iloc[i1][2]="""+"""x_"""+str(variable.name)+"""
    print(df1.iloc[i1][2])

    """
 
    strnewrowl1=strnewrowl1+""",'"""+str(variable.name)+"""':x_"""+str(variable.name)
    strnewrowl2=strnewrowl2+""",'"""+str(variable.name)+"""':x_"""+str(variable.name)
    strfigs=strfigs+"""

#plt.figure()
#plt.title('target')
#df4.plot(x='Iteration',y='"""+str(variable.name)+"""')
#plt.savefig(dirforresults+'iter_"""+str(variable.name)+"""'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#df4.plot(x="Score",y='"""+str(variable.name)+"""', kind="box")
#plt.savefig(dirforresults+'targ_"""+str(variable.name)+"""'+'.pdf',format="pdf", bbox_inches="tight")
#plt.figure()
#ax=sb.regplot(x=df5['"""+str(variable.name)+"""'],y=df5['Score'],ci=95)
#plt.savefig(dirforresults+'p95_"""+str(variable.name)+"""'+'.pdf',format="pdf", bbox_inches="tight")

"""

strdf4=strdf4+"""})"""
strfundef=strfundef[:-1]
strfundef=strfundef+"""):"""
strnewrowl2=strnewrowl2+"""}"""
strnewrowl1=strnewrowl1+"""}"""
strbounds=strbounds[:-1]
strbounds=strbounds+"""}"""
strbody=strbody+"""

    li1=[]
    for s in range(sims):
      i1=(df1[df1['parameter_name']=='output_directory'].index[0])
      df1.iloc[i1][2]="results_"+str(s)
      df1.to_csv('"""+str(json.inputcsv)+"""', index=False)
      try:
        p = subprocess.Popen("make run", start_new_session=True, shell=True, bufsize=1)
        p.wait(timeout="""+str(json.timeout)+""")
      except:
        os.killpg(os.getpgid(p.pid), signal.SIGKILL)
        return float("-10000000000000000")
    for s in range(sims):
      dfffp=pd.read_csv('"""+str(json.exampledir)+"""/results_'+str(s)+'/stats.csv',comment='#')
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
        new_row="""+str(strnewrowl1)+"""
        print(new_row)
        if (new_row['Score']>=bestiteration):
          bestiteration=new_row['Score']
          df3.to_csv(dirforresults+"bestiteration.csv",index=False)
        lendf4=len(df4)
        df4.loc[lendf4]=new_row
        print(df4)
        if lendf4>="""+str(json.saveresultsperiter)+""":
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
        new_row="""+str(strnewrowl2)+"""
        print(new_row)
        if (new_row['Score']>=bestiteration):
          bestiteration=new_row['Score']
          df3.to_csv(dirforresults+"bestiteration.csv",index=False)
        lendf4=len(df4)
        df4.loc[lendf4]=new_row
        print(df4)
        if lendf4>="""+str(json.saveresultsperiter)+""":
          df4[1:].to_csv(dirforresults+'output.csv',index=False,mode="a",header=False)
          df4=pd.DataFrame({'Iteration':[0],'Score':None,'can_apoptose_prob':None,'can_apoptose_prob_inc':None,'can_apoptose_time_window':None,'can_apoptose_time_window_to_delete':None,'can_divide_prob':None,'can_divide_prob_inc':None,'can_divide_time_window':None,'can_divide_max':None,'can_grow_prob':None,'diameter_rate':None,'can_migrate_half_range':None})
          pickle.dump(optimizer, open(dirforresults+'model.pkl', 'wb'))
        ctiterations=ctiterations+1
        return score

    """
strfinal="""

optimizer = BayesianOptimization(
    f=black_box_function,
    pbounds=pbounds,
    random_state=1,
    #allow_duplicate_points=True,
)
optimizer.maximize(
    init_points=0,
    n_iter="""+str(json.iterations)+""",
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

"""
strfinalf="""

print(optimizer.max)
pickle.dump(optimizer, open(dirforresults+'model.pkl', 'wb'))

"""
#the last lines are for test

strfinalfile=str1+str2+strdf4+strdf4init+strlgenerators+strlgenerators2+strfundef+strbody+strbounds+strfinal+strfigs+strfinalf
print(strfinalfile)

ff=open(""+str(json.filename)+"","w")
ff.write(strfinalfile)
ff.close()
f.close()
