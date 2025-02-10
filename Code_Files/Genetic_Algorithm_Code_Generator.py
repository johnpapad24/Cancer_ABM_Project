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

f=open('/your/dir/Genetic_Algorithm_Example_Conf.json')

json=Prodict.from_dict(json.load(f))

strfigs=""

varindex=0
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


import operator
import random

import numpy

from deap import base
from deap import benchmarks
from deap import creator
from deap import tools

"""

strvars="""
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
"""

strindivreg="""toolbox.register("individual", tools.initCycle, creator.Individual, ("""

strcheckmutant="""def checkmutant(mutant):
"""

strnnrow="""      nnew_row={'String':(ctiterations+1),'Generation':generations,'Score':ind.fitness.values[0]"""


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
generations=0
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
      abs(((df2.iloc[:,7][i]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][i]+df3.iloc[:,8][i]+df3.iloc[:,9][i])))

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
      pow(abs(((df2.iloc[:,7][i]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][i]+df3.iloc[:,8][i]+df3.iloc[:,9][i]))),2)
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
    abs(((df2.iloc[:,7][i]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][i]+df3.iloc[:,8][i]+df3.iloc[:,9][i])))

    i=i+st
  if ((i-st)<steps):
    sum=sum+abs((df2.iloc[:,3][steps])-(df3.iloc[:,3][steps]))+ \\
    abs((df2.iloc[:,4][steps])-(df3.iloc[:,4][steps]))+ \\
    abs((df2.iloc[:,5][steps])-(df3.iloc[:,5][steps]))+ \\
    abs((df2.iloc[:,6][steps])-(df3.iloc[:,6][steps]))+ \\
    abs(((df2.iloc[:,7][i]+df2.iloc[:,8][steps]+df2.iloc[:,9][steps])-(df3.iloc[:,7][i]+df3.iloc[:,8][steps]+df3.iloc[:,9][steps])))

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
    pow(abs(((df2.iloc[:,7][i]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][i]+df3.iloc[:,8][i]+df3.iloc[:,9][i]))),2)
    i=i+st
  if ((i-st)<steps):
    sum=sum+pow(abs((df2.iloc[:,3][steps])-(df3.iloc[:,3][steps])),2)+ \\
    pow(abs((df2.iloc[:,4][steps])-(df3.iloc[:,4][steps])),2)+ \\
    pow(abs((df2.iloc[:,5][steps])-(df3.iloc[:,5][steps])),2)+ \\
    pow(abs((df2.iloc[:,6][steps])-(df3.iloc[:,6][steps])),2)+ \\
    pow(abs(((df2.iloc[:,7][i]+df2.iloc[:,8][steps]+df2.iloc[:,9][steps])-(df3.iloc[:,7][i]+df3.iloc[:,8][steps]+df3.iloc[:,9][steps]))),2)
  return math.sqrt(sum)

"""


strnewrow="""{'String':(ctiterations+1),'Generation':generations,'Score':"""
strnewrowl2=strnewrow+"""score"""
strnewrowl1=strnewrow+"""score"""



strdf4="""

df4=pd.DataFrame({'String':[0],'Generation':None,'Score':None"""

strdf4init="""

df4.to_csv(dirforresults+'output.csv',index=False, mode="w")


"""
strfundef="""

def black_box_function(individual):

"""


strbody="""
    global ctiterations
    global df4
    df1=pd.read_csv('"""+str(json.inputcsv)+"""',comment='#')
    i1=(df1[df1['parameter_name']=='number_of_time_steps'].index[0])
    df1.iloc[i1][2]=steps

"""
strbounds="""

pbounds={"""
cii=0
for variable in json.variables:
    variable=Prodict.from_dict(variable)
    strdf4=strdf4+""",'"""+str(variable.name)+"""':None"""
    if(variable.static==False):
        
        if(variable.type=="float"):
            strvars=strvars+ """toolbox.register("""+'"x_'+str(variable.name)+'"'+""", random.uniform, """+str(variable.lbound)+""", """+str(variable.ubound)+""")
"""
            strnnrow=strnnrow+""",'"""+str(variable.name)+"""':ind["""+str(cii)+"""]"""
            strindivreg=strindivreg+"""toolbox.x_"""+str(variable.name)+","
            strcheckmutant=strcheckmutant+""" if(mutant["""+str(cii)+"""]>="""+str(variable.ubound)+"""):
  mutant["""+str(cii)+"""]="""+str(variable.ubound)+"""
 if(mutant["""+str(cii)+"""]<="""+str(variable.lbound)+"""):
  mutant["""+str(cii)+"""]="""+str(variable.lbound)+"""

"""
            strbody=strbody+"""
    i1=(df1[df1['parameter_name']=='"""+str(variable.realname)+"""'].index[0])
    df1.iloc[i1][2]=individual["""+str(cii)+"""]
    print(df1.iloc[i1][2])

    """
        if(variable.type=="int"):
          strnnrow=strnnrow+""",'"""+str(variable.name)+"""':ind["""+str(cii)+"""]"""
          strvars=strvars+ """toolbox.register("""+'"x_'+str(variable.name)+'"'+""", random.randint, """+str(variable.lbound)+""", """+str(variable.ubound)+""")
"""
          strindivreg=strindivreg+"""toolbox.x_"""+str(variable.name)+","
          strcheckmutant=strcheckmutant+""" if(mutant["""+str(cii)+"""]>="""+str(variable.ubound)+"""):
  mutant["""+str(cii)+"""]="""+str(variable.ubound)+"""
 if(mutant["""+str(cii)+"""]<="""+str(variable.lbound)+"""):
  mutant["""+str(cii)+"""]="""+str(variable.lbound)+"""
 if not isinstance(mutant["""+str(cii)+"""], int):
  mutant["""+str(cii)+"""]=int(round(mutant["""+str(cii)+"""]))

"""
          strbody=strbody+"""
    i1=(df1[df1['parameter_name']=='"""+str(variable.realname)+"""'].index[0])
    df1.iloc[i1][2]=individual["""+str(cii)+"""]
    print(df1.iloc[i1][2])
    """
        if(variable.type=="bool"):
            strbody=strbody+"""
    i1=(df1[df1['parameter_name']=='"""+str(variable.realname)+"""'].index[0])
    if(individual["""+str(cii)+"""]>=("""+str(variable.ubound)+"""/2):
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
 
    strnewrowl1=strnewrowl1+""",'"""+str(variable.name)+"""':individual["""+str(cii)+"""]"""
    strnewrowl2=strnewrowl2+""",'"""+str(variable.name)+"""':individual["""+str(cii)+"""]"""
    cii=cii+1
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
strindivreg=strindivreg[:-1]
strindivreg=strindivreg+"""))"""
strindivreg=strindivreg+"""

toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.1)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("evaluate",black_box_function)
pop=toolbox.population("""+str(json.iterationspergen)+""")
print("-- Generation %i --" % generations)
fitnesses = list(map(toolbox.evaluate, pop))
print(tuple(zip(pop, fitnesses)))
for ind, fit in zip(pop, fitnesses):
  ind.fitness.values = fit
while generations<="""+str(json.maxgenerations)+""":
  CXPB, MUTPB = """+str(json.crossoverprob)+""", """+str(json.mutationprob)+"""
  fits = [ind.fitness.values[0] for ind in pop]
  generations = generations + 1
  print("-- Generation %i --" % generations)
  offspring = toolbox.select(pop, len(pop))
  offspring = list(map(toolbox.clone, offspring))
  for child1, child2 in zip(offspring[::2], offspring[1::2]):
    if random.random() < CXPB:
      toolbox.mate(child1, child2)
      child1=checkmutant(child1)
      child2=checkmutant(child2)
      del child1.fitness.values
      del child2.fitness.values
  for mutant in offspring:
    if random.random() < MUTPB:
      toolbox.mutate(mutant)
      mutant=checkmutant(mutant)
      del mutant.fitness.values
  validind = [ind for ind in offspring if ind.fitness.valid]
  for ind in validind:
"""+strnnrow+"""}

      llendf4=len(df4)
      df4.loc[llendf4]=nnew_row
      ctiterations=ctiterations+1
      print('here    here')
      print(df4)
  invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
  fitnessesinvld = map(toolbox.evaluate, invalid_ind)
  for ind, fit in zip(invalid_ind, fitnessesinvld):
        ind.fitness.values = fit
  pop[:] = offspring
  lengthl = len(pop)
  meanl = sum(fits) / lengthl
  sum2 = sum(xf*xf for xf in fits)
  std = abs(sum2 / lengthl - meanl**2)**0.5
  print("  Min %s" % min(fits))
  print("  Max %s" % max(fits))
  print("  Avg %s" % meanl)
  print("  Std %s" % std)




print(df4)
df5=df4.copy(deep=True)
df5['Score']=df5['Score'].replace(float(-10000000000000000),np.nan)
df5['Score']=df5['Score'].replace('',np.nan)
df5=df5.dropna()
print(df5)
df4[1:].to_csv(dirforresults+'output.csv',index=False, mode="a",header=False)

"""



strdf4=strdf4+"""})"""
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
          """+strdf4.strip()+"""
      
        ctiterations=ctiterations+1
        return score,
    
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
          """+strdf4.strip()+"""
        ctiterations=ctiterations+1
        return score,

    """
strfinal="""

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
strcheckmutant=strcheckmutant+"""
 return mutant"""
#the last lines are for test

strfinalfile=str1+str2+strdf4+strdf4init+strcheckmutant+strlgenerators+strlgenerators2+strfundef+strbody+strvars+strindivreg
print(strfinalfile)

ff=open(""+str(json.filename)+"","w")
ff.write(strfinalfile)
ff.close()
f.close()
