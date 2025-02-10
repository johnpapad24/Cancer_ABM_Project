

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



level=2
steps=11
sample=8
dirforresults='/home/john/vasvav-invitro_neuro-fd03fdfd4ce4/res1/'
ignoref=False
onlymean=False
sims=3
samplelistt=[0, 3, 5, 7, 11]
withsamplelist=True
bestiteration=float(-10000000000000000)
ctiterations=0
generations=0
os.chdir('/home/john/vasvav-invitro_neuro-fd03fdfd4ce4/examples/HeLa_cells')
pd.set_option('mode.chained_assignment',None)
warnings.simplefilter(action='ignore', category=FutureWarning)

df2=[]
dfm=pd.DataFrame()
targetcsvs=['/home/john/restreated/t01.csv', '/home/john/restreated/t02.csv', '/home/john/restreated/t03.csv', '/home/john/restreated/t04.csv', '/home/john/restreated/t05.csv', '/home/john/restreated/t06.csv', '/home/john/restreated/t07.csv', '/home/john/restreated/t08.csv']
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



df4=pd.DataFrame({'String':[0],'Generation':None,'Score':None,'can_apoptose_prob':None,'can_apoptose_prob_inc':None,'can_apoptose_time_window':None,'can_apoptose_time_window_to_delete':None,'can_divide_prob':None,'can_divide_prob_inc':None,'can_divide_time_window':None,'can_divide_max':None,'can_grow_prob':None,'diameter_rate':None})

df4.to_csv(dirforresults+'output.csv',index=False, mode="w")


def checkmutant(mutant):
 if(mutant[0]>=0.95):
  mutant[0]=0.95
 if(mutant[0]<=0.1):
  mutant[0]=0.1

 if(mutant[1]>=0.1):
  mutant[1]=0.1
 if(mutant[1]<=0.0001):
  mutant[1]=0.0001

 if(mutant[2]>=11):
  mutant[2]=11
 if(mutant[2]<=0):
  mutant[2]=0
 if not isinstance(mutant[2], int):
  mutant[2]=int(round(mutant[2]))

 if(mutant[3]>=11):
  mutant[3]=11
 if(mutant[3]<=0):
  mutant[3]=0
 if not isinstance(mutant[3], int):
  mutant[3]=int(round(mutant[3]))

 if(mutant[4]>=0.95):
  mutant[4]=0.95
 if(mutant[4]<=0.1):
  mutant[4]=0.1

 if(mutant[5]>=0.1):
  mutant[5]=0.1
 if(mutant[5]<=0.0001):
  mutant[5]=0.0001

 if(mutant[6]>=11):
  mutant[6]=11
 if(mutant[6]<=0):
  mutant[6]=0
 if not isinstance(mutant[6], int):
  mutant[6]=int(round(mutant[6]))

 if(mutant[7]>=5):
  mutant[7]=5
 if(mutant[7]<=0):
  mutant[7]=0
 if not isinstance(mutant[7], int):
  mutant[7]=int(round(mutant[7]))

 if(mutant[8]>=0.95):
  mutant[8]=0.95
 if(mutant[8]<=0.1):
  mutant[8]=0.1

 if(mutant[9]>=19.5):
  mutant[9]=19.5
 if(mutant[9]<=0.5):
  mutant[9]=0.5


 return mutant

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
    abs(((df2.iloc[:,7][i]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][i]+df3.iloc[:,8][i]+df3.iloc[:,9][i])))

    i=i+st
  if ((i-st)<steps):
    sum=sum+abs((df2.iloc[:,3][steps])-(df3.iloc[:,3][steps]))+ \
    abs((df2.iloc[:,4][steps])-(df3.iloc[:,4][steps]))+ \
    abs((df2.iloc[:,5][steps])-(df3.iloc[:,5][steps]))+ \
    abs((df2.iloc[:,6][steps])-(df3.iloc[:,6][steps]))+ \
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
    sum=sum+pow(abs((df2.iloc[:,3][i])-(df3.iloc[:,3][i])),2)+ \
    pow(abs((df2.iloc[:,4][i])-(df3.iloc[:,4][i])),2)+ \
    pow(abs((df2.iloc[:,5][i])-(df3.iloc[:,5][i])),2)+ \
    pow(abs((df2.iloc[:,6][i])-(df3.iloc[:,6][i])),2)+ \
    pow(abs(((df2.iloc[:,7][i]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][i]+df3.iloc[:,8][i]+df3.iloc[:,9][i]))),2)
    i=i+st
  if ((i-st)<steps):
    sum=sum+pow(abs((df2.iloc[:,3][steps])-(df3.iloc[:,3][steps])),2)+ \
    pow(abs((df2.iloc[:,4][steps])-(df3.iloc[:,4][steps])),2)+ \
    pow(abs((df2.iloc[:,5][steps])-(df3.iloc[:,5][steps])),2)+ \
    pow(abs((df2.iloc[:,6][steps])-(df3.iloc[:,6][steps])),2)+ \
    pow(abs(((df2.iloc[:,7][i]+df2.iloc[:,8][steps]+df2.iloc[:,9][steps])-(df3.iloc[:,7][i]+df3.iloc[:,8][steps]+df3.iloc[:,9][steps]))),2)
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
      abs(((df2.iloc[:,7][i]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][i]+df3.iloc[:,8][i]+df3.iloc[:,9][i])))

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
      pow(abs(((df2.iloc[:,7][i]+df2.iloc[:,8][i]+df2.iloc[:,9][i])-(df3.iloc[:,7][i]+df3.iloc[:,8][i]+df3.iloc[:,9][i]))),2)
    i=i+1
  return float(sum)




def black_box_function(individual):


    global ctiterations
    global df4
    df1=pd.read_csv('/home/john/vasvav-invitro_neuro-fd03fdfd4ce4/examples/HeLa_cells/input.csv',comment='#')
    i1=(df1[df1['parameter_name']=='number_of_time_steps'].index[0])
    df1.iloc[i1][2]=steps


    i1=(df1[df1['parameter_name']=='cancer_cell/can_apoptose/probability'].index[0])
    df1.iloc[i1][2]=individual[0]
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_apoptose/probability_increment_with_age'].index[0])
    df1.iloc[i1][2]=individual[1]
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_apoptose/time_window'].index[0])
    df1.iloc[i1][2]=individual[2]
    print(df1.iloc[i1][2])
    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_apoptose/time_window/to_delete'].index[0])
    df1.iloc[i1][2]=individual[3]
    print(df1.iloc[i1][2])
    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_divide/probability'].index[0])
    df1.iloc[i1][2]=individual[4]
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_divide/probability_increment_with_age'].index[0])
    df1.iloc[i1][2]=individual[5]
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_divide/time_window'].index[0])
    df1.iloc[i1][2]=individual[6]
    print(df1.iloc[i1][2])
    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_divide/max'].index[0])
    df1.iloc[i1][2]=individual[7]
    print(df1.iloc[i1][2])
    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_grow/probability'].index[0])
    df1.iloc[i1][2]=individual[8]
    print(df1.iloc[i1][2])

    
    i1=(df1[df1['parameter_name']=='cancer_cell/can_grow/diameter_rate'].index[0])
    df1.iloc[i1][2]=individual[9]
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
        new_row={'String':(ctiterations+1),'Generation':generations,'Score':score,'can_apoptose_prob':individual[0],'can_apoptose_prob_inc':individual[1],'can_apoptose_time_window':individual[2],'can_apoptose_time_window_to_delete':individual[3],'can_divide_prob':individual[4],'can_divide_prob_inc':individual[5],'can_divide_time_window':individual[6],'can_divide_max':individual[7],'can_grow_prob':individual[8],'diameter_rate':individual[9]}
        print(new_row)
        if (new_row['Score']>=bestiteration):
          bestiteration=new_row['Score']
          df3.to_csv(dirforresults+"bestiteration.csv",index=False)
        lendf4=len(df4)
        df4.loc[lendf4]=new_row
        print(df4)
        if lendf4>=200:
          df4[1:].to_csv(dirforresults+'output.csv',index=False,mode="a",header=False)
          df4=pd.DataFrame({'String':[0],'Generation':None,'Score':None,'can_apoptose_prob':None,'can_apoptose_prob_inc':None,'can_apoptose_time_window':None,'can_apoptose_time_window_to_delete':None,'can_divide_prob':None,'can_divide_prob_inc':None,'can_divide_time_window':None,'can_divide_max':None,'can_grow_prob':None,'diameter_rate':None})
      
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
        new_row={'String':(ctiterations+1),'Generation':generations,'Score':score,'can_apoptose_prob':individual[0],'can_apoptose_prob_inc':individual[1],'can_apoptose_time_window':individual[2],'can_apoptose_time_window_to_delete':individual[3],'can_divide_prob':individual[4],'can_divide_prob_inc':individual[5],'can_divide_time_window':individual[6],'can_divide_max':individual[7],'can_grow_prob':individual[8],'diameter_rate':individual[9]}
        print(new_row)
        if (new_row['Score']>=bestiteration):
          bestiteration=new_row['Score']
          df3.to_csv(dirforresults+"bestiteration.csv",index=False)
        lendf4=len(df4)
        df4.loc[lendf4]=new_row
        print(df4)
        if lendf4>=200:
          df4[1:].to_csv(dirforresults+'output.csv',index=False,mode="a",header=False)
          df4=pd.DataFrame({'String':[0],'Generation':None,'Score':None,'can_apoptose_prob':None,'can_apoptose_prob_inc':None,'can_apoptose_time_window':None,'can_apoptose_time_window_to_delete':None,'can_divide_prob':None,'can_divide_prob_inc':None,'can_divide_time_window':None,'can_divide_max':None,'can_grow_prob':None,'diameter_rate':None})
        ctiterations=ctiterations+1
        return score,

    
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("x_can_apoptose_prob", random.uniform, 0.1, 0.95)
toolbox.register("x_can_apoptose_prob_inc", random.uniform, 0.0001, 0.1)
toolbox.register("x_can_apoptose_time_window", random.randint, 0, 11)
toolbox.register("x_can_apoptose_time_window_to_delete", random.randint, 0, 11)
toolbox.register("x_can_divide_prob", random.uniform, 0.1, 0.95)
toolbox.register("x_can_divide_prob_inc", random.uniform, 0.0001, 0.1)
toolbox.register("x_can_divide_time_window", random.randint, 0, 11)
toolbox.register("x_can_divide_max", random.randint, 0, 5)
toolbox.register("x_can_grow_prob", random.uniform, 0.1, 0.95)
toolbox.register("x_diameter_rate", random.uniform, 0.5, 19.5)
toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.x_can_apoptose_prob,toolbox.x_can_apoptose_prob_inc,toolbox.x_can_apoptose_time_window,toolbox.x_can_apoptose_time_window_to_delete,toolbox.x_can_divide_prob,toolbox.x_can_divide_prob_inc,toolbox.x_can_divide_time_window,toolbox.x_can_divide_max,toolbox.x_can_grow_prob,toolbox.x_diameter_rate))

toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.1)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("evaluate",black_box_function)
pop=toolbox.population(250)
print("-- Generation %i --" % generations)
fitnesses = list(map(toolbox.evaluate, pop))
print(tuple(zip(pop, fitnesses)))
for ind, fit in zip(pop, fitnesses):
  ind.fitness.values = fit
while generations<=500:
  CXPB, MUTPB = 0.75, 0.24
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
      nnew_row={'String':(ctiterations+1),'Generation':generations,'Score':ind.fitness.values[0],'can_apoptose_prob':ind[0],'can_apoptose_prob_inc':ind[1],'can_apoptose_time_window':ind[2],'can_apoptose_time_window_to_delete':ind[3],'can_divide_prob':ind[4],'can_divide_prob_inc':ind[5],'can_divide_time_window':ind[6],'can_divide_max':ind[7],'can_grow_prob':ind[8],'diameter_rate':ind[9]}

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

