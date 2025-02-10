import os
import pandas as pd
import xlrd
import re

excelfilepath="/home/john/Downloads/HeLa_cell_cycle_analysis_-_01Dec2023.xlsx"
excelfilepath2="/home/john/Downloads/HeLa_cells_viability_after_DAC_-_01Dec2023.xlsx"
xls = pd.ExcelFile(excelfilepath)
scale=1000
df2=pd.read_excel(excelfilepath2,header=None,sheet_name="normalized")
df2=df2.dropna(how="all")
df2=df2.fillna(method="ffill", axis=1)
df2=df2.reset_index(drop=True)
i=0
while(not(df2[0].iloc[i]=="values in %")):
      i=i+1
print(i)
df2=df2[:].iloc[:i]
df2=df2.tail(-1)
i=1
j=1
df2=df2.reset_index(drop=True)
print(df2)

print(df2[0].iloc[2])
ustr=0
utstr=0
uttstrt=0
for sheetname in xls.sheet_names:

     df1=pd.read_excel(excelfilepath,header=None, sheet_name=sheetname)
     df1=df1.dropna(how="all")
     df1=df1.fillna(method="ffill", axis=1)
     df1=df1.tail(-1)
     df1=df1.reset_index(drop=True)
     print(df1)
     treatedstr="current_time, N_vessels, N_cells, N_cells_pheno_0, N_cells_pheno_1, N_cells_pheno_1_Ap, N_cells_pheno_1_G1, N_cells_pheno_1_Sy, N_cells_pheno_1_G2, N_cells_pheno_1_Di, N_cells_pheno_1_Tr, N_cell_protrusions"
     untreatedstr="current_time, N_vessels, N_cells, N_cells_pheno_0, N_cells_pheno_1, N_cells_pheno_1_Ap, N_cells_pheno_1_G1, N_cells_pheno_1_Sy, N_cells_pheno_1_G2, N_cells_pheno_1_Di, N_cells_pheno_1_Tr, N_cell_protrusions"
     udf1=df1.copy(deep=True)
     udf2=df2.copy(deep=True)
     ign=0
     for col in udf1:
          if(ign==0): 
               ign=1
          else:
               if (udf1[col].iloc[0]!="Untreated"):
                    udf1=udf1.drop(col, axis=1)
     for col in udf2:
        if(ign==0): 
               ign=1
        else:
               if (udf2[col].iloc[0]!="Untreated"):
                    udf2=udf2.drop(col, axis=1)
    
     udf2=udf2.T.reset_index(drop=True).T
     print(udf1)
     print(udf2)

     maxm=len(udf2.axes[0])
     maxl=len(udf2.axes[1])
     maxi=len(udf1.axes[0])
     maxj=len(udf1.axes[1])
     print ("\n\n"+str(maxl)+" "+str(maxm)+" "+str(maxi)+" "+str(maxj))
     i=0
     j=1
     l=0
     m=1
     while(True):
             print(udf2[0].iloc[m])
             ctime=udf1[0].iloc[i]
             print(ctime)
             cctime=re.findall(r'\d+', ctime)[0]
             print(cctime)
             totalcells=float(udf2[l].iloc[m])*scale
             g0cells=round((float(udf1[j].iloc[i+1])/100)*totalcells)
             scells=round((float(udf1[j].iloc[i+2])/100)*totalcells)
             g2cells=round((float(udf1[j].iloc[i+3])/100)*totalcells)
             apcells=totalcells-g0cells-scells-g2cells
             untreatedstr=untreatedstr+"\n"+cctime+","+"0"+","+str(int(totalcells))+","+"0,"+str(int(totalcells))+","+str(int(apcells))+","+str(int(g0cells))+","+str(int(scells))+","+str(int(g2cells))+",0,0,0"
             m=m+1
             i=i+4
             #print(m)
            # print(i)
             
             if(i>=maxi or m>=maxm):
                ustr=ustr+1
                #print(untreatedstr)
                f1=open("untreated_"+str(ustr)+".csv","w")
                f1.write(untreatedstr)
                f1.close()
                untreatedstr="current_time, N_vessels, N_cells, N_cells_pheno_0, N_cells_pheno_1, N_cells_pheno_1_Ap, N_cells_pheno_1_G1, N_cells_pheno_1_Sy, N_cells_pheno_1_G2, N_cells_pheno_1_Di, N_cells_pheno_1_Tr, N_cell_protrusions"
                l=l+1
                m=1
                i=0
                if(l>=maxl):
                     l=0
                     j=j+1
                     if(j>=maxj): break

                     #if(l==2 and sheetname.__contains__("1")):
                  #l=l+2
     tdf1=df1.copy(deep=True)
     tdf2=df2.copy(deep=True)
     tdf22=df2.copy(deep=True)
     ign=0
     for col in udf1:
          if(ign==0): 
               ign=1
          else:
               if (tdf1[col].iloc[0]!="Treated"):
                    tdf1=tdf1.drop(col, axis=1)
     for col in tdf2:
        if(ign==0): 
               ign=1
        else:
               if (tdf2[col].iloc[0]!="Treated (0.3 μΜ)"):
                    tdf2=tdf2.drop(col, axis=1)
               if (tdf22[col].iloc[0]!="Treated (1 μΜ)"):
                    tdf22=tdf22.drop(col, axis=1)
    
     tdf2=tdf2.T.reset_index(drop=True).T
     tdf22=tdf22.T.reset_index(drop=True).T
     tdf1=tdf1.T.reset_index(drop=True).T
     print(tdf1)
     print(tdf2)
     print(tdf22)

     maxm=len(tdf2.axes[0])
     maxl=len(tdf2.axes[1])
     maxi=len(tdf1.axes[0])
     maxj=len(tdf1.axes[1])
     print ("\n\n"+str(maxl)+" "+str(maxm)+" "+str(maxi)+" "+str(maxj))
     i=0
     j=1
     l=0
     m=1
     while(True):
            if "1" in sheetname:
                print(tdf22[l].iloc[m])
                ctime=tdf1[0].iloc[i]
                print(ctime)
                cctime=re.findall(r'\d+', ctime)[0]
                print(cctime)
                totalcells=float(tdf22[l].iloc[m])*scale
                g0cells=round((float(tdf1[j].iloc[i+1])/100)*totalcells)
                scells=round((float(tdf1[j].iloc[i+2])/100)*totalcells)
                g2cells=round((float(tdf1[j].iloc[i+3])/100)*totalcells)
                apcells=totalcells-g0cells-scells-g2cells
                treatedstr=treatedstr+"\n"+cctime+","+"0"+","+str(int(totalcells))+","+"0,"+str(int(totalcells))+","+str(int(apcells))+","+str(int(g0cells))+","+str(int(scells))+","+str(int(g2cells))+",0,0,0"
                m=m+1
                i=i+4
                #print(m)
                #print(i)
                #print(j)
                if(i>=maxi or m>=maxm):
                    utstr=utstr+1
                    print("tstr: "+str(utstr))
                    f1=open("treated_"+sheetname+"_"+str(utstr)+".csv","w")
                    f1.write(treatedstr)
                    f1.close()
                    #print(treatedstr)
                    treatedstr="current_time, N_vessels, N_cells, N_cells_pheno_0, N_cells_pheno_1, N_cells_pheno_1_Ap, N_cells_pheno_1_G1, N_cells_pheno_1_Sy, N_cells_pheno_1_G2, N_cells_pheno_1_Di, N_cells_pheno_1_Tr, N_cell_protrusions"
                    l=l+1
                    m=1
                    i=0
                    if(l>=maxl):
                     l=0
                     j=j+1
                     if(j>=maxj):
                          break
            else:
                print(tdf2[l].iloc[m])
                ctime=tdf1[0].iloc[i]
                print(ctime)
                cctime=re.findall(r'\d+', ctime)[0]
                print(cctime)
                totalcells=float(tdf2[l].iloc[m])*scale
                g0cells=round((float(tdf1[j].iloc[i+1])/100)*totalcells)
                scells=round((float(tdf1[j].iloc[i+2])/100)*totalcells)
                g2cells=round((float(tdf1[j].iloc[i+3])/100)*totalcells)
                apcells=totalcells-g0cells-scells-g2cells
                treatedstr=treatedstr+"\n"+cctime+","+"0"+","+str(int(totalcells))+","+"0,"+str(int(totalcells))+","+str(int(apcells))+","+str(int(g0cells))+","+str(int(scells))+","+str(int(g2cells))+",0,0,0"
                m=m+1
                i=i+4
                #print(m)
                #print(i)
                #print(j)
                if(i>=maxi or m>=maxm):
                    uttstrt=uttstrt+1
                    print("tstr: "+str(uttstrt))
                    f1=open("treated_"+sheetname+"_"+str(uttstrt)+".csv","w")
                    f1.write(treatedstr)
                    f1.close()
                    #print(treatedstr)
                    treatedstr="current_time, N_vessels, N_cells, N_cells_pheno_0, N_cells_pheno_1, N_cells_pheno_1_Ap, N_cells_pheno_1_G1, N_cells_pheno_1_Sy, N_cells_pheno_1_G2, N_cells_pheno_1_Di, N_cells_pheno_1_Tr, N_cell_protrusions"
                    l=l+1
                    m=1
                    i=0
                    if(l>=maxl):
                     l=0
                     j=j+1
                     if(j>=maxj):
                          break