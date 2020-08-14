# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 16:46:11 2019

@author: André Then

Finds all combinations of 5 features which unambiguously classify the 
proteinogenic amino acids. All features passed in an excel-file are considered.
See example excel-file in the repository for the required format.
"""

import pandas as pd

#function which identifies the solutions by recursion. Considers all features
#in searchdf. maxgroupsize is the maximum tolerated amount of elements after
#the first two separations took place. So for separating with 5 features this
#parameter has to be initially set to 8. Classification with more features
#requires this parameter to be set to 8*(2^n) where n is the amount of features
#added (e.g. n=1 for six feature classification).
def findsets(searchdf, maxgroupsize, sets=[], IDs=[], adjustindex=0):
    if sets==[]:
        for index, row in searchdf.iterrows():
            for index2, row2 in searchdf[index+1:].iterrows():
                Lintersection=set(row["leftset"])&set(row2["leftset"])
                Ldifference=set(row["leftset"])-set(row2["leftset"])
                if len(list(Lintersection)) > maxgroupsize or len(list(Ldifference)) > maxgroupsize:
                    continue
                Rintersection=set(row["rightset"])&set(row2["rightset"])
                Rdifference=set(row["rightset"])-set(row2["rightset"])
                if len(list(Rintersection)) > maxgroupsize or len(list(Rdifference)) > maxgroupsize:
                    continue
                IDlist=[row["ID"],row2["ID"]]
                setlist=[Lintersection,Ldifference,Rintersection,Rdifference]
                findsets(searchdf, maxgroupsize/2, setlist, IDlist, index2 + 1)
    else:
       setlist=sets
       for index, row in searchdf[adjustindex:].iterrows():
           newsetlist=[]
           for group in setlist:
               Suitable=False
               Lintersection=group&set(row["leftset"])
               Rintersection=group-set(row["leftset"])
               if len(list(Lintersection)) > maxgroupsize or len(list(Rintersection)) > maxgroupsize:
                   break
               newsetlist.append(Lintersection)
               newsetlist.append(Rintersection)
               Suitable=True
           if Suitable==True:
               IDlist=IDs+[row["ID"]]
               if maxgroupsize/2==0.5 or [len(x) for x in newsetlist].count(1)==20:
                   allsolutions.append(IDlist)
               else:
                   findsets(searchdf, maxgroupsize/2, newsetlist, IDlist, index + 1)

#renumerates the index to make sure the dataframe is numerated correctly from 
#0 to length(df) - 1.                 
def renameIndex(dfToRename):
    i=0
    for index, row in dfToRename.iterrows():
        dfToRename.rename(index={index:i}, inplace=True)
        i+=1
    return dfToRename
           
#red input excel-file
df=pd.read_excel("D:\Macha\AAindexInterpretable4.xlsx")
#get a list of all the amino acid one letter codes.
AAcols=list(df.columns[2:])
df["median"]=df.median(axis=1)
allleftsets=[]
allrightsets=[]
#build left and right sets for each feature based on if the numerical value
#for the amino acids is smaller or larger than the median.
for index, row in df.iterrows():
   leftset=[]
   rightset=[]
   for AA in AAcols:
       if row[AA] <= row["median"]:
           leftset.append(AA)
       else:
           rightset.append(AA)
   allleftsets.append(leftset)
   allrightsets.append(rightset)
   
       
df["leftset"]=allleftsets
df["rightset"]=allrightsets
#Drop features which result in separations not able to participate in an optimal
#solution (5 features).
for index, row in df.iterrows():
    if (len(row["leftset"]) < 4) or (len(row["leftset"]) > 16):
        df=df.drop(index=index)

df = renameIndex(df)

allsolutions=list()
findsets(df, 8)
    
#['R', 'Q', 'H', 'I', 'L', 'K', 'M', 'F', 'W', 'Y']
#['R', 'D', 'C', 'Q', 'M', 'F', 'S', 'T', 'W', 'Y']
#['R', 'N', 'D', 'Q', 'E', 'H', 'K', 'P', 'S', 'Y']
#['A', 'D', 'C', 'Q', 'E', 'H', 'L', 'K', 'M', 'F']
#['N', 'C', 'E', 'H', 'I', 'L', 'F', 'W', 'Y', 'V']