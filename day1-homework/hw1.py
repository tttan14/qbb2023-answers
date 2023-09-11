#!/usr/bin/env python

import numpy as np

#read csv file
f=open("/Users/cmdb/Desktop/swc-python/data/inflammation-01.csv",'r')
lines=f.readlines()
nlines=[]
fla=[]
fla_data=[]
for i in lines:
    nlines=i.rstrip() #get rid of spaces
    fla=nlines.split(',')
    fla_data.append(fla)
    #fla_data.append(int(fla_data[i]))
pa_5=fla_data[4]
print("Question1: ", end="")
print(pa_5[0],pa_5[9],pa_5[-1])
nfla=[]

for i in range(len(fla_data)): #turn data all into integers
    exfla=[]
    for j in range(len(fla_data[i])):
        exfla.append(int(fla_data[i][j]))
    nfla.append(exfla)          

avg=0
aver_pa=[]
nfla=np.array(nfla)
for i in nfla:
    avg=0
    avg=np.mean(i)
    aver_pa.append(avg) #calculate average flare-ups for each patient and put them in one array

print("Question2: "+str(aver_pa[0:10]))
print("Question3 max: "+str(np.max(aver_pa)))
print("Question3 min: "+str(np.min(aver_pa)))


pa_1=nfla[0]
pa_5=nfla[4]
diff_15=np.subtract(pa_1,pa_5) #calculate difference between pa1 and pa5
print("Question4: "+str(diff_15))

f.close()

