#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#EXERCISE 1
#########################################################################
#1.1

dnms=pd.read_csv('aau1043_dnm.csv',sep=',',header=0) 
proband=dnms.loc[:,'Proband_id'] 

##########################################################################
#1.2

a=0
b=0
deNovoCount={}
for i in proband:
	deNovoCount[i]=[a,b]

for i in range(len(proband)):
	if dnms.loc[i,'Phase_combined']=='mother':
			deNovoCount[dnms.loc[i,'Proband_id']][0]+=1
	elif dnms.loc[i,'Phase_combined']=='father':
			deNovoCount[dnms.loc[i,'Proband_id']][1]+=1

##########################################################################
#1.3

deNovoCountDF = pd.DataFrame.from_dict(deNovoCount, orient = 'index', columns = ['maternal_dnm', 'paternal_dnm'])

##########################################################################
#1.4

dwage=pd.read_csv('aau1043_parental_age.csv',sep=',',header=0,index_col=0)

##########################################################################
#1.5 

merged=pd.concat([deNovoCountDF,dwage],axis=1,join='inner')




