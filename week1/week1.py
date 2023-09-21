#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.formula.api as smf
import scipy.stats as sps


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

#EXERCISE 2
#########################################################################
#2.1

fig, ax = plt.subplots()
ax.scatter(merged.loc[:,'maternal_dnm'],merged.loc[:,'Mother_age'])
ax.set_xlabel("maternal de novo mutations")
ax.set_ylabel("maternal age")
ax.set_title( "maternal de novo mutations vs. maternal age" )
#fig.savefig( "ex2_a.png" )
#plt.show()

fig1, ax1 = plt.subplots()
ax1.scatter(merged.loc[:,'paternal_dnm'],merged.loc[:,'Father_age'])
ax1.set_xlabel("paternal de novo mutations")
ax1.set_ylabel("paternal age")
ax1.set_title( "paternal de novo mutations vs. paternal age" )
#fig1.savefig( "ex2_b.png" )
#plt.show()

#########################################################################
#2.2

merged_1=pd.concat([merged.loc[:,'maternal_dnm'],merged.loc[:,'Mother_age']],axis=1,join='inner')
model1 = smf.ols(formula='Mother_age ~1 + maternal_dnm',data=merged_1)
results1 = model1.fit()
print(results1.summary())
print('maternal pvalues: ',results1.pvalues['maternal_dnm'])

#########################################################################
#2.3

merged_2=pd.concat([merged.loc[:,'paternal_dnm'],merged.loc[:,'Father_age']],axis=1,join='inner')
model2 = smf.ols(formula='Father_age ~1 + paternal_dnm',data=merged_2)
results2 = model2.fit()
print(results2.summary())
print('paternal pvalues: ',results2.pvalues['paternal_dnm'])
#########################################################################
#2.5

fig2, ax2 = plt.subplots()
#merged.loc[:,'paternal_dnm'].plot.bar()
ax2.hist(merged.loc[:,'paternal_dnm'],label='paternal de novo mutations',alpha=0.5)
ax2.hist(merged.loc[:,'maternal_dnm'],label='maternal de novo mutations',alpha=0.5)
ax2.legend()
ax2.set_xlabel("de novo mutation counts")
ax2.set_ylabel("Occurrences")
#fig2.savefig( "ex2_c.png" )
plt.show()

#########################################################################
#2.6

mean1=merged_1.loc[:,'maternal_dnm'].mean()
mean2=merged_2.loc[:,'paternal_dnm'].mean()

std1=merged_1.loc[:,'maternal_dnm'].std()
std2=merged_2.loc[:,'paternal_dnm'].std()

print(sps.ttest_ind_from_stats(mean1, std1, len(merged_1.loc[:,'maternal_dnm']), mean2, std2, len(merged_2.loc[:,'paternal_dnm']), equal_var=True, alternative='two-sided'))








