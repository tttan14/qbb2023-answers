#!/usr/bin/env python

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
from matplotlib import pyplot as plt
import seaborn as sns


# read in data
counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

# read in metadata
metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)

# normalize
counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]

# log
counts_df_logged = np.log2(counts_df_normed + 1)

# merge with metadata
full_design_df = pd.concat([counts_df_logged, metadata], axis=1)

#1.1
# list1 = pd.DataFrame(full_design_df.loc['GTEX-113JC',:])

# trueRows = [ i  != 0.0 for i in list1['GTEX-113JC'] ]
# list1 = list1.loc[trueRows, :]
# trueRows = [ i  != '50-59' for i in list1['GTEX-113JC'] ]
# list1 = list1.loc[trueRows, :]
# list1_log = np.array(list1['GTEX-113JC'].tolist())
# list1_log = np.log2(list1_log)
# list1['logged'] = list1_log

# fig,ax = plt.subplots()
# ax.hist(list1['logged'])
# ax.set_ylabel("occurrences")
# ax.set_xlabel("logged normalized counts")

# plt.show()

# list2 = pd.DataFrame(full_design_df.loc[:,'MXD4'])
# list2['SEX'] = full_design_df.loc[:,'SEX']
# trueRows = [ i  == 1 for i in list2['SEX'] ]
# male = list2.loc[trueRows, :]
# trueRows = [ i  == 2 for i in list2['SEX'] ]
# female = list2.loc[trueRows, :]


# fig,ax = plt.subplots()
# ax.hist(male['MXD4'], bins = 70, alpha = 0.3, label = 'Male')
# ax.hist(female['MXD4'],bins = 70, alpha = 0.3, label = 'Female')
# ax.legend()
# ax.set_ylabel("occurrences")
# ax.set_xlabel("logged normalized counts")
# plt.show()

#3
# age = full_design_df.loc[:,'AGE'].tolist()
# age = set(age)
# age = list(age)
# list3 = pd.DataFrame(index = age,columns = ['counts'])

# trueRows = [ i  == '50-59' for i in full_design_df.loc[:,'AGE'] ]
# g50 = full_design_df.loc[trueRows, :]
# trueRows = [ i  == '20-29'for i in full_design_df.loc[:,'AGE'] ]
# g20 = full_design_df.loc[trueRows, :]
# trueRows = [ i  == '60-69'for i in full_design_df.loc[:,'AGE'] ]
# g60 = full_design_df.loc[trueRows, :]
# trueRows = [ i  == '30-39'for i in full_design_df.loc[:,'AGE'] ]
# g30 = full_design_df.loc[trueRows, :]
# trueRows = [ i  == '70-79'for i in full_design_df.loc[:,'AGE'] ]
# g70 = full_design_df.loc[trueRows, :]
# trueRows = [ i  == '40-49'for i in full_design_df.loc[:,'AGE'] ]
# g40 = full_design_df.loc[trueRows, :]
# list3.loc['50-59','counts'] = len(g50)
# list3.loc['60-69','counts'] = len(g60)
# list3.loc['30-39','counts'] = len(g30)
# list3.loc['70-79','counts'] = len(g70)
# list3.loc['20-29','counts'] = len(g20)
# list3.loc['40-49','counts'] = len(g40)

# ax = list3.plot.bar(rot=0)

# plt.xlabel('age group')
# plt.ylabel('counts')

# # Show the plot
# plt.show()

#4
list4 = pd.DataFrame(full_design_df.loc[:, 'LPXN'])
idlist = list4.index.tolist()

for i in idlist:
	age = full_design_df.loc[i,'AGE']
	list4.loc[i,'AGE'] = age
	sex = full_design_df.loc[i,'SEX']
	if sex == 1:
		list4.loc[i,'SEX'] = 'Male'
	else:
		list4.loc[i,'SEX'] = 'Female'


plt.figure(figsize=(10, 6))
sns.violinplot(x='AGE', y='LPXN', hue='SEX', data=list4, split=True)

plt.show()


