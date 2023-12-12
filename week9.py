#!/usr/bin/env python

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats import multitest
from pydeseq2 import preprocessing
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt

# read in data
# counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

# # read in metadata
# metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)

# counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]
# counts_df_normed = np.log2(counts_df_normed + 1)

# full_design_df = pd.concat([counts_df_normed, metadata], axis=1)

# # model = smf.ols(formula = 'Q("DDX11L1") ~ SEX', data=full_design_df)
# # results = model.fit()

# # slope = results.params[1]
# # pval = results.pvalues[1]

# geneName = counts_df.columns.tolist()
# header_list=['geneName','slope','p-val','rejected','p_corrected']
# regGene = pd.DataFrame(columns=header_list)
# ind=0

# for i in geneName:
# 	model = smf.ols(formula = 'Q(i) ~ SEX', data=full_design_df)
# 	results = model.fit()

# 	slope = results.params[1]
# 	pval = results.pvalues[1]	

# 	regGene.loc[ind,'geneName']=i
# 	regGene.loc[ind,'slope']=slope
# 	regGene.loc[ind,'p-val']=pval

# 	ind+=1

# pval=np.array(regGene['p-val'].tolist())
# rejected,corrected_p = statsmodels.stats.multitest.fdrcorrection(pval, alpha=0.10, method='fdr_bh', is_sorted=False)

# for i in range(len(rejected)):

# 	regGene.loc[i,'rejected']=rejected[i]
# 	regGene.loc[i,'p_corrected']=corrected_p[i]

# #regGene.to_csv("regGene.csv")
# #regGene=pd.read_csv('regGene.csv')


# trueRows = [ i  == True for i in regGene['rejected'] ]

# fdr_list = regGene.loc[trueRows, :]	
# fdr_list.to_csv('fdr_list1.csv')
fdr_list=pd.read_csv('fdr_list.csv')
# print(max(fdr_list['p_corrected']))
# trueRows = [ i  != 'NaN'  for i in fdr_list['p_corrected'] ]
# fdr_list = fdr_list.loc[trueRows, :]	
# print(fdr_list)

#2
# dds = DeseqDataSet(counts=counts_df, metadata=metadata, design_factors="SEX")
# dds.deseq2()
# stat_res = DeseqStats(dds)
# stat_res.summary()
# results = stat_res.results_df

# results.to_csv("results.csv")
results=pd.read_csv('results.csv')
results = results.rename(columns={'Unnamed: 0': 'geneName'})
results1 = results.dropna()
common_gene = fdr_list[fdr_list['geneName'].isin(results1['geneName'])]['geneName']
# overlap_fdr = (len(common_gene)/len(fdr_list['geneName']))*100
# overlap_res = (len(common_gene)/len(results1['geneName']))*100
# print('Jaccard index step1: ',overlap_fdr)
# print('Jaccard index step2: ',overlap_res)



#plot
#sig_change = results[(df['padj'] < fdr_threshold) & (abs(df['log2FoldChange']) > 1)]
padj=np.array(results1['padj'].tolist())
neglog10p=np.log10(padj)*-1
fc=np.array(results1['log2FoldChange'].tolist())
common_list = common_gene.tolist()
trueR = [i in common_list for i in results1['geneName']]
common_gene = results1.loc[trueR,:]

trueRows = [ i  > 1  for i in common_gene['log2FoldChange'] ]
fc_filteredres = common_gene.loc[trueRows, :]

padj2=np.array(fc_filteredres['padj'].tolist())
neglog10p2=np.log10(padj2)*-1
fc2=np.array(fc_filteredres['log2FoldChange'].tolist())


fig,ax=plt.subplots()
ax.scatter(fc,neglog10p)
ax.scatter(fc2,neglog10p2,color='orange')

ax.set_xlabel("log2FoldChange")
ax.set_ylabel("-log10(p_value)")

plt.show()