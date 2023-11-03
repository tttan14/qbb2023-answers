#!/usr/bin/env python

import sys

#from model_peaks import load_bedgraph, bin_array
import numpy as np
#import scipy.stats
import matplotlib.pyplot as plt
import pandas as pd
import math


# pcalist=pd.read_csv('plink.eigenvec',sep=' ',header=None)
# list1=pcalist.loc[:,2]
# list2=pcalist.loc[:,3]

# fig, ax = plt.subplots()
# ax.scatter(list1,list2)
# ax.set_xlabel("similarity among individuals")
# ax.set_ylabel("similarity among individuals")
# plt.show()

# frelist=pd.read_csv('plink.frq',delim_whitespace=True)
# fig, ax = plt.subplots()
# ax.hist(frelist.loc[:,'MAF'])
# ax.set_xlabel("allele frequency")
# ax.set_ylabel("counts")
# plt.show()
pvalGS=pd.read_csv('GS_gwas_results.assoc.linear',delim_whitespace=True)

pval1=pd.read_csv('phenotype_gwas_results.assoc.linear',delim_whitespace=True)
# fig, (ax1,ax2) = plt.subplots(2)
# a=np.array(pval1.loc[:,'P'])
# a=np.log10(a)
# b=np.array(pvalGS.loc[:,'P'])
# b=np.log10(b)
# ax1.scatter(range(len(pval1.loc[:,'P'])),a*-1,c=np.where(pval1.loc[:,'P']>(10**(-5)), 'g', 'b')) 
# ax2.scatter(range(len(pvalGS.loc[:,'P'])),b*-1,c=np.where(pval1.loc[:,'P']>(10**(-5)), 'g', 'b'))
# ax1.set_xlabel("SNP(CB1908_IC50)")
# ax1.set_ylabel("-log10(P-value)")
# ax2.set_xlabel("SNP(GS451_IC50)")
# ax2.set_ylabel("-log10(P-value)")
# plt.show()

geno=pd.read_csv('genotypes.vcf',delim_whitespace=True, skiprows=27)
pmin=pvalGS.loc[:,'P'].min()
row_num = pvalGS[pvalGS.loc[:,'P'] == pmin].index 
row1=pvalGS.iloc[2650065]
genotypes=geno[geno.loc[:,'ID']=='rs7257475']

val=pd.read_csv('GS451_IC50.txt',delim_whitespace=True)

homo=np.array()
hete=np.array()
for i in len(val.loc[:,'GS451_IC50']):
	if genotypes.loc[i,'FID'] 
