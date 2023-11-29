#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# file=open('/Users/cmdb/qbb2023-answers/week5/annotated.vcf','r')
# print(file)
# n1=pd.DataFrame()
# for line in open('annotated.vcf'):
#     if line.startswith('#'):
#         continue
#     fields = line.rstrip('\n').split('\t')
#     fields.
# n1=pd.DataFrame(columns=[])
vcf_path = "/Users/cmdb/qbb2023-answers/week5/annotated.vcf"
vcf_df = pd.read_csv(vcf_path, sep='\t', comment='#', header=None)
# Set column names
vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '63','62','35','31','27','24','23','39','11','09']

DP_new=[]
samples=['63','62','35','31','27','24','23','39','11','09']
for sample in samples:
    for j in range(len(vcf_df['INFO'])):
        parts=vcf_df.loc[j,sample].split(':')
        if parts[2]!='.':
            DP_new.append(int(parts[2]))
fig, ax = plt.subplots(nrows=2, ncols=2)
ax[0,0].hist(DP_new,bins=100)
#ax.barh(bins[:-1], counts, height=np.diff(bins), color='blue', edgecolor='black')
#plt.title('Distribution of Read Depth Across Variants')
ax[0,0].set_xlabel("read depth")
ax[0,0].set_ylabel("counts")

GQ_new=[]
for sample in samples:
    for j in range(len(vcf_df['INFO'])):
        parts=vcf_df.loc[j,sample].split(':')
        if parts[1]!='.':
            GQ_new.append(float(parts[1]))

ax[0,1].hist(GQ_new,bins=100)
ax[0,1].set_xlabel("genotype quality")
ax[0,1].set_ylabel("counts")

AF_new=[]

for j in range(len(vcf_df['INFO'])):
    parts=vcf_df.loc[j,'INFO'].split(';')
    AF=parts[3]
    if ',' not in AF:
        AF_new.append(float(AF[3:]))
ax[1,0].hist(AF_new)
ax[1,0].set_xlabel("allele frequency")
ax[1,0].set_ylabel("counts")

types=[]
for j in range(len(vcf_df['INFO'])):
    parts=vcf_df.loc[j,'INFO'].split(';')
    nparts=parts[-1][4:10]
    if '|' == nparts[1]:
        types.append(nparts[0])
    elif '|' == nparts[2]:
        types.append(nparts[0:1])
    elif '|' == nparts[3]:
        types.append(nparts[0:2])
    elif '|' == nparts[4]:
        types.append(nparts[0:3])
# fig,ax=plt.subplots()
ax[1,1].hist(types)
ax[1,1].set_xlabel("predicted effects")
#ax[1,1].set_xticklabels(ax[1,1].get_xticks(), rotation=45)
ax[1,1].set_ylabel("counts")
plt.show()
