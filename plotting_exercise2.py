#!/usr/bin/env python

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
from matplotlib import pyplot as plt
import seaborn as sns

# read in data
meta_df = pd.read_csv("age_gaps.csv", index_col = 0)

# Distribution of age difference
# age_diff = meta_df['age_difference']
# fig1, ax1 = plt.subplots()
# ax1.hist(age_diff)
# ax1.set_xlabel('age_difference')
# ax1.set_ylabel('Frequency')
# ax1.set_title('Distribution of age difference')
# plt.tight_layout()

# # Distribution of couple numbers
# No_couple = meta_df['couple_number']
# fig2, ax2 = plt.subplots()
# ax2.hist(No_couple)
# ax2.set_xlabel('couple_number')
# ax2.set_ylabel('Frequency')
# ax2.set_title('Distribution of couple_number')
# plt.tight_layout()

# # Gendar
# diff_gendar = meta_df[meta_df['character_1_gender'] != meta_df['character_2_gender']]
# diff_gendar['diff_gender'] = 'heterosexual'
# same_gendar = meta_df[meta_df['character_1_gender'] == meta_df['character_2_gender']]
# same_gendar['diff_gender'] = 'homosexual'

# combined_df = pd.concat([diff_gendar, same_gendar])
# plt.figure(figsize=(8, 6))
# sns.violinplot(data=combined_df, x='diff_gender', y='age_difference')
# plt.title('Age gap distribution in homosexual vs. heterosexual couples ')
# plt.xlabel('Couple gender')
# plt.ylabel('Age gaps')
# plt.show()

# # Relation between max_actor_age and age_difference
# idlist = meta_df.index.tolist()
# for i in idlist:

# 	age1 = meta_df.loc[i,'actor_1_age'].tolist()
# 	age2 = meta_df.loc[i,'actor_2_age'].tolist()
# 	meta_df.loc[i,'max_actor_age']  = max(age1,age2)

# fig1, ax1 = plt.subplots()
# ax1.scatter(meta_df['age_difference'],meta_df['max_actor_age'])
# ax1.set_xlabel('age difference')
# ax1.set_ylabel('larger actor age within couple')
# ax1.set_title('Age difference vs. larger actor age within couple')
# plt.tight_layout()

# plt.show()

# Relation between moview release year and age difference

plt.bar(meta_df['release_year'], meta_df['age_difference'])
plt.xlabel('movie release year')
plt.ylabel('age difference')

plt.show()

