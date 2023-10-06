#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sps
import statsmodels.formula.api as smf
import sys
from fasta import readFASTA

################################################
# seq1='TGTTACGG'
# seq2='GGTTGACTA'

# print(sys.argv[1])

# print('Arguments passed: ',sys.argv[1:])

# match_score=float(sys.argv[1])
# mismatch_score=float(sys.argv[2])
# gap_penalty=float(sys.argv[3])


# # lis_of_list=[[1,2,3],[4,5,6],[7,8,9]]
# # my_2d_array=np.array(lis_of_list)
# # #print(my_2d_array.shape[0])

# # my_2d_array[1,1]=314

# # for row in my_2d_array:
# # 	for val in row:
# # 		print(val)


# # for row_index in range(0,my_2d_array.shape[0]):
# # 	for col_index in range(0, my_2d_array[1]):
# # 		print(my_2d_array[row_index,col_index])
# # 		#looping through the array

# F_matrix=np.zeros((len(seq1)+1,len(seq2)+1))

# for i in range(len(seq1)+1):
# 	F_matrix[i,0]=i*gap_penalty

# for j in range(len(seq2)+1):
# 	F_matrix[0,j]=j*gap_penalty

# for i in range(1, F_matrix.shape[0]):
# 	for j in range(1,F_matrix.shape[1]):
# 		if seq1[i-1]==seq2[j-1]:
# 			d=F_matrix[i-1,j-1]+match_score
# 		else:
# 			d=F_matrix[i-1,j-1]+mismatch_score
# 		h=F_matrix[i,j-1]+gap_penalty
# 		v=F_matrix[i-1,j]+gap_penalty

# 		F_matrix[i,j]=max(d,h,v)

#print(F_matrix)

####################################################
fname = sys.argv[1]
fs = open(fname)

fname1=sys.argv[2]
score_matrix=pd.read_csv(fname1, index_col=0, delim_whitespace = True)


gap_penalty=int(sys.argv[3])

input_sequences = readFASTA(fs)

seq1_id, sequence1 = input_sequences[0][0],input_sequences[0][1]
seq2_id, sequence2 = input_sequences[1][0],input_sequences[1][1]

F_matrix=np.zeros((len(sequence1)+1,len(sequence2)+1))
TB_matrix=np.zeros((len(sequence1)+1,len(sequence2)+1),str)

for i in range(len(sequence1)+1):
	F_matrix[i,0]=i*gap_penalty

for j in range(len(sequence2)+1):
	F_matrix[0,j]=j*gap_penalty

for i in range(1, F_matrix.shape[0]):
	for j in range(1,F_matrix.shape[1]):
		if sequence1[i-1]==sequence2[j-1]:
			d=F_matrix[i-1,j-1]+score_matrix.loc[str(sequence1[i-1]),str(sequence2[j-1])]
		else:
			d=F_matrix[i-1,j-1]+score_matrix.loc[str(sequence1[i-1]),str(sequence2[j-1])]
		h=F_matrix[i,j-1]+gap_penalty
		v=F_matrix[i-1,j]+gap_penalty

		F_matrix[i,j]=max(d,h,v)
		if max(d,h,v)==d:
			TB_matrix[i,j]='d'
		elif max(d,h,v)==h:
			TB_matrix[i,j]='h'
		else:
			TB_matrix[i,j]='v'

seq1=str()
seq2=str()

gap_count=0
i=len(sequence1)-1
j=len(sequence2)-1
alignment_score=0

if TB_matrix[i,j]=="d":
	seq1=seq1+sequence1[i]
	seq2=seq2+sequence1[i]
	alignment_score+=F_matrix[i,j]
while len(seq1)!=len(sequence1):
		if TB_matrix[i,j]=="d":
			i-=1
			j-=1
			seq1=seq1+sequence1[i]
			seq2=seq2+sequence1[i]
			alignment_score+=F_matrix[i,j]

		elif TB_matrix[i,j]=='h':
			gap_count+=1
			j-=1
			seq1=seq1+'-'
			seq2=seq2+sequence2[j]
			alignment_score+=F_matrix[i,j]

		elif TB_matrix[i,j]=='v':
			gap_count+=1
			i-=1
			seq1+=sequence1[i]
			seq2+='-'
			alignment_score+=F_matrix[i,j]

print('sequence1: ',seq1)
print('sequence2: ',seq2)
print('gap_count: ',gap_count)
print('alignment_score: ',alignment_score)

fs.close()



