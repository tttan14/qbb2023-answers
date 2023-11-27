#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sps
import statsmodels.formula.api as smf
import sys
from fasta import readFASTA

################################################
####################################################
# fname = sys.argv[1]
# fs = open(fname)

# fname1=sys.argv[2]
# score_matrix=pd.read_csv(fname1, index_col=0, delim_whitespace = True)


# gap_penalty=int(sys.argv[3])

# input_sequences = readFASTA(fs)

# seq1_id, sequence1 = input_sequences[0][0],input_sequences[0][1]
# seq2_id, sequence2 = input_sequences[1][0],input_sequences[1][1]

# F_matrix=np.zeros((len(sequence1)+1,len(sequence2)+1))
# TB_matrix=np.zeros((len(sequence1)+1,len(sequence2)+1),str)

# def initialize_matrices(m, n):
#     F_matrix = [[0] * (n + 1) for _ in range(m + 1)]
#     TB_matrix = [[''] * (n + 1) for _ in range(m + 1)]  # Renamed the traceback matrix to TB
#     return F_matrix, TB_matrix

# for i in range(1, F_matrix.shape[0]):
# 	for j in range(1,F_matrix.shape[1]):
# 		if sequence1[i-1]==sequence2[j-1]:
# 			d=F_matrix[i-1,j-1]+score_matrix.loc[str(sequence1[i-1]),str(sequence2[j-1])]
# 		else:
# 			d=F_matrix[i-1,j-1]+score_matrix.loc[str(sequence1[i-1]),str(sequence2[j-1])]
# 		h=F_matrix[i,j-1]+gap_penalty
# 		v=F_matrix[i-1,j]+gap_penalty

# 		F_matrix[i,j]=max(d,h,v)
# 		if max(d,h,v)==d:
# 			TB_matrix[i,j]='d'
# 		elif max(d,h,v)==h:
# 			TB_matrix[i,j]='h'
# 		else:
# 			TB_matrix[i,j]='v'

# seq1=str()
# seq2=str()

# gap_count=0
# # i=len(sequence1)-1
# # j=len(sequence2)-1
# alignment_score=0

# def populate_matrices(sequence1, sequence2, score_matrix, gap_penalty, F_matrix, TB_matrix):
#     m, n = len(sequence1), len(sequence2)

#     for i in range(1, m + 1):
#         for j in range(1, n + 1):
#             match = F_matrix[i - 1][j - 1] + score_matrix.loc[sequence1[i - 1], sequence2[j - 1]]
#             gap_seq1 = F_matrix[i - 1][j] + gap_penalty
#             gap_seq2 = F_matrix[i][j - 1] + gap_penalty

#             max_score = max(match, gap_seq1, gap_seq2)

#             F_matrix[i][j] = max_score

#             if max_score == match:
#                 TB_matrix[i][j] = 'd'  # Diagonal (match)
#             elif max_score == gap_seq1:
#                 TB_matrix[i][j] = 'u'  # Up (gap in sequence 1)
#             else:
#                 TB_matrix[i][j] = 'l'  # Left (gap in sequence 2)

#     return F_matrix, TB_matrix

# # while len(seq1)!=len(sequence1):
# # 		if TB_matrix[i,j]=="d":
# # 			i-=1
# # 			j-=1
# # 			seq1=seq1+sequence1[i]
# # 			seq2=seq2+sequence1[i]
# # 			alignment_score+=F_matrix[i,j]

# # 		elif TB_matrix[i,j]=='h':
# # 			gap_count+=1
# # 			j-=1
# # 			seq1=seq1+'-'
# # 			seq2=seq2+sequence2[j]
# # 			alignment_score+=F_matrix[i,j]

# # 		elif TB_matrix[i,j]=='v':
# # 			gap_count+=1
# # 			i-=1
# # 			seq1+=sequence1[i]
# # 			seq2+='-'
# # 			alignment_score+=F_matrix[i,j]
# def traceback_alignment(TB_matrix, sequence1, sequence2):
#     m, n = len(TB_matrix) - 1, len(TB_matrix[0]) - 1
#     alignment_seq1, alignment_seq2 = "", ""

#     while m > 0 or n > 0:
#         if m > 0 and n > 0 and TB_matrix[m][n] == 'd':
#             alignment_seq1 = sequence1[m - 1] + alignment_seq1
#             alignment_seq2 = sequence2[n - 1] + alignment_seq2
#             m -= 1
#             n -= 1
#         elif m > 0 and TB_matrix[m][n] == 'u':
#             alignment_seq1 = sequence1[m - 1] + alignment_seq1
#             alignment_seq2 = '-' + alignment_seq2
#             m -= 1
#         elif n > 0 and TB_matrix[m][n] == 'l':
#             alignment_seq1 = '-' + alignment_seq1
#             alignment_seq2 = sequence2[n - 1] + alignment_seq2
#             n -= 1

#     return alignment_seq1, alignment_seq2

# F_matrix, TB_matrix = initialize_matrices(len(sequence1), len(sequence2))
# F_matrix, TB_matrix = populate_matrices(sequence1, sequence2, score_matrix, gap_penalty, F_matrix, TB_matrix)

# aligned_seq1, aligned_seq2 = traceback_alignment(TB_matrix, sequence1, sequence2)

#     # Count gaps and calculate alignment score
# gaps_seq1 = aligned_seq1.count('-')
# gaps_seq2 = aligned_seq2.count('-')
# alignment_score = F_matrix[len(sequence1)][len(sequence2)]

# print('sequence1: ',aligned_seq1)
# print('sequence2: ',aligned_seq2)
# print('gap_count: ',gap_count)
# print('alignment_score: ',alignment_score)

# fs.close()

#!/usr/bin/env python


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

for i in range(len(sequence1)+1):
	TB_matrix[i,0]=i*''

for j in range(len(sequence2)+1):
	TB_matrix[0,j]=j*''

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
j=len(sequence2)
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





