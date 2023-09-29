#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sps
import statsmodels.formula.api as smf

#EXERCISE 1
#########################################################

# def simulate_coverage(coverage,genome_len,read_len,figname):
	
# 	coverage_array=np.zeros(genome_len)
# 	#generate an array of zeros as genome
# 	#0 coverage at the start

# 	num_reads=int(coverage*genome_len/read_len)
# 	#calculate number of reads needed for the coverage

# 	low=0 #start position of genome
# 	high=genome_len-read_len 
# 	#end position of genome given length of reads
# 	#last piece read map to the last few bp of at end of genome

# 	start_positions=np.random.randint(low=0,high=high+1,size=num_reads) 
# 	#high value is exclusive

# 	#loop through genome positions
# 	for start in start_positions:
# 		coverage_array[start:start+read_len]+=1

# 	#create a np array of range of coverage array
# 	x=np.arange(0,max(coverage_array)+1)
# 	#x represents all the possible coverage for each bp

# 	#count how many nonzeros to give number of zero coverage
# 	sim_0cov=genome_len-np.count_nonzero(coverage_array)
# 	sim_0cov_pct=(sim_0cov/genome_len)*100
# 	print(f'In the simulation, there are {sim_0cov} bases with 0 coverage')
# 	print(f'This is {sim_0cov_pct}% of the genome.')

# 	#Get poisson distribution
# 	y_poisson=sps.poisson.pmf(x,mu=coverage)*genome_len
# 	#to get similar sizes --> times genome size
# 	#generates an array for probability of having each possible coverage for each bps(x)
# 	#Poisson is giving what we should expect without simulation
# 	#It's good when distribution matches the Poisson expection

# 	#Get normal distribution
# 	#cannot determine probability for single event(value) --> represents likelihoods
# 	#for a range of events happening
# 	y_normal=sps.norm.pdf(x,loc=coverage,scale=np.sqrt(coverage))*genome_len

# 	fig,ax=plt.subplots()
# 	ax.hist(coverage_array,bins=x,label='Simulation')
# 	ax.plot(x,y_poisson,label='Poisson')
# 	ax.plot(x,y_normal,label='Normal')
# 	ax.set_xlabel('Coverage')
# 	ax.set_ylabel('Frequency (bp)')
# 	ax.legend()
# 	fig.tight_layout()
# 	fig.savefig(figname)
# 	plt.show()

# #when increase coverage --> increase bases covered 
# simulate_coverage(3,1_000_000,100,'ex1_3x_cov.png')
# simulate_coverage(10,1_000_000,100,'ex1_10x_cov.png')
# simulate_coverage(30,1_000_000,100,'ex1_30x_cov.png')

#EXERCISE 1
#########################################################
k=3

graph=set()

reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

for read in reads: 
	for j in range(len(read)-k):
		kmer1=read[j:j+k]
		kmer2=read[j+1: j+1+k]
		graph.add(kmer1+'->'+kmer2)

f=open('digraph.txt','w')
f.write('digraph {\n')
for i in graph:
	f.write(i+'\n')

f.write('}\n')
f.close()



