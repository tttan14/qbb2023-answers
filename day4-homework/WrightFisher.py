#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#####################################################################
#Sudo code helps thinking


#Get a starting frequency and population size
#Input these as parameters for a function


#While allele frequency is between 0 and 1:
#	Make empty list to store allele frequency
#	Get # of successes for next generation from numpy binomial distribution
#	--> convert # of successes to new allele frequency
#	Store allele frequency into a list
#	Return the list


#Return a list of allele frequency at each time pt
#The # of generations to fixation is the length of the returned list


#######################################################################


def wright_fisher(start_AlFreq, pop_size):
	
	new_AlFreq_list=[start_AlFreq] #new list
	AlFreq=0
	
	while start_AlFreq>0 and start_AlFreq<1:
		numSuc=np.random.binomial(pop_size*2,start_AlFreq) #generate # of success
		AlFreq=numSuc/(pop_size*2) #calculate new allele frequency
		#print(AlFreq)
		new_AlFreq_list.append(AlFreq)
		#print(new_AlFreq)
		start_AlFreq=AlFreq

	return new_AlFreq_list

print(wright_fisher(0.5,100))
