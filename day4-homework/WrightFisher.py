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
#Ex1 pt1

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

al_list=wright_fisher(0.5,1000)
#print("Number of generation to fixation:",len(al_list))

#Ex1 pt2

fig,ax=plt.subplots()

n=1
x_traject=[]
for j in range(len(al_list)):
	n+=1
	x_traject.append(n)

ax.plot(x_traject,al_list)

ax.set_xlabel("Generations")
ax.set_ylabel("Allele frequency")
fig.savefig("WrightFisher-Ex1.png")

plt.show()

######################################################################

#Ex2 pt1
# for i in range(30):
# 	x_traject=[]
# 	x_traject=wright_fisher(0.5,100)
# 	y_traject=[]
# 	n=1
# 	for j in range(len(x_traject)):
# 		n+=1
# 		y_traject.append(n)

# 	ax.plot(y_traject,x_traject)
#ax.set_xlabel("Generations")
#ax.set_ylabel("Allele frequency")
#fig.savefig("WrightFisher-Ex2.png")

#plt.show()

##############################################################

#Ex2 pt2 histogram
# x_traject=[]
# for i in range(1000):
#  	y_traject=wright_fisher(0.5,100)
#  	for j in y_traject:
#  		x_traject.append(len(y_traject))	
# plt.hist(x_traject)


# ax.set_xlabel("Times to fixation")
# ax.set_ylabel("Number of occurrences")
# fig.savefig("WrightFisher-Ex2pt2.png")

# plt.show()


##############################################################
#Ex3 pt1

# pop_input=[60,70,80,90,100]
# num_fix=[]
# aver_fix=0
# aver_plot=[]
# for i in range(5):
# 	al_freq=[]
# 	num_fix=[]
# 	sum=0
# 	for j in range(60):
# 		al_freq.append(wright_fisher(0.5,pop_input[i]))
# 		num_fix.append(len(al_freq[0])) #num fix one time
# 		al_freq=[]

# 	num_fix=np.array(num_fix) #num to fix for 60 times
# 	aver_fix=np.mean(num_fix) #average of the 60 times for one pop size
# 	aver_plot.append(aver_fix) #list of 5 pop size with average num to fix 


# ax.scatter(pop_input,aver_plot)
# ax.set_xlabel("Population sizes")
# ax.set_ylabel("Average number to fixations (60 trials)")
# fig.savefig("WrightFisher-Ex3pt1.png")
# plt.show()

###########################################################
#Ex3 pt2

# salf_input=[0.1,0.2,0.3,0.4,0.5]
# num_fix=[]
# aver_fix=0
# aver_plot=[]
# for i in range(5):
# 	al_freq=[]
# 	num_fix=[]
# 	sum=0
# 	for j in range(10):
# 		al_freq.append(wright_fisher(salf_input[i],1000))
# 		num_fix.append(len(al_freq[0])) #num fix one time
# 		al_freq=[]

# 	num_fix=np.array(num_fix) #num to fix for 60 times
# 	aver_fix=np.mean(num_fix) #average of the 60 times for one pop size
# 	aver_plot.append(aver_fix) #list of 5 pop size with average num to fix 

# ax.scatter(salf_input,aver_plot)
# ax.set_xlabel("Allele frequencies")
# ax.set_ylabel("Average number to fixations (10 trials)")
# fig.savefig("WrightFisher-Ex3pt2.png")
# plt.show()




