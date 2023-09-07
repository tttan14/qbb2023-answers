#!/usr/bin/env python

import sys

fname=sys.argv[1]
data=open(fname)
lines=data.readlines()

patientID="10"

def mean_pa(pa_id): #function that calculates the mean for specific patient
	fla_data=pa_dict[str(pa_id)] #dictionary uses []
	sum_pa=0
	for i in fla_data:
		sum_pa+=i
	return sum_pa/len(fla_data)

n=1
flare_data=[]
pa_dict={}
for i in lines:
	flare_data=[]
	i=i.rstrip('\n')
	clean=i.split(',')
	for j in range(len(clean)):
		flare_data.append(int(clean[j])) #making record integers
	key=str(n)
	score=flare_data
	pa_dict[key]=score #match patientID to flare-ups record
	n+=1


print(mean_pa(patientID))