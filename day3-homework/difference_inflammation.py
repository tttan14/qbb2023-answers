#!/usr/bin/env python

import sys

fname=sys.argv[1]
data=open(fname)
lines=data.readlines()

patientID1="10"
patientID2="22"


def diff_pa(pa_id1,pa_id2):
	diff_data=[]
	pa1_data=pa_dict[str(pa_id1)]
	pa2_data=pa_dict[str(pa_id2)]
	for i in range(len(pa1_data)):
		diff_data.append(pa2_data[i]-pa1_data[i])
	return diff_data



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

print(diff_pa(patientID1,patientID2))