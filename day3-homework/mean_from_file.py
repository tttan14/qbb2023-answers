#!/usr/bin/env python

import sys

def meanlist(input_list):
	sum=0
	for i in input_list:
		sum+=i
	return sum/len(input_list)

fname=sys.argv[1]
data=open(fname)
lines=data.readlines()

intlist=[]
for i in lines:
	i=i.rstrip('\n')
	clean=i.split()
	intlist.append(int(clean[0]))

print(meanlist(intlist))
