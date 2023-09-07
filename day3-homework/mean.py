#!/usr/bin/env python

a=[1,3,4,6,7,43,7,5,89,32]

def meanlist(input_list):
	sum=0
	for i in input_list:
		sum+=i
	return sum/len(input_list)

print(meanlist(a))