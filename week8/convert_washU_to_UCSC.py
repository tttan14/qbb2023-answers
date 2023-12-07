#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
	baitmap_file = '/Users/cmdb/qbb2023-answers/week8/h19_chr20and21.baitmap'
	washU_txt_file = '/Users/cmdb/qbb2023-answers/week8/output_washU_text.txt'
	output_bed_file = '/Users/cmdb/qbb2023-answers/week8/out_put_bed.bed'
	header_line= 'track type=interact name="pCHIC" description="Chromatin interactions" useScore=on maxHeightPixels=200:100:50 visibility=full\n'
	with open(output_bed_file, 'w') as output_file:
		output_file.write(header_line)

	bait_map = read_baitmap(baitmap_file)
	washU=read_baitmap(washU_txt_file)
	max_strength = float(max(washU[2].tolist()))

	header_bed= ["chrom", "chromStart", "chromEnd", "name", "score", "value", "exp", "color", "sourceChrom", "sourceStart", "sourceEnd", "sourceName", "sourceStrand", "targetChrom", "targetStart", "targetEnd", "targetName", "targetStrand"]
	emptybed = pd.DataFrame(columns=header_bed)
	ind=0
	with open(washU_txt_file, 'r') as washU_file:
		for line in washU_file:

			fields = line.rstrip().split('\t')

			chrom=fields[0].rstrip().split(',')[0]
			fragment1=fields[0].rstrip().split(',')
			fragment2=fields[1].rstrip().split(',')

			start1=int(fragment1[1])
			end1=int(fragment1[2])
			start2=int(fragment2[1])
			end2=int(fragment2[2])

			emptybed.loc[ind,'chrom']=chrom
			emptybed.loc[ind,'chromStart']=min(start2,start1)
			emptybed.loc[ind,'chromEnd']=max(end2,end1)
			emptybed.loc[ind,'name']='.'
			emptybed.loc[ind,'exp']='.'
			emptybed.loc[ind,'color']='0'
			emptybed.loc[ind,'sourceStrand']='+'
			emptybed.loc[ind,'value']=float(fields[2])

			strength=float(fields[2])
			score = int(strength/max_strength*1000)
			emptybed.loc[ind,'score']=score



			#determine bait
			if int(start1) in bait_map[1].tolist() and int(end1) in bait_map[2].tolist():
				emptybed.loc[ind,'sourceStart']=start1
				emptybed.loc[ind,'sourceEnd']=end1
				bait_in = bait_map.index[bait_map[1] == int(start1)]
				bait_in = int(bait_in[0])
				emptybed.loc[ind,'sourceChrom']=bait_map.loc[bait_in,0]
				emptybed.loc[ind,'sourceName']=bait_map.loc[bait_in,3]

				if int(start2) in bait_map[1].tolist() and int(end2) in bait_map[2].tolist() :
					emptybed.loc[ind,'targetStart']=start2
					emptybed.loc[ind,'targetEnd']=end2
					bait_in = bait_map.index[bait_map[1] == int(start2)]
					bait_in = int(bait_in[0])
					emptybed.loc[ind,'targetChrom']=bait_map.loc[bait_in,0]
					emptybed.loc[ind,'targetName']=bait_map.loc[bait_in,3]
					emptybed.loc[ind,'targetStrand']='+'	
				
				else:
					emptybed.loc[ind,'targetStart']=start2
					emptybed.loc[ind,'targetEnd']=end2
					emptybed.loc[ind,'targetChrom']=chrom
					emptybed.loc[ind,'targetName']='.'	
					emptybed.loc[ind,'targetStrand']='-'					
			elif int(start2) in bait_map[1].tolist() and int(end2) in bait_map[2].tolist():
				emptybed.loc[ind,'targetStart']=start1
				emptybed.loc[ind,'targetEnd']=end1
				emptybed.loc[ind,'targetChrom']=chrom
				emptybed.loc[ind,'targetName']='.'	
				emptybed.loc[ind,'targetStrand']='-'
				emptybed.loc[ind,'sourceStart']=start2
				emptybed.loc[ind,'sourceEnd']=end2
				bait_in = bait_map.index[bait_map[1] == int(start2)]
				bait_in = int(bait_in[0])
				emptybed.loc[ind,'sourceChrom']=bait_map.loc[bait_in,0]
				emptybed.loc[ind,'sourceName']=bait_map.loc[bait_in,3]



			ind+=1

	emptybed.to_csv('/Users/cmdb/qbb2023-answers/week8/out_put_bed.bed',sep='\t',mode='a',header=False, index=False)

	#top pro-pro
	pro_list=emptybed[emptybed['targetStrand']=='+']
	sorted_pro = pro_list.sort_values(by='score', ascending=False)
	
	#top pro-enh
	enh_list=emptybed[emptybed['targetStrand']=='-']
	sorted_enh = enh_list.sort_values(by='score', ascending=False)	


def read_baitmap(baitmap_file):
	bait_map = pd.read_csv(baitmap_file, sep='\t',header=None)
    
	return bait_map


if __name__ == "__main__":
	main()