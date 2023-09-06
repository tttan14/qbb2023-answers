#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

transcripts = np.loadtxt( "all_annotated.csv", delimiter=",", usecols=0, dtype="<U30", skiprows=1 )
print( "transcripts: ", transcripts[0:5] )

samples = np.loadtxt( "all_annotated.csv", delimiter=",", max_rows=1, dtype="<U30" )[2:]
print( "samples: ", samples[0:5] )

data = np.loadtxt( "all_annotated.csv", delimiter=",", dtype=np.float32, skiprows=1, usecols=range(2, len(samples) + 2) )
print( "data: ", data[0:5, 0:5] )

# Find row with transcript of interest
for i in range(len(transcripts)):
    if transcripts[i] == 'FBtr0073461': #Find data related to sisA
        row = i

# Find columns with samples of interest
cols = []
for i in range(len(samples)):
    if "female" in samples[i]:
        cols.append(i)
#male
colm=[]
for i in range(len(samples)):
	if "male" in samples[i] and "fe" not in samples[i]:
		colm.append(i)



# Subset data of interest
expression = data[row, cols]
expressionm=data[row,colm]
expressionm=np.array(expressionm)
m2expression=2*expressionm


# Prepare data
x = samples[cols]
y = expression
y1=expressionm
y2=m2expression

x=["10",'11','12','13','14A','14B','14C','14D']

# Plot data
fig, ax = plt.subplots()

ax.plot( x, y, label="female")
ax.plot(x,y1,label="male")
ax.plot(x,y2,label="2*male")

ax.set_title( "sisA(FBtr0073461)" )
ax.set_xlabel("developmental stage")
ax.set_ylabel("mRNA abundance (RPKM)")
ax.legend()


plt.tight_layout()
fig.savefig( "sisA-f+m+2m(annotated).png" )
plt.show()






