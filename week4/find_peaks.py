#!/usr/bin/env python

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from scipy.stats import poisson

def load_bedgraph(fname, target, chromstart, chromend):
    # Create array to hold tag counts
    coverage = np.zeros(chromend - chromstart, int)

    #Read the file in line by line
    for line in open(fname):
        # Break the line into individual fields
        chrom, start, end, score = line.rstrip().split('\t')
        # Check if the data fall in our target region
        if chrom != target:
            continue
        start = int(start)
        end = int(end)
        if start < chromstart or end >= chromend:
            continue
        # Add tags to our array
        coverage[start-chromstart:end-chromend] = int(score)
    return coverage

def bin_array(data, binsize):
    # Create array to hold scores
    binned = np.zeros(data.shape[0], data.dtype)

    # For each position in the window, add to the score array
    for i in range(binsize):
        binned[i:data.shape[0] - binsize + i] += data[binsize//2:-binsize//2]
    return binned

def main():
    # Load file names and fragment width
    forward_sname, reverse_sname = sys.argv[1:3]
    frag_width = sys.argv[3]
    forward_cname, reverse_cname = sys.argv[4:6]

    # Define what genomic region we want to analyze
    chrom = "chr2R"
    chromstart = 10000000
    chromend =  12000000
    chromlen = chromend - chromstart

    # Load the sample bedgraph data, reusing the function we already wrote
    sample_forward = load_bedgraph(forward_sname, chrom, chromstart, chromend)
    sample_reverse = load_bedgraph(reverse_sname, chrom, chromstart, chromend)
 
    # Combine tag densities, shifting by our previously found fragment width 198
    combined_samp = sample_forward[0] + sample_reverse[:-198]

#     # Load the control bedgraph data, reusing the function we already wrote
    control_forward = load_bedgraph(forward_cname, chrom, chromstart, chromend)
    control_reverse = load_bedgraph(reverse_cname, chrom, chromstart, chromend)

#     # Combine tag densities
    combined_ctrl = control_forward[0]+ control_reverse[:-198] 

#     # Adjust the control to have the same coverage as our sample
    adj_index = float(np.sum(combined_samp)/np.sum(combined_ctrl))
    combined_ctrl = combined_ctrl * adj_index

#     # Create a background mean using our previous binning function and a 1K window
#     # Make sure to adjust to be the mean expected per base
    binned_ctrl = bin_array(combined_ctrl, 1000)/1000


#     # Find the mean tags/bp and make each background position the higher of the
#     # the binned score and global background score
    global_score= np.mean(combined_ctrl)
    found_tb=[]
    for i in range(len(binned_ctrl)):
        found_tb.append(max(binned_ctrl[i], global_score))
    found_tb=np.array(found_tb)

#     # Score the sample using a binsize that is twice our fragment size
#     # We can reuse the binning function we already wrote
    binned_samp = bin_array(combined_samp, 2*198)

    # Find the p-value for each position (you can pass a whole array of values
    # and and array of means). Use scipy.stats.poisson for the distribution.
    # Remeber that we're looking for the probability of seeing a value this large
    # or larger
    # Also, don't forget that your background is per base, while your sample is
    # per 2 * width bases. You'll need to adjust your background
    found_tb_adj=(found_tb*2)/1000
    p_values =  1 - poisson.cdf(binned_samp, found_tb_adj)

#     # Transform the p-values into -log10
#     # You will also need to set a minimum pvalue so you doen't get a divide by
#     # zero error. I suggest using 1e-250
    mini= 1e-250
    log_pval = -np.log10(np.maximum(p_values, mini))

    # Write p-values to a wiggle file
    # The file should start with the line
    # "fixedStep chrom=CHROM start=CHROMSTART step=1 span=1" where CHROM and
    # CHROMSTART are filled in from your target genomic region. Then you have
    # one value per line (in this case, representing a value for each basepair).
    # Note that wiggle files start coordinates at 1, not zero, so add 1 to your
    # chromstart. Also, the file should end in the suffix ".wig"
    #write_wiggle(log_pval,chrom,chromstart,str(forward_sname)[0:7]+'.wig')

    # Write bed file with non-overlapping peaks defined by high-scoring regions 
    write_bed(log_pval,chrom,chromstart,chromend,198,str(forward_sname)[0:7]+'.bed')

def write_wiggle(pvalues, chrom, chromstart, fname):
    output = open(fname, 'w')
    print(f"fixedStep chrom={chrom} start={chromstart + 1} step=1 span=1",
          file=output)
    for i in pvalues:
        print(i, file=output)
    output.close()

def write_bed(scores, chrom, chromstart, chromend, width, fname):
    chromlen = chromend - chromstart
    output = open(fname, 'w')
    while np.amax(scores) >= 10:
        pos = np.argmax(scores)
        start = pos
        while start > 0 and scores[start - 1] >= 10:
            start -= 1
        end = pos
        while end < chromlen - 1 and scores[end + 1] >= 10:
            end += 1
        end = min(chromlen, end + width - 1)
        print(f"{chrom}\t{start + chromstart}\t{end + chromstart}", file=output)
        scores[start:end] = 0
    output.close()


if __name__ == "__main__":
    main()