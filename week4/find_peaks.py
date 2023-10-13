#!/usr/bin/env python

import sys

from model_peaks import load_bedgraph, bin_array
import numpy
import scipy.stats
import matplotlib.pyplot as plt


def main():
    # Load bedgraph datasets
    forward_fname, reverse_fname, out_fname = sys.argv[1:4]

    # Define our variables and target region
    target = "chr2R"
    chromstart = 10000000
    chromend =  12000000
    chromlen = chromend - chromstart
    binsize = 200
    number_of_peaks = 200

    # Load in the bedgraph data
    forward = load_bedgraph(forward_fname, target, 0, chromlen)
    reverse = load_bedgraph(reverse_fname, target, 0, chromlen)

    # Combine forward and reverse tags
    combined = forward + reverse

    # Bin tag density using a sliding window
    scores = bin_array(combined, binsize)

    # Identify the top N peak positions
    peaks = find_peaks(scores, number_of_peaks, binsize)

    # Find combined tag densities over the selected peaks by strand
    reverse_curve = find_profile(reverse, peaks, binsize * 4)
    forward_curve = find_profile(forward, peaks, binsize * 4)

    # Find best shift to match forward and reverse peaks
    correlations = find_correlations(reverse_curve, forward_curve)

    # Plot results
    fig, ax = plt.subplots(2, 1, figsize=(5, 10))

    # Plot unshifted strand profiles
    X = numpy.arange(-reverse_curve.shape[0] // 2, reverse_curve.shape[0] // 2)
    ax[0].plot(X, reverse_curve, color='blue', label='reverse')
    ax[0].plot(X, forward_curve, color='red', label='foward')
    ax[0].set_xlabel("Distance from combined peak (bp)")
    ax[0].set_ylabel("Mean tag count")
    ax[0].set_title("Peak offset by strand")

    # Plot correlation curve
    ax[1].plot(numpy.arange(correlations.shape[0]), correlations, color='black')
    ax[1].set_xlabel("Fragment size")
    ax[1].set_ylabel("Correlation")
    ax[1].set_title('Correlation by offset')

    # Mark correlation peak
    best = numpy.argmax(correlations)
    bestval = numpy.amax(correlations)
    ax[1].plot([best, best], [0, bestval], color='black')
    plt.tight_layout()
    plt.savefig(out_fname)
    plt.close()

    # Report ideal shift for downstream analysis
    print("Best offset was {}".format(best))


def load_bedgraph(fname, target, chromstart, chromend):
    # Create array to hold tag counts
    coverage = numpy.zeros(chromend - chromstart, int)

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
    binned = numpy.zeros(data.shape[0], data.dtype)

    # For each position in the window, add to the score array
    for i in range(binsize):
        binned[i:data.shape[0] - binsize + i] += data[binsize//2:-binsize//2]
    return binned

def find_peaks(data, target, binsize):
    order = numpy.argsort(data)[::-1]
    return order[:target]

def find_profile(data, peaks, binsize):
    # Create array to hold the combined profile
    results = numpy.zeros(binsize, float)
    for i in range(peaks.shape[0]):
        # For each peak, get the tag density +/- half the binsize
        results += data[peaks[i]-binsize//2:peaks[i]+binsize//2]
    # Return average tag density
    return results / peaks.shape[0]

def find_correlations(reverse, forward):
    # Since we made the profiles 4 times bigger than our estimated fragment size
    # we will search a size from 0 to 2 times our estimate
    width = reverse.shape[0] // 2
    # Create an array to hold correlations
    corrs = numpy.zeros(width, float)
    for i in range(width):
        # For each shift, find correlation of overlapping regions 
        corrs[i] = numpy.corrcoef(reverse[i:], forward[:forward.shape[0]-i])[0, 1]
    return corrs
    # Load the sample bedgraph data, reusing the function we already wrote
    
    # Combine tag densities, shifting by our previously found fragment width

    # Load the control bedgraph data, reusing the function we already wrote
    
    # Combine tag densities
    
    # Adjust the control to have the same coverage as our sample

    # Create a background mean using our previous binning function and a 1K window
    # Make sure to adjust to be the mean expected per base

    # Find the mean tags/bp and make each background position the higher of the
    # the binned score and global background score

    # Score the sample using a binsize that is twice our fragment size
    # We can reuse the binning function we already wrote

    # Find the p-value for each position (you can pass a whole array of values
    # and and array of means). Use scipy.stats.poisson for the distribution.
    # Remeber that we're looking for the probability of seeing a value this large
    # or larger
    # Also, don't forget that your background is per base, while your sample is
    # per 2 * width bases. You'll need to adjust your background

    # Transform the p-values into -log10
    # You will also need to set a minimum pvalue so you doen't get a divide by
    # zero error. I suggest using 1e-250

    # Write p-values to a wiggle file
    # The file should start with the line
    # "fixedStep chrom=CHROM start=CHROMSTART step=1 span=1" where CHROM and
    # CHROMSTART are filled in from your target genomic region. Then you have
    # one value per line (in this case, representing a value for each basepair).
    # Note that wiggle files start coordinates at 1, not zero, so add 1 to your
    # chromstart. Also, the file should end in the suffix ".wig"

    # Write bed file with non-overlapping peaks defined by high-scoring regions 

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
    while numpy.amax(scores) >= 10:
        pos = numpy.argmax(scores)
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