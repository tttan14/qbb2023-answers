#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from model_peaks import load_bedgraph, bin_array

def parse_bedgraph(file_path):
    columns = ["chrom", "start", "end", "methylation_score","coverage"]
    df = pd.read_csv(file_path, sep='\t', header=None, names=columns)
    return df

def compare(file1,file2):
    #bisulfite
    pos1=file1.loc[:,'start']
    pos2=file2.loc[:,'start']
    common_loc=pos1[pos1.isin(pos2)]

    diff=[]
    a=0
    b=1
    for i in common_loc:
        #print(i)
        ind1=int(pos1.index.get_loc(pos1[pos1 ==i].index[0]))
        ind2=int(pos2.index.get_loc(pos2[pos2 ==i].index[0]))
        #print(ind1)
        if file2.loc[ind2,'methylation_score'] != file1.loc[ind1,'methylation_score']:
            a=float(file2.loc[ind2,'methylation_score'])-float(file1.loc[ind1,'methylation_score'])
            diff.append(a)
        #print(diff)
        a=0
    diff=np.array(diff)
    return diff



def main():
    bis_file=parse_bedgraph("bisulfite.cpg.chr2.bedgraph")
    nanopore_file=parse_bedgraph("ONT.cpg.chr2.bedgraph")
    #unique,shared=compare_stats(bis_file,nanopore_file)
    pos_bis=bis_file.loc[:,'start']
    pos_nano=nanopore_file.loc[:,'start']
    common_val=pos_bis[pos_bis.isin(pos_nano)]
    unique_bis=pos_bis[~pos_bis.isin(common_val)]
    unique_nano=pos_nano[~pos_nano.isin(common_val)]  
    #plt.violinplot(common_val, positions=[1], widths=0.8, showmedians=True, showextrema=False)

    print('unique site number in bisulfite: ',len(unique_bis))
    print('unique site number in bisulfite percentage: ',(len(unique_bis)/len(pos_bis))*100,'%')
    print('unique site number in nanopore: ',len(unique_nano))
    print('unique site number in nanopore percentage: ',(len(unique_nano)/len(pos_nano))*100,'%')
    print('shared site number: ',len(common_val))
    print('shared site number as percentage with bisulfite: ',(len(common_val)/len(pos_bis))*100,'%')
    print('shared site number as percentage with nanopore: ',(len(common_val)/len(pos_nano))*100,'%')

    #Plot
    fig,ax = plt.subplots(nrows=1, ncols=3)
    ax[0].hist(bis_file.loc[:,'coverage'],bins=100,alpha=0.5,color='orange',label='bisulfite')
    ax[0].hist(nanopore_file.loc[:,'coverage'],bins=100,alpha=0.5,color='blue',label='nanopore')
    ax[0].set_xlabel("coverage")
    ax[0].set_ylabel("counts")
    ax[0].set_title('bisulfite vs nanopore coverage') 
    ax[0].legend()

    bis_score=bis_file.loc[:,'methylation_score']
    nano_score=nanopore_file.loc[:,'methylation_score']
    min_length = min(len(nano_score), len(bis_score))
    nano_scorepad = np.pad(nano_score, (0, max(0, len(bis_score) - len(nano_score))), mode='constant')
    hist, xedges, yedges = np.histogram2d(bis_score, nano_scorepad, bins=100)
    adj_hist=np.log10(hist+1)
    im = ax[1].imshow(adj_hist.T, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], origin='lower', cmap='viridis')
    cbar = plt.colorbar(im)
    cbar.set_label('Frequency')
    ax[1].set_xlabel('bisulfite')
    ax[1].set_ylabel('nanopore')

    p_corr = np.corrcoef(bis_score, nano_scorepad)[0, 1]
    ax[1].set_title(f'bisulfite vs nanopore methylation score (Pearson R: {p_corr:.3f})')
    
    nornano_file=parse_bedgraph("normal.ONT.chr2.bedgraph")
    tunano_file=parse_bedgraph("tumor.ONT.chr2.bedgraph")
    norbis_file=parse_bedgraph("normal.bisulfite.chr2.bedgraph")
    tubis_file=parse_bedgraph("tumor.bisulfite.chr2.bedgraph")

    diffnano=compare(norbis_file,tubis_file)
    diffbis=compare(nornano_file,tunano_file)

    ax[2].violinplot(diffnano, positions=[1], widths=0.8, showmedians=True, showextrema=False)
    ax[2].violinplot(diffbis, positions=[2], widths=0.8, showmedians=True, showextrema=False)
    ax[2].set_xticks([1, 2])
    ax[2].set_xticklabels(['Nanopore', 'Bisulfite'])
    ax[2].set_xlabel('Approach')
    ax[2].set_ylabel('Methylation Change')
    corr_diff = np.corrcoef(diffnano, diffbis)[0, 1]
    ax[2].set_title('bisulfite vs nanopore methylation changes (tumor/normal)\nPearson R: {:.2f}'.format(corr_diff))
    plt.show()








if __name__ == "__main__":
    main()

    # bismark_only, nanopore_only, shared, shared_percentage = calculate_comparison_stats(bismark_file_path, nanopore_file_path)

    # print(f"only in normal file: {bismark_only}")
    # print(f"only in nanopore file: {nanopore_only}")
    # print(f"shared sites: {shared}")
    # print(f"% shared sites: {shared_percentage:.2f}%")
