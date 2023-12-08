#!/usr/bin/env python

import sys

import scanpy as sc
import numpy
import matplotlib.pyplot as plt


# def main():
#     sc.settings.verbosity = 3
#     sc.logging.print_header()
#     adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/',
#                             var_names='gene_symbols', cache=True)

#     adata.var_names_make_unique()

#     sc.tl.pca(adata, svd_solver='arpack')
#     raw = adata.copy()
#     #sc.pl.pca(adata, title='Unfiltered', save="_unfiltered.pdf")

#     print("# cells, # genes before filtering:", adata.shape)
#     sc.pp.filter_cells(adata, min_genes=200)
#     sc.pp.filter_genes(adata, min_cells=3)
#     print("# cells, # genes after filtering:", adata.shape)

#     adata.var['mt'] = adata.var_names.str.startswith('MT-')
#     sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
#                                log1p=False, inplace=True)
#     #sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#     #         jitter=0.4, multi_panel=True)
    
#     adata = adata[adata.obs.n_genes_by_counts < 2500, :]
#     adata = adata[adata.obs.pct_counts_mt < 5, :]
#     print("# cells, # genes after MT filtering:", adata.shape)
#     #sc.pl.highest_expr_genes(adata, n_top=20)

#     sc.pp.normalize_total(adata, target_sum=1e4)
#     sc.pp.log1p(adata)

#     sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3,
#                                 min_disp=0.5)
#     #sc.pl.highly_variable_genes(adata)
#     adata.write("filtered_data.h5")

#     adata.raw = adata
#     adata = adata[:, adata.var.highly_variable]
#     print("# cells, # genes after variability filtering:", adata.shape)

#     sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#     sc.pp.scale(adata, max_value=10)

#     sc.tl.pca(adata, svd_solver='arpack')
#     #sc.pl.pca_variance_ratio(adata, log=True)
#     #sc.pl.pca(adata, color='CST3')

#     #fig, ax = plt.subplots(1, 2, figsize=(10, 5))
#     #sc.pl.pca(raw, ax=ax[0], title="Uniltered", show=False)
#     #sc.pl.pca(adata, ax=ax[1], title="Filtered", show=False)
#     #plt.tight_layout()
#     #plt.savefig("pca.pdf")
#     #plt.close()

#     adata.write('variable_data.h5')

# main()


# Read the 10x dataset filtered down to just the highly-variable genes
adata = sc.read_h5ad("variable_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 

sc.pp.neighbors(adata, 10,40)
sc.tl.leiden(adata)
sc.tl.umap(adata,maxiter = 900)
sc.tl.tsne(adata)

# fig, axes = plt.subplots(ncols=2)
# sc.pl.umap(adata,color='leiden',ax = axes[0],title='UMAP - Leiden Clusters', show=False,legend_loc='right margin')
# sc.pl.tsne(adata, color='leiden', ax=axes[1], title='t-SNE - Leiden Clusters', show=False)

# plt.show()

# wilcoxon_adata = sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', use_raw=True, copy=True)

# logreg_adata = sc.tl.rank_genes_groups(adata, groupby='leiden', method='logreg', use_raw=True, copy=True)

# fig, ax = plt.subplots(figsize=(10, 8))
# sc.pl.rank_genes_groups(wilcoxon_adata, n_genes=25, ax=ax, title='Wilcoxon Rank-Sum', sharey=False, show=False, use_raw=True)
# plt.show()

# fig, ax = plt.subplots(figsize=(10, 8))
# sc.pl.rank_genes_groups(logreg_adata, n_genes=25, ax=ax, title='Logistic Regression', sharey=False, show=False, use_raw=True)

# plt.show()

#3
leiden = adata.obs['leiden']
umap = adata.obsm['X_umap']
tsne = adata.obsm['X_tsne']
adata = sc.read_h5ad('filtered_data.h5')
adata.obs['leiden'] = leiden
adata.obsm['X_umap'] = umap
adata.obsm['X_tsne'] = tsne

adata.write('filtered_clustered_data.h5')

adata = sc.read_h5ad("filtered_clustered_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 

#marker_gene_list={'cluster0':['RPS12','LDHB','RPS25'],'cluster4':['COTL1','FCER1G','LST1'],'cluster3':['CCL5','NKG7','B2M']}

name_list=['T-cell','Myeloid','B-cell','T-cell ','Myeloid ','T-cell  ',6,7]

# fig, axes = plt.subplots(ncols=3)
# sc.pl.umap(adata,color='MS4A1',ax = axes[0],title='MS4A1', show=False)
# sc.pl.umap(adata,color='CST3',ax = axes[1],title='CST3', show=False)
# sc.pl.umap(adata,color='CD3D',ax = axes[2],title='CD3D', show=False)
# plt.show()


adata.rename_categories('leiden',name_list)

fig, axes = plt.subplots(ncols=2)
sc.pl.umap(adata,color='leiden',ax = axes[0],title='UMAP - Leiden Clusters', show=False,legend_loc='on data')
sc.pl.tsne(adata, color='leiden', ax=axes[1], title='t-SNE - Leiden Clusters', show=False,legend_loc = 'on data')

plt.show()