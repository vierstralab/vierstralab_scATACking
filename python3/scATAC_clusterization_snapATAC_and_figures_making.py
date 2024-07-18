import sys
import pandas as pd
import numpy as np
from scipy.sparse import csr_array
import anndata as ad
from matplotlib import pyplot as plt
from sklearn import preprocessing
import snapatac2 as snap
from scipy.sparse import csr_matrix
from matplotlib import colors
from scipy.sparse import load_npz
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

ids = str(sys.argv[1])
counts_file = sys.argv[2]
input_file = sys.argv[3]
barcodes_map = sys.argv[4]
index_mapping = sys.argv[5]
fragments_numbers = [2000,4000]


meta_data = pd.read_csv(counts_file, sep = '\t', header = None)
meta_data['ID'] = meta_data[0].apply(lambda x: x[:x.find('_')])
meta_data = meta_data.rename(columns={1:'#_of_detectesd_scATAC_fragments'})

data = load_npz(input_file)
chunks_classes = np.loadtxt(index_mapping, dtype = 'str')
cells = np.loadtxt(barcodes_map, dtype = 'str')

adata = ad.AnnData(data)
adata.var_names = cells
adata.obs_names = chunks_classes
data_csr = csr_matrix(adata.X)
adata_csr = ad.AnnData(data_csr)
adata_csr.var_names = cells
adata_csr.obs_names = chunks_classes

for i in fragments_numbers:
   adata_csr =  adata_csr[:,meta_data.loc[meta_data['#_of_detectesd_scATAC_fragments'] >= i, 0]]
   adata_csr_transponated = adata_csr.copy().T
   snap.pp.scrublet(adata_csr_transponated, features = None, n_jobs=20)
   snap.pp.filter_doublets(adata_csr_transponated)
   adata_selection = snap.pp.select_features(adata_csr_transponated, filter_upper_quantile=0.05, filter_lower_quantile=0.05, inplace = False)
   adata_csr_selected = adata_csr_transponated[:,adata_selection]
   snap.pp.select_features(adata_csr_selected, filter_upper_quantile=0.05, filter_lower_quantile=0.05)
   snap.tl.spectral(adata_csr_selected)
   snap.tl.umap(adata_csr_selected)
   snap.pp.knn(adata_csr_selected)
   snap.tl.leiden(adata_csr_selected, resolution = 0.5)
   snap.pl.umap(adata_csr_selected, color='leiden', interactive=False, height=750, width=1000, out_file=ids+'_'+str(i)+'_fragments_leiden_plot.pdf')

   adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] = meta_data.loc[meta_data[0].isin(adata_csr_selected.obs.index), '#_of_detectesd_scATAC_fragments'].values
   adata_csr_selected.obs['_of_detectesd_scATAC_fragments_groups'] = ''
   adata_csr_selected.obs.loc[adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] <= 3000, '#_of_detectesd_scATAC_fragments_groups'] = '<= 3000 fragments per cell'
   adata_csr_selected.obs.loc[(adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] <= 5000)&(adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] > 3000), '#_of_detectesd_scATAC_fragments_groups'] = '> 3000 and <= 5000 fragments per cell'
   adata_csr_selected.obs.loc[(adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] <= 7500)&(adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] > 5000), '#_of_detectesd_scATAC_fragments_groups'] = '> 5000 and <= 7500 fragments per cell'
   adata_csr_selected.obs.loc[(adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] <= 10000)&(adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] > 7500), '#_of_detectesd_scATAC_fragments_groups'] = '> 7500 and <= 10000 fragments per cell'
   adata_csr_selected.obs.loc[(adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] <= 20000)&(adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] > 10000), '#_of_detectesd_scATAC_fragments_groups'] = '> 10000 and <= 20000 fragments per cell'
   adata_csr_selected.obs.loc[adata_csr_selected.obs['#_of_detectesd_scATAC_fragments'] > 20000, '#_of_detectesd_scATAC_fragments_groups'] = '> 20000 fragments per cell'
   snap.pl.umap(adata_csr_selected, color='_of_detectesd_scATAC_fragments_groups', interactive=False, height=750, width=1200, out_file=ids+'_'+str(i)+'_fragments_number_of_fragments_plot.pdf')

fig, ax = plt.subplots()
ax.hist(meta_data.loc[meta_data['#_of_detectesd_scATAC_fragments'] >= 2000,'#_of_detectesd_scATAC_fragments'], bins = 100, orientation = "horizontal")
ax.set_xscale('log')
plt.ylabel('Number of detectesd scATAC fragments')
plt.xlabel('Number of cells in log scale')
plt.savefig(ids+'_cells_more_than_2000_fragments_hist_log.pdf')

fig, ax = plt.subplots()
ax.hist(meta_data.loc[(meta_data['#_of_detectesd_scATAC_fragments'] >= 2000)&(meta_data['#_of_detectesd_scATAC_fragments'] <= 10000),'#_of_detectesd_scATAC_fragments'], bins = 100, orientation = "horizontal")
ax.set_xscale('log')
plt.ylabel('Number of detectesd scATAC fragments')
plt.xlabel('Number of cells in log scale')
plt.savefig(ids+'_cells_more_than_2000_and_less_than_10000_fragments_hist_log.pdf')

fig, ax = plt.subplots()
ax.hist(meta_data.loc[(meta_data['#_of_detectesd_scATAC_fragments'] >= 2000)&(meta_data['#_of_detectesd_scATAC_fragments'] <= 10000),'#_of_detectesd_scATAC_fragments'], bins = 100, orientation = "horizontal")
plt.ylabel('Number of detectesd scATAC fragments')
plt.xlabel('Number of cells')
plt.savefig(ids+'_cells_more_than_2000_and_less_than_10000_fragments_hist.pdf')

