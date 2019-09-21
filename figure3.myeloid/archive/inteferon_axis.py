"""
	interferon_axis.py
	Contains script for computing distances between the monocytes in the lupus dataset
	to monocytes in the IFN beta stim dataset.
"""


import scanpy.api as sc
import numpy as np
import scipy.spatial.distance as dist

if __name__ == '__main__':

	lupus_mono_adata = sc.read('/ye/yelabstore3/Richard/PlottingScripts/ez_scanpy/V6/CLUESImmVarMonoDC.V6.combat.refined.h5ad')
	ifn_adata = sc.read('/netapp/home/mincheol/parameter_estimation/inteferon_data/interferon.norm.lupus.h5ad')
	ifn_mono_adata = ifn_adata[(ifn_adata.obs['stim'] == 'stim') & ifn_adata.obs['cell'].isin(['CD14+ Monocytes', 'FCGR3A+ Monocytes']), :]

	classical_mc = ifn_mono_adata[ifn_mono_adata.obs.cell == 'CD14+ Monocytes', :].X.mean(axis=0)
	nonclassical_mc = ifn_mono_adata[ifn_mono_adata.obs.cell == 'FCGR3A+ Monocytes', :].X.mean(axis=0)
	mean_signal = np.vstack([classical_mc, nonclassical_mc])

	# Compute transcriptomic correlations
	corr = dist.cdist(lupus_mono_adata.X, mean_signal, metric='correlation')
	np.save('/ye/yelabstore3/mincheol/lupus_ifn_corr_mean.npy', corr)

	# Compute Euclidean distances
	l2_dist = dist.cdist(lupus_mono_adata.X, mean_signal)
	np.save('/ye/yelabstore3/mincheol/lupus_ifn_l2_mean.npy', l2_dist)

	# # Compute transcriptomic correlations
	# corr = dist.cdist(lupus_mono_adata.X, ifn_mono_adata.X, metric='correlation')
	# np.save('/ye/yelabstore3/mincheol/lupus_ifn_corr.npy', corr)

	# # Compute Euclidean distances
	# l2_dist = dist.cdist(lupus_mono_adata.X, ifn_mono_adata.X)
	# np.save('/ye/yelabstore3/mincheol/lupus_ifn_l2.npy', l2_dist)
