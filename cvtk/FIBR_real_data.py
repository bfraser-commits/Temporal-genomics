#!/usr/bin/env python
# coding: utf-8

###### Working python script to run cvtk on FIBR temporal data

# ## Setup
import os
import sys
import glob
nb_dir = os.path.split(os.getcwd())[0]
if nb_dir not in sys.path:
    sys.path.append(nb_dir)

import re
import pickle
from collections import Counter
from functools import reduce

import pandas as pd
import scipy
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl


import sys
sys.path.append("/lustre/home/bf299/cvtk")

from cvtk.cvtk import TemporalFreqs, TiledTemporalFreqs
from cvtk.cov import stack_temporal_covariances
import cvtk.variant_files as vf
from cvtk.gintervals import GenomicIntervals
from cvtk.pca import FreqPCA
from cvtk.plots import rep_plot_pca, correction_diagnostic_plot
from cvtk.utils import integerize
from cvtk.utils import extract_empirical_nulls_diagonals, extract_temporal_cov_diagonals
from cvtk.cov import stack_replicate_covariances, stack_temporal_covs_by_group
from cvtk.variant_files import VCFFile


# # Varianta Data Loading

# ### Load in VCF data


vcf=VCFFile('../merged_output_unphase.vcf.gz')

# Remove fixed sites — those that are not polymorphic in any samples / timepoints. These just needlessly shrink the covariance towards zero.
# this doesn't change the number of SNPs, not really necessary

# ### Sample Data
#change this so that it's now rep and gen and they are integers
#using pop and individual
#CA=1 TA=2 LL=3 UL=4
samples = pd.read_csv("merged_popmap_int.txt", header=None, names = ('pop', 'individual'))
sample_map = {k:v for k, v in zip(samples['individual'], samples['pop'])}

subpop_indices = defaultdict(list)
for i, k in enumerate(vcf.samples):
    subpop_indices[sample_map[k.decode()]].append(i)


# From this, we can map the `vcf.geno_mat` table to subpopulation counts. 

counts_mat=vcf.geno_mat.count_alleles_subpops(subpop_indices)
vcf.mat = np.stack(list(counts_mat.values()))
vcf.subpops = subpop_indices.keys()

# Now we count the number of diploids in each sample.
ndiploids = [Counter(sample_map.values())[k] for k in vcf.subpops]


#parse IDs so that they fit the time and replicate format 
def parse_samples(x):
    rep, gen = x.split('_')
    gen = int(gen)
    rep = int(rep)
    return (rep, gen)

design = [parse_samples(x) for x in vcf.subpops]

freq_mat_all = vcf.calc_freqs()

print("number of loci: ", freq_mat_all.shape[1])

# With the frequencies calculated, now we filter out all non-segregating sites.
vcf.remove_fixed()
freq_mat = vcf.calc_freqs()
print("number of loci: ", freq_mat.shape[1])
print("loci not segregating removed: ", freq_mat_all.shape[1] - freq_mat.shape[1])

#this doesn't change number of SNPs so could remove but will move forward

# ## Replicate Covariance Analysis

#should we use 50kb or 1 Mb
tile_width = 1e6
gi = vcf.build_gintervals()
gi.infer_seqlens()


tiles = GenomicIntervals.from_tiles(gi.seqlens, width=tile_width)


#TemporalFreqs and TiledTemporalFreqs take the desing as a tuple (timepoint, replicate)
d = TiledTemporalFreqs(tiles, freqs=freq_mat, depths=vcf.N, diploids=ndiploids, samples=design, gintervals=gi)

#print things out to see if it's loaded all the data properly
#d.freqs.shape
print("d.freqs.shape: ", d.freqs.shape)

#d.samples
print("d.samples:",  d.samples)

#depths
print("d.depths:", d.depths.shape)

##genome-wide covariances
#first just a straight full covariance matrix
gw_covs= d.calc_cov()

print("gw_shape:" , gw_covs.shape)

# build labels for each row/column: replicate-major ordering
R, T = d.R, d.T
labels = []
for r in range(R):
    base = r * d.ntimepoints
    for t_idx in range(T):
        before = d.samples[base + t_idx]
        after  = d.samples[base + t_idx + 1]
        labels.append(f"rep{r}-delta-{before}-{after}")

header = '\t'.join(labels)
np.savetxt('gw_covs_unphased.tsv', gw_covs, delimiter='\t', header=header, comments='', fmt='%.8g')

##we could try lader to remove LG12 or LGs with significant selection
#autosomes = list(set(gi.intervals.keys()) - set('LG12'))

## bootstrap covariance matrix
#just set the progress_bar to false
#set average replicates to false bootstraps the entire-replicate covariance matrix
#the off diaganol is the covariances
covs_cis = d.bootstrap_cov(B=5000, average_replicates=False, progress_bar=False)
print("cov_cis shape:", covs_cis.shape)

#this produces a 3D matrices an upper, lower, and mean covs_cis
mean_cov = covs_cis[1]  
np.savetxt('gw_covs_mean_unphased.tsv', mean_cov, delimiter='\t', header=header, comments='', fmt='%.8g')
low_cov = covs_cis[0]
np.savetxt('gw_covs_low_unphased.tsv', low_cov, delimiter='\t', header=header, comments='', fmt='%.8g')
high_cov = covs_cis[2]
np.savetxt('gw_covs_hi_unphased.tsv', high_cov, delimiter='\t', header=header, comments='', fmt='%.8g')


#do those supplemental figures with to test bias correction

#calculate the temporal-replicate covariances per tile
tile_covs = d.calc_cov_by_tile(progress_bar=False)

#reshape the matrices extracting the temporal replicates to get temporal replicates
#this creates figs  
tile_temp_covs = stack_temporal_covs_by_group(tile_covs, d.R, d.T)

print("shape of tile_temp_covs:", tile_temp_covs.shape)

#next is the bias correction plots
diagnostics = d.correction_diagnostics()

df, models, xpreds, ypreds = diagnostics

# Save full diagnostics table
df.to_csv('correction_diagnostics_uphased.csv', index=False, na_rep='NA', float_format='%.6g')

