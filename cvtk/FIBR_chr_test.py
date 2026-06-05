#!/usr/bin/env python
# coding: utf-8

###### Working python script to run cvtk on FIBR temporal data with minimal changes

# ## Setup
import os
import sys
from collections import Counter, defaultdict

import pandas as pd
import numpy as np

# Ensure the path to cvtk is maintained
sys.path.append("/lustre/home/bf299/cvtk")

from cvtk.variant_files import VCFFile
from cvtk.gintervals import GenomicIntervals
from cvtk.cvtk import TiledTemporalFreqs

# Argument parser for command-line inputs
import argparse
def parse_arguments():
    parser = argparse.ArgumentParser(description="Run CVTK covariance analysis on FIBR data.")
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="Path to the input VCF file (e.g., merged_output_unphase.vcf.gz)."
    )
    parser.add_argument(
        "--popmap",
        type=str,
        required=True,
        help="Path to the population map file (e.g., merged_popmap_int.txt)."
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path for the output covariance file (e.g., gw_covs_unphased.tsv)."
    )
    return parser.parse_args()

def main():
    # Parse command-line arguments
    args = parse_arguments()
    VCFFilePath = args.vcf
    popmapPath = args.popmap
    outputPath = args.output

    # Load the VCF data
    print(f"Loading VCF file: {VCFFilePath}")
    vcf = VCFFile(VCFFilePath)

    # Load sample/population data from the population map file
    print(f"Loading population map: {popmapPath}")
    samples = pd.read_csv(popmapPath, header=None, names=('pop', 'individual'))
    sample_map = {k: v for k, v in zip(samples['individual'], samples['pop'])}

    # Map individual IDs to subpopulation indices
    subpop_indices = defaultdict(list)
    for i, k in enumerate(vcf.samples):
        subpop_indices[sample_map[k.decode()]].append(i)

    # Calculate allele counts for subpopulations
    counts_mat = vcf.geno_mat.count_alleles_subpops(subpop_indices)
    vcf.mat = np.stack(list(counts_mat.values()))
    vcf.subpops = subpop_indices.keys()
    
    # Count diploids in each sample
    ndiploids = [Counter(sample_map.values())[k] for k in vcf.subpops]

    # Parse design matrix: (replicate, generation)
    def parse_samples(x):
        rep, gen = x.split('_')
        gen = int(gen)
        rep = int(rep)
        return (rep, gen)

    design = [parse_samples(x) for x in vcf.subpops]

    # Calculate frequencies
    freq_mat_all = vcf.calc_freqs()
    print("number of loci: ", freq_mat_all.shape[1])

    # Filter out non-segregating sites
    vcf.remove_fixed()
    freq_mat = vcf.calc_freqs()
    print("number of loci (after filtering): ", freq_mat.shape[1])
    print("loci not segregating removed: ", freq_mat_all.shape[1] - freq_mat.shape[1])

    # Define genomic intervals and divide into tiles
    tile_width = 1e6  # 1 Mb tile width
    gi = vcf.build_gintervals()
    gi.infer_seqlens()
    tiles = GenomicIntervals.from_tiles(gi.seqlens, width=tile_width)

    # Create TiledTemporalFreqs object
    d = TiledTemporalFreqs(
        tiles, freqs=freq_mat, depths=vcf.N, diploids=ndiploids,
        samples=design, gintervals=gi
    )

    # Print data dimensions for verification
    #print("d.freqs.shape: ", d.freqs.shape)
    #print("d.samples:", d.samples)
    #print("d.depths:", d.depths.shape)

   #ave_depth = d.depths.mean(axis=(0,1))
    #print("ave_depth: min, mean, median, max =", ave_depth.min(), ave_depth.mean(), np.median(ave_depth), ave_depth.max())

    #vcf.N seems to be fine
    #arr = np.asarray(vcf.N)
    #print("type:", type(vcf), "  vcf.N type:", type(vcf.N))
    #print("vcf.N shape:", arr.shape)
    #print("vcf.N dtype:", arr.dtype)
    #print("first 20 values:", arr.flatten()[:20])
    #print("unique values (first 20):", np.unique(arr)[:20])


    # Calculate genome-wide covariance matrix
    gw_covs = d.calc_cov()
    #print("gw_shape:", gw_covs.shape)

    # Build labels for rows/columns of covariance matrix
    R, T = d.R, d.T
    labels = []
    for r in range(R):
        base = r * d.ntimepoints
        for t_idx in range(T):
            before = d.samples[base + t_idx]
            after = d.samples[base + t_idx + 1]
            labels.append(f"rep{r}-delta-{before}-{after}")
    
    # Write genome-wide covariance matrix to output file
    header = '\t'.join(labels)
    #print(f"Saving genome-wide covariance matrix to: {outputPath}")
    #np.savetxt(outputPath, gw_covs, delimiter='\t', header=header, comments='', fmt='%.8g')

    covs_cis = d.bootstrap_cov(B=5000, average_replicates=False, progress_bar=False)
    print("cov_cis shape:", covs_cis.shape)

    #this produces a 3D matrices an upper, lower, and mean covs_cis
    mean_cov = covs_cis[1]
    filename_mean = f"{outputPath}_gw_covs_mean_unphased.tsv"
    np.savetxt(filename_mean, mean_cov, delimiter='\t', header=header, comments='', fmt='%.8g')
    low_cov = covs_cis[0]
    filename_low = f"{outputPath}_gw_covs_low_unphased.tsv"
    np.savetxt(filename_low, low_cov, delimiter='\t', header=header, comments='', fmt='%.8g')
    high_cov = covs_cis[2]
    filename_high = f"{outputPath}_gw_covs_high_unphased.tsv"
    np.savetxt(filename_high, high_cov, delimiter='\t', header=header, comments='', fmt='%.8g')


if __name__ == "__main__":
    main()
