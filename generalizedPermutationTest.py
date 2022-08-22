### Rauf Salamzade; Kalan Lab; MMI UW-Madison
### This script performs a mean-based permutation test for looking at whether 
### homolog groups are enriched/depleted in copy-count/presence in a focal
### vs. comparator/other genome set. The format of the first two arguments
### should be a genome identifier per line (should match the header of the
### OrthoFinder Orthogroups.tsv result file.

import os
import sys
import statistics
import random
from scipy import stats 
import numpy as np

def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)

def permute(focal_vals, other_vals, S=100000):
    lean = 'NA'
    if statistics.mean(focal_vals) < statistics.mean(other_vals):
        lean = 'Other'
    elif statistics.mean(focal_vals) > statistics.mean(other_vals):
        lean = 'Focal'
    res = stats.permutation_test((np.array(focal_vals), np.array(other_vals)), statistic, vectorized=True, n_resamples=S)
    stat = res.statistic
    emp_pval = res.pvalue
    return([emp_pval, lean])

if len(sys.argv) < 2:
    print('Usage: python generalizedPermutationTest.py <focal_clade_listing_file> <comparator_clade_listing_file> <orthofinder_orthogroups.tsv_file>')
    sys.exit(1)
focal_listing_file = sys.argv[1]
other_listing_file = sys.argv[2] # 
orthogroups_tsv = sys.argv[3] # Orthogroups.tsv file from OrthoFinder2 

focal_clade = set([x.strip() for x in open(focal_listing_file).readlines()])
other_clade = set([x.strip() for x in open(other_listing_file).readlines()])

focal_inds = set([])
other_inds = set([])
print('\t'.join(['homolog_group', 'empirical_pvalue', 'lean', 'clade1_proportion_with_og', 'clade2_proportion_with_og', 'clade1_median_copy_count', 'clade2_median_copy_count']))
with open(orthogroups_tsv) as otf:
    for i, line in enumerate(otf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: 
            tot_samps = 0
            for j, val in enumerate(ls[1:]):
                if val in focal_clade: 
                    focal_inds.add(j)
                elif val in other_clade:
                    other_inds.add(j)
                tot_samps += 1
        else:
            og = ls[0]
            foc_samp_genes = []
            focal_copy_counts = []
            other_copy_counts = []
            for j, val in enumerate(ls[1:]):
                gene_count = len([x for x in val.split(', ') if x.strip() != ''])
                if j in focal_inds:
                    focal_copy_counts.append(gene_count)
                elif j in other_inds:
                    other_copy_counts.append(gene_count)
            
            focal_samples_with_og = sum([1 for x in focal_copy_counts if x > 0])
            other_samples_with_og = sum([1 for x in other_copy_counts if x > 0])
            focal_proportion_with_og = float(focal_samples_with_og)/(len(focal_inds))
            other_proportion_with_og = float(other_samples_with_og)/(len(other_inds))
            focal_copy_count = "NA"
            other_copy_count = "NA"
            if sum(focal_copy_counts) > 0:
                focal_copy_count = statistics.median([x for x in focal_copy_counts if x > 0])
            if sum(other_copy_counts) > 0:
                other_copy_count = statistics.median([x for x in other_copy_counts if x > 0])
            emp_pval, lean = permute(focal_copy_counts, other_copy_counts)
            print('\t'.join([str(x) for x in [og, emp_pval, lean, focal_proportion_with_og, other_proportion_with_og, focal_copy_count, other_copy_count]]))
