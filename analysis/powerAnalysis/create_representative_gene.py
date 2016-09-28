import copy
import re
import StringIO

import numpy as np
import pandas as pd
import yaml

# Global constants
N_CHROM = 6

# Set default arguments

# Config file to use
path_cfg = 'config/H_1-combined.yml'

# Prefix for data files
prefix_data = '../XuData/'

# Index of ORFs
path_gene_index = '../geneInfo/geneIndex.txt'

# Basepairs per cluster
bp_per_cluster = 147

# Length to extract as promoter
length_promoter = 1000

# Length of ORF in representative gene
length_representative_orf = 2500

# Anchor nucleosome for representative gene construction
# This cluster's center is computed directly from the observed ORFs
# Other clusters are positioned based upon observed spacings
anchor = 1

# Output filename
path_output = 'powerAnalysis/data/representative_gene.txt'


def clean_config(cfg):
    '''
    Clean configuration dictionary parsed from YAML by stripping newlines and
    whitespace.
    '''
    cfg = copy.copy(cfg)
    for k, v in cfg.iteritems():
        if type(v) is str:
            cfg[k] = v.strip()
        elif type(v) is dict:
            cfg[k] = clean_config(v)

    return(cfg)


def parse_gene_index(gene_index_file, delimiter='\t'):
    '''
    Parse gene index file to numpy record array with clean format
    Each row becomes a record
    '''
    # Setup input file if needed
    if type(gene_index_file) is not file:
        gene_index_file = open(gene_index_file, "rb")
 
    # Read all text as string
    gene_index_txt = gene_index_file.read()
    gene_index_file.close()
    
    # Strip out quotes (np.genfromtxt does not remove them)
    gene_index_txt = re.sub('"|\'', '', gene_index_txt)
    
    # Setup record array for gene info
    gene_index = np.recfromtxt(StringIO.StringIO(gene_index_txt),
                               delimiter=delimiter, names=True)
    
    return gene_index


def gene_to_slice(gene, promoter_length=600):
    '''
    Obtain slice for extraction from chromosome based on gene record.
    Record must have 'start' and 'stop' fields.
    '''
    if gene['start'] < gene['stop']:
        return slice(gene['start']-promoter_length, gene['stop']+1, 1)
    else:
        return slice(gene['start']+promoter_length, gene['stop']-1, -1)


def compute_ranks(x):
    '''
    Compute ranks of cluster centers relative to TSS
    '''
    ranks = np.argsort(x) - np.sum(x<=0)
    if np.all(x != 0):
        ranks[ranks >= 0] += 1
    return ranks


def compute_prop(x):
    '''
    Convert x to proportion; x -> x/sum(x)
    '''
    return x / np.sum(x)


def main():
    # Load configuration
    with open(path_cfg, 'rb') as f:
        cfg = yaml.load(f)
        cfg = clean_config(cfg)

    # Load posterior summaries
    clusters = []
    pattern_clusters = cfg['mcmc_output']['cluster_pattern'].format(**cfg)
    for chrom in xrange(1, N_CHROM+1, 1):
        with open(pattern_clusters % chrom, 'rb') as f:
            data = pd.read_table(f, sep=' ', header=0)
            data['chrom'] = chrom
            clusters.append(data)

    clusters = pd.concat(clusters)

    # Load and parse gene index
    with open(path_gene_index, 'rb') as f:
        gene_index = parse_gene_index(f)
    
    # Label each ORF and its promoter
    orfs = {}
    for gene in gene_index:
        # Skip missing chromosomes; primarily for testing
        if gene['chromosome'] > len(clusters):
            continue

        gene_slice = gene_to_slice(gene, promoter_length=length_promoter)
        systematic = gene['systematic']
        chrom = gene['chromosome']
        
        # Find clusters to extract
        query = ((clusters.chrom == chrom) &
                 (clusters.center >= gene_slice.start) & 
                 (clusters.center < gene_slice.stop))

        # Extract and label clusters
        orfs[systematic] = clusters[query].copy()
        orfs[systematic]['orf'] = systematic
        orfs[systematic]['orf_length'] = gene_slice.stop - gene_slice.start

        # Make cluster centers relative to TSS
        orfs[systematic]['center'] -= gene['start']
        if gene['stop'] < gene['start']:
            orfs[systematic]['center'] *= -1

    # Combine ORFs into a single DataFrame
    orfs = pd.concat(orfs)

    # Append ranks by ORF
    grouped = orfs.groupby('orf')
    orfs['rank'] = grouped['center'].transform(compute_ranks)
    
    # Append relative occupancy by ORF
    orfs['relative_occupancy'] = grouped['occupancy'].transform(compute_prop)
    orfs['relative_occupancy'] *= (orfs['orf_length'] /
                                   np.mean(orfs['orf_length'].unique()))

    # Reindex orfs dataframe
    orfs.set_index(['orf', 'rank'], drop=False, inplace=True)

    # Append spacing between clusters by ORF
    grouped = orfs.groupby('orf')
    orfs['spacing'] = grouped['center'].diff()

    # Summarize cluster properties by cluster order from TSS
    grouped = orfs[orfs['rank'] != 0].groupby('rank')

    spacing_by_rank = grouped['spacing'].agg({'mean' : np.mean, 'std' : np.std})
    center_by_rank = grouped['center'].agg({'mean' : np.mean, 'std' : np.std})
    ro_by_rank = grouped['relative_occupancy'].agg({'mean' : np.mean,
                                                    'std' : np.std})
    structure_by_rank = grouped['structure'].agg({'mean' : np.mean,
                                                  'std' : np.std})

    # Determine cluster centers for representative gene
    rep_centers = pd.Series(np.zeros(center_by_rank.shape[0]-1),
                            index=center_by_rank.index[1:])

    # Find anchor in rep_centers index
    ind_anchor = np.searchsorted(rep_centers.index, anchor, side='left')

    # Start with anchor
    rep_centers += center_by_rank['mean'][anchor]

    # Set remaining centers based upon spacing
    adj = spacing_by_rank['mean'][ind_anchor+1:1:-1].cumsum()[::-1]
    rep_centers[:ind_anchor] -= adj.values

    adj = spacing_by_rank['mean'][ind_anchor+2:].cumsum()
    rep_centers[ind_anchor+1:] += adj

    # Restrict to desired range of distances from TSS
    rep_centers = rep_centers[(rep_centers > -length_promoter) &
                              (rep_centers <= length_representative_orf)]

    # Build dataframe of relevant quantities
    rep_gene = pd.DataFrame({'center' : rep_centers,
        'relative_occupancy' : ro_by_rank['mean'][rep_centers.index]})
    rep_gene['relative_occupancy'] /= rep_gene['relative_occupancy'].sum()
    rep_gene['structure'] = structure_by_rank['mean'][rep_centers.index]
    
    # Save representative gene
    rep_gene.to_csv(path_output)

    return rep_gene


if __name__ == '__main__':
    rep_gene = main()

