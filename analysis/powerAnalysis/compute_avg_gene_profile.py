import copy
import re
import StringIO

import numpy as np
import yaml

# Global constants
N_CHROM = 16

# Set default arguments

# Config file to use
path_cfg = 'config/H_1-combined.yml'

# Prefix for data files
prefix_data = '../XuData/'

# Index of ORFs
path_gene_index = '../geneInfo/geneIndex.txt'

# Length to extract as promoter
length_promoter = 1000

# Length of ORF in representative gene
length_representative_orf = 2500

# Output filename
path_output = 'powerAnalysis/data/avg_gene_profile.txt'


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


def main():
    # Load configuration
    with open(path_cfg, 'rb') as f:
        cfg = yaml.load(f)
        cfg = clean_config(cfg)

    # Load posterior summaries
    summaries = []
    pattern_summaries = cfg['mcmc_output']['summary_pattern'].format(**cfg)
    for chrom in xrange(1, N_CHROM+1, 1):
        with open(pattern_summaries % chrom, 'rb') as f:
            summaries.append(np.genfromtxt(f, names=True))

    # Load and parse gene index
    with open(path_gene_index, 'rb') as f:
        gene_index = parse_gene_index(f)
    
    # Extract each ORF with its promoter
    orfs = {}
    for gene in gene_index:
        # Skip missing chromosomes; primarily for testing
        if gene['chromosome'] > len(summaries):
            continue

        gene_slice = gene_to_slice(gene, promoter_length=length_promoter)
        systematic = gene['systematic']
        orfs[systematic] = summaries[gene['chromosome'] - 1][gene_slice]

    # Iterate through orfs and build mean relative occupancy profile
    profile = np.zeros(length_promoter + length_representative_orf)
    n = 0.
    for orf in orfs.values():
        # Determine amount of profile to update
        length_extract = min(profile.size, orf.size)
        
        # Extract relevant portion of ORF and normalize
        profile_orf = orf['b'][:length_extract]
        profile_orf /= np.sum(profile_orf)
        
        # Update profile
        n += 1.
        profile *= 1. - 1. / n
        profile[:length_extract] += profile_orf / n
    
    # Output profile
    np.savetxt(path_output, profile)

if __name__ == '__main__':
    profile = main()
