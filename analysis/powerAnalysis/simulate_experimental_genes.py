#!/usr/bin/python

import gc
import itertools
import re
import StringIO

import numpy as np
import pandas as pd
from scipy import stats

# Constants

# Number of chromosomes
N_CHROM = 16

# RNG seed
np.random.seed(20131109)

# Data paths

# Index of ORFs
path_gene_index = '../geneInfo/geneIndex.txt'

# Path for representative gene information
path_rep_gene = 'powerAnalysis/data/representative_gene.txt'

# Path for template
path_template = 'templates/template_H_1-combined.txt'

# Path for reads
path_reads = '../XuData/H_1-combined_chrom.txt'

# Path for regions
path_regions = '../XuData/regionDefs/regions_H_1-combined_800.txt'

# Pattern for parameter information
pattern_param = 'results/mcmc_params_H_1-combined_chrom%02d.txt'

# Pattern for output
pattern_output = 'powerAnalysis/data/simChrom_%s.%s'


# Representative gene design parameters

# Length to extract as promoter
length_promoter = 1000

# Length of ORF in representative gene
length_representative_orf = 2500

length_rep_gene = length_promoter + length_representative_orf + 1

# Truncation quantile for log-normal background variation
bkg_quantile = 1. - 5. / 147.


# Experimental design

# Number of replicates
n_replicates = 10

# Quantiles for coverage levels
design_coverage_quantiles = np.arange(0.05, 1., 0.1)

# Spacings for alternative positioning
design_alt_pos_spacing = np.array((0, 5, 10, 15, 20, 25, 30, 35, 40, 45),
                                  dtype=int)

# Relative magnitudes of alternative positions (alternative / primary)
design_alt_pos_magnitude = np.arange(0., 1.1, 0.1)

design_factors = {'coverage_quantile' : design_coverage_quantiles,
                  'alt_pos_spacing' : design_alt_pos_spacing,
                  'alt_pos_magnitude' : design_alt_pos_magnitude}


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
    # Initialize design matrix
    n_levels = [len(factor) for factor in design_factors.values()]
    n_settings = np.prod(n_levels)

    design_data = [(key, np.zeros(n_settings, dtype=factor.dtype)) for key,
                   factor in design_factors.iteritems()]
    design_data = dict(design_data)

    design = pd.DataFrame(data=design_data)

    # Iterate through possible combinations of factors to fill design matrix
    for i, setting in enumerate(itertools.product(*[design_factors[s] for s in
                                                    design.columns])):
        design.ix[i] = setting

    # Load representative gene
    rep_gene = pd.read_csv(path_rep_gene, index_col=0)

    # Load template
    template = np.loadtxt(path_template)
    
    # Load reads data
    reads = []
    with open(path_reads, 'rb') as f:
            chrom = 0
            for line in f:
                    chrom += 1
                    reads.append(pd.DataFrame(
                        {'chrom' : chrom,
                         'y' : np.fromstring(line.strip(), sep=',')}))
    
    # Load regions data
    regions = []
    with open(path_regions, 'rb') as f:
            chrom = 0
            for line in f:
                    chrom += 1
                    regions.append(pd.DataFrame(
                        {'chrom' : chrom,
                         'region_id' : np.fromstring(line.strip(), sep=' ',
                                                     dtype=int)}))
    
    # Add region information to reads
    for reads_df, regions_df in itertools.izip(reads, regions):
        reads_df['region_id'] = regions_df['region_id']
    
    # Combine regions into single DataFrame
    regions = pd.concat(regions)

    # Load parameter posterior summaries
    params = []
    for chrom in xrange(1, N_CHROM + 1):
        with open(pattern_param % chrom, 'rb') as f:
            params.append(pd.read_table(f, sep=' '))
        params[chrom-1]['chrom'] = chrom
        params[chrom-1].index = pd.MultiIndex.from_arrays(
            [params[chrom-1]['chrom'], params[chrom-1]['region_id']])
    
    params = pd.concat(params)
    
    # Load and parse gene index
    with open(path_gene_index, 'rb') as f:
        gene_index = parse_gene_index(f)

    # Label each ORF and its promoter
    orfs = []
    for gene in gene_index:
        # Skip missing chromosomes; primarily for testing
        if gene['chromosome'] > len(reads):
            continue

        gene_slice = gene_to_slice(gene, promoter_length=length_promoter)
        systematic = gene['systematic']
        chrom = gene['chromosome']

        # Setup ORF entry
        orf = reads[chrom-1][gene_slice]
        orf['systematic'] = systematic
        orfs.append(orf)
    
    # Combine ORFs into a single dataframe
    orfs = pd.concat(orfs)
    
    # Combine reads into single DataFrame
    reads = pd.concat(reads)

    # Compute ORF properties
    coverage_by_orf = orfs.groupby('systematic')['y'].mean()
    coverage_quantiles = [coverage_by_orf.quantile(p) for p in
                          design_factors['coverage_quantile']]
    coverage_quantiles = dict(zip(design_factors['coverage_quantile'],
                                  coverage_quantiles))

    # Compute parameter summaries by coverage quantile
    n_quantiles = np.size(design_factors['coverage_quantile'])
    coverage_by_region = reads.groupby(['chrom', 'region_id']).mean()
    coverage_by_region_quantiles = [coverage_by_region.quantile(p) for p in
                                    np.linspace(0., 1., n_quantiles + 1)]
    params['coverage'] = coverage_by_region['y']
    
    param_medians = [params[(params['coverage'] >
                             coverage_by_region_quantiles[i]) &
                             (params['coverage'] <=
                              coverage_by_region_quantiles[i+1])].median()
                     for i in xrange(n_quantiles)]
    mu = dict((q, medians['mu_postmean']) for (q, medians) in
              itertools.izip(design_factors['coverage_quantile'],
                             param_medians))
    sigma = dict((q, medians['sigma_postmean']) for (q, medians) in
                 itertools.izip(design_factors['coverage_quantile'],
                                param_medians))

    # Compute cluster center indices from TSS-relative positions
    primary_positions = (length_promoter +
                         np.round(rep_gene['center']).astype(int))

    # Cleanup before large memory allocations
    del reads, regions, params, orfs
    gc.collect()
    
    # Initialize chromosome vectors
    theta = np.zeros((n_replicates, design.shape[0] * length_rep_gene))
    b = np.zeros_like(theta)
    region_types = np.zeros(theta.shape, dtype=int)

    # Setup DataFrame for ground truth
    n_primary_positions = (design.shape[0] * np.size(primary_positions) *
                           n_replicates)
    n_alt_positions = (2 * np.size(primary_positions) *
                       np.sum((design['alt_pos_magnitude'] > 0) *
                              (design['alt_pos_spacing'] > 0)) *
                       n_replicates)
    n_positions = n_primary_positions + n_alt_positions
    groundtruth = pd.DataFrame({'rep' : np.zeros(n_positions, dtype=int),
                                'pos' : np.zeros(n_positions, dtype=int),
                                'b' : np.zeros(n_positions, dtype=float),
                                'type' : np.zeros(n_positions, dtype=int),
                                'offset' : np.zeros(n_positions, dtype=int),
                                'magnitude' : np.zeros(n_positions,
                                                       dtype=float),
                                'coverage' : np.zeros(n_positions, dtype=float),
                                'eff.magnitude' : np.zeros(n_positions,
                                                           dtype=float)
                               })
    
    row_groundtruth = 0
    for rep in xrange(n_replicates):
        for region_id in design.index:
            # Extract design information
            coverage_quantile = design.get_value(region_id, 'coverage_quantile')
            alt_pos_spacing = design.get_value(region_id, 'alt_pos_spacing')
            alt_pos_magnitude = design.get_value(region_id, 'alt_pos_magnitude')

            # Compute limits of gene
            gene_slice = slice(region_id * length_rep_gene,
                               (region_id + 1) * length_rep_gene)

            # Pull out relevant region of theta as a view
            theta_gene = theta[rep, gene_slice]

            # Set region type
            region_types[rep, gene_slice] = region_id

            # Compute total expected reads for gene
            n_reads = (coverage_quantiles[coverage_quantile] *
                       length_rep_gene)

            # Draw background noise from truncated log-normal
            log_bkg_max = (mu[coverage_quantile] + sigma[coverage_quantile] *
                           stats.norm.ppf(bkg_quantile))
            theta[rep, gene_slice] = (mu[coverage_quantile] +
                                      sigma[coverage_quantile] *
                                      np.random.randn(length_rep_gene))
            while theta[rep, gene_slice].max() > log_bkg_max:
                truncated = np.where(theta[rep, gene_slice] > log_bkg_max)[0]
                theta_gene[truncated] = (mu[coverage_quantile] +
                                         sigma[coverage_quantile] *
                                         np.random.randn(truncated.size))

            # Add primary and alternative positions
            if alt_pos_spacing > 0 and alt_pos_magnitude > 0:
                # Have alternative positions
                # Locate alternative positions
                alt_positions = np.r_[primary_positions - alt_pos_spacing,
                                      primary_positions + alt_pos_spacing]

                # Compute occupancy for primary and alternative positions
                position_theta = (
                    np.log(rep_gene['relative_occupancy']) +
                    np.log(n_reads - np.exp(theta_gene).sum() +
                           np.exp(theta_gene[primary_positions]).sum() +
                           np.exp(theta_gene[alt_positions]).sum()))

                # Allocate occupancy to primary and alternative positions
                primary_theta = (position_theta -
                                 np.log(1. + 2. * alt_pos_magnitude))
                alt_theta = primary_theta + np.log(alt_pos_magnitude)
                alt_theta = np.r_[alt_theta, alt_theta]

                # Assign occupancies to primary and alternative positions
                theta_gene[primary_positions] = primary_theta
                theta_gene[alt_positions] = alt_theta
            else:
                # Only have primary positions
                # Compute occupancy for primary positions
                primary_theta = (
                    np.log(rep_gene['relative_occupancy']) +
                    np.log(n_reads - np.exp(theta_gene).sum() +
                           np.exp(theta_gene[primary_positions]).sum()))

                # Assign occupancies to primary positions
                theta_gene[primary_positions] = primary_theta
                
            # Add positions to groundtruth

            # Always add primary positions to groundtruth
            primary_effective_magnitude = 1. / (1 + 2. * alt_pos_magnitude *
                                                (alt_pos_spacing > 0))
            groundtruth_slice = slice(row_groundtruth, row_groundtruth +
                                      np.size(primary_positions))
            groundtruth['rep'][groundtruth_slice] = rep
            groundtruth['pos'][groundtruth_slice] = (
                primary_positions.values + gene_slice.start)
            groundtruth['b'][groundtruth_slice] = np.exp(primary_theta)
            groundtruth['type'][groundtruth_slice] = 0 # Primary = 0
            groundtruth['magnitude'][groundtruth_slice] = 1.
            groundtruth['offset'][groundtruth_slice] = 0
            groundtruth['coverage'][groundtruth_slice] = (
                coverage_quantiles[coverage_quantile])
            groundtruth['eff.magnitude'][groundtruth_slice] = (
                primary_effective_magnitude)
            
            row_groundtruth += np.size(primary_positions)

            # Add alternative positions if needed
            if alt_pos_spacing > 0 and alt_pos_magnitude > 0:
                groundtruth_slice = slice(row_groundtruth, row_groundtruth +
                                          np.size(alt_positions))
                groundtruth['rep'][groundtruth_slice] = rep
                groundtruth['pos'][groundtruth_slice] = (
                    alt_positions + gene_slice.start)
                groundtruth['b'][groundtruth_slice] = np.exp(alt_theta)
                groundtruth['type'][groundtruth_slice] = 1 # Alternative = 0
                groundtruth['magnitude'][groundtruth_slice] = alt_pos_magnitude
                groundtruth['offset'][groundtruth_slice] = alt_pos_spacing
                groundtruth['coverage'][groundtruth_slice] = (
                    coverage_quantiles[coverage_quantile])
                groundtruth['eff.magnitude'][groundtruth_slice] = (
                    alt_pos_magnitude / (1 + 2. * alt_pos_magnitude))

                row_groundtruth += np.size(alt_positions)
        
    # Exponentiate to determine coefficients
    b = np.exp(theta)
    
    # Convolve coefficients with template to determine Poisson rates
    lmbda = np.array([np.convolve(b_rep, template, mode='same') for b_rep in
                      b])

    # Simulate read counts from Poisson distribution
    y = np.random.poisson(lmbda)

    # Generate null reads
    null = np.empty_like(y)
    p_null = np.ones(length_rep_gene) / (length_rep_gene + 0.)

    for rep in xrange(n_replicates):
        for region_id in design.index:
            # Compute limits of gene
            gene_slice = slice(region_id * length_rep_gene,
                               (region_id + 1) * length_rep_gene)

            # Permute reads within gene for null. This amounts to a multinomial
            # draw, not permutation of the counts in y; working at the
            # individual read level, not the read count level.
            null[rep, gene_slice]  = np.random.multinomial(
                y[rep, gene_slice].sum(), p_null)

    # Save output

    # Log-coefficients (theta)
    with open(pattern_output % ('theta', 'npy'), 'wb') as f:
        np.save(f, theta)

    # Coefficients (beta)
    with open(pattern_output % ('b', 'npy'), 'wb') as f:
        np.save(f, b)
    
    # Expected counts (lambda)
    with open(pattern_output % ('lambda', 'npy'), 'wb') as f:
        np.save(f, lmbda)
    
    # Region information
    with open(pattern_output % ('regions', 'txt'), 'wb') as f:
        np.savetxt(f, region_types, delimiter=' ', fmt='%d')
    
    # Simulated read counts
    with open(pattern_output % ('y', 'txt'), 'wb') as f:
        np.savetxt(f, y, delimiter=',', fmt='%d')
    
    # Null read counts
    with open(pattern_output % ('null', 'txt'), 'wb') as f:
        np.savetxt(f, null, delimiter=',', fmt='%d')

    # Design matrix
    with open(pattern_output % ('design', 'txt'), 'wb') as f:
        design.to_csv(f, index=False, float_format='%0.3f')
    
    # Groundtruth positions
    with open(pattern_output % ('groundtruth', 'txt'), 'wb') as f:
        groundtruth.to_csv(f, index=False)

    return design, region_types, theta, b, lmbda, y, groundtruth


if __name__ == '__main__':
    design, region_types, theta, b, lmbda, y, groundtruth = main()

