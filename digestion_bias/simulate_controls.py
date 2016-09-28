"""Uses dinucleotide statistics from get_dinucleotide_statistics.py to simulate
control reads such that:
1. The read length distributions matches the observed distribution.
2. The joint distribution of cut frequencies from both ends of the reads matches
   the observed distribution.
3. The number of read centers per segment (as estimated from the observed data)
   is conserved.
4. The fragment centers are distributed as uniformly as possible within each
   segment, subject to the above constraints.

Testing code:

import argparse
import itertools
import os
import multiprocessing as mp
import numpy as np
import pandas as pd
import simulate_controls

args = {'reference_sequence': '../Scerevisiae_chr.fna',
        'read_counts': '../XuData/titration/y_Sample_NucDG_8_ACTTGA.txt',
        'regions': '../XuData/titration/regions_Sample_NucDG_8_ACTTGA.txt',
        'lengths': '../XuData/titration/lengthDist_Sample_NucDG_8_ACTTGA.txt',
        'dinucleotide_statistics':
        './joint_dinucleotide_statistics_Sample_NucDG_8_ACTTGA.csv',}

ref = simulate_controls.load_reference_sequence(args['reference_sequence'])
reads = simulate_controls.load_read_counts(args['read_counts'])
regions = simulate_controls.load_regions(args['regions'])
stats = simulate_controls.load_dinucleotide_statistics(
    args['dinucleotide_statistics'])
p_length = simulate_controls.load_length_distribution(args['lengths'])

num_reads = simulate_controls.get_num_reads_by_region(reads, regions)
dn_dist = simulate_controls.DinucleotideDistributions(stats)

counts, cuts, lengths = simulate_controls.draw_chromosome(
    num_reads[1], regions[1], ref['Scchr01'].seq, p_length, dn_dist,
    num_processes=8)
"""

#!/usr/bin/python

import argparse
import itertools
import os
import multiprocessing as mp
import numpy as np
import pandas as pd
from Bio import Seq
from Bio import SeqIO

DINUCLEOTIDES = ['AA', 'AT', 'AG', 'AC',
                 'TA', 'TT', 'TG', 'TC',
                 'GA', 'GT', 'GG', 'GC',
                 'CA', 'CT', 'CG', 'CC']
CUT_INDEX = pd.MultiIndex.from_product([DINUCLEOTIDES, DINUCLEOTIDES],
                                       names=['fwd', 'rev'])
DELIMITER = ','


class DiscreteRV(object):

    """DiscreteDistribution encapsulates a single random variable with arbitrary
    support."""

    def __init__(self, p, values=None):
        """
        Construct a DiscreteRV object.

        Args:
            p: Vector of probabilities per item. This does not need to be
                normalized.
            x: Optional vector of values that the random variable can take.
                Defaults to range(0, len(p)).
        """

        if values is None:
            values = range(0, len(p))
        self.__values = pd.Series(values)
        self.__p = pd.Series(p / np.sum(p))
        self.__p.index = pd.Index(self.__values)
        self.__cdf = self.__p.cumsum()

    def sample(self, include_p=False):
        """Sample from the distribution of a discrete random variable with
        arbitrary support.

        Args:
            include_p: Include probability of drawing resulting sample along
                with the draw?

        Returns:
            x: A draw from the specified discrete distribution
            p: The probablity of drawing this value if include_p is True
        """
        i = min(len(self.__cdf), np.searchsorted(
            self.__cdf.values, np.random.uniform(), 'right'))
        if include_p:
            return self.__values.iloc[i], self.__p.iloc[i]
        return self.__values.iloc[i]

    def pmf(self, val=None):
        """Compute the probability mass function of the random variable at the
        specified value. Returns 0 if the value is not in the variable's
        support.

        Args:
            val: Value at which to compute the PMF. If None, a copy of self.__p
                is returned.

        Returns:
            p: PMF evaluated at x.
        """
        if val is None:
            return self.__p.copy()
        if not val in self.__values.values:
            return 0.
        return self.__p.loc[val]

    def support(self):
        """Return the set of values the random variable can take on."""
        return self.__values


class DinucleotideDistributions(object):

    """Encapsulates marginal and conditional distributions associated with
    cut dinucleotides from paired end data."""

    def __init__(self, stats):
        """Builds marginal and conditional dinucleotide distributions needed for
        rejection sampling from paired-end dinucleotide statistics.

        Args:
            stats: A pandas DataFrame as produced by
                load_dinucleotide_statistics

        Returns:
            fwd: DiscreteRV for the marginal distribution of the forward read
                dinucleotide.
            rev_given_fwd: Dictionary of DiscreteRV's for the conditional
                distribution of reverse read cut dinucleotide given their mate's
                cut dinucleotides.
            rev_given_fwd_control: Dictionary of DiscreteRV's for the
                conditional distribution of reverse read cut dinucleotide given
                their mate's cut dinucleotides from local controls.
        """
        p_fwd = stats.groupby(level=0)['p_obs'].sum()
        self.p_fwd = DiscreteRV(p_fwd.values, p_fwd.index.values)
        self.p_rev_given_fwd = self.conditional_from_dataframe(stats, 'p_obs')
        self.p_rev_given_fwd_control = self.conditional_from_dataframe(
            stats, 'p_control')

    @staticmethod
    def conditional_from_dataframe(dataframe, column='p'):
        """Build a dictionary of DiscreteRV from a DataFrame with a
        multiindex."""
        conditional = {}
        for given in dataframe.index.levels[0]:
            conditional[given] = DiscreteRV(
                p=dataframe.loc[given][column],
                values=dataframe.loc[given].index.values)
        return conditional

    def p_accept(self, fwd, rev):
        """Compute acceptance probability for rejection sampler using
        p_rev_given_fwd_control as proposal and p_rev_given_fwd as target.

        Args:
            fwd: Forward cut dinucleotide to condition on.
            rev: Proposed reverse cut dinucleotide.

        Returns:
            p_accept: Float acceptance probability.
        """
        proposal = self.p_rev_given_fwd_control[fwd]
        target = self.p_rev_given_fwd[fwd]
        env = (target.pmf() / proposal.pmf()).max()
        p_accept = target.pmf(rev) / (env * proposal.pmf(rev))
        if not np.alltrue(np.isfinite(p_accept)):
            return 0.0
        return p_accept


def load_length_distribution(path, sep=' '):
    """Load length distribution into a DiscreteRV.

    Args:
        path: Path to a delimited file of length frequencies

    Returns:
        length_rv: DiscreteRV with the distribution given by the provided length
            frequencies.
    """
    length_dist = pd.read_table(path, sep=sep, names=('length', 'frequency'))
    return DiscreteRV(p=length_dist['frequency'], values=length_dist['length'])


def load_read_counts(path, sep=','):
    """Load read counts.

    Args:
        path: String path to ragged sep-delimited file of regions
        sep: Separator used in delimited file

    Returns:
        reads: A dictionary keyed by chromosome containing a pandas Series of
            read counts per basepair.
    """
    reads = {}
    with open(path, 'rb') as inputs:
        chrom = 0
        for line in inputs:
            chrom += 1
            counts = np.fromstring(line, sep=sep)
            reads[chrom] = pd.Series(counts)
    return reads


def load_dinucleotide_statistics(path):
    """Load dinucleotide statistics from csv created by
    get_dinucleotide_statistics.

    Args:
        path: String path to csv

    Returns:
        stats: A pandas DataFrame of statistics indexed by dinucleotide1 and
        dinucleotide2 with observed and control frequencies and probabilities.
    """
    stats = pd.read_csv(path, sep=DELIMITER, index_col=(0, 1))
    stats['p_obs'] = stats['cut_frequency'] / stats['cut_frequency'].sum()
    stats['p_control'] = stats['control_frequency'] / \
        stats['control_frequency'].sum()
    return stats


def load_regions(path, sep=' '):
    """Load region definitions as produced by raw counts.

    Args:
        path: String path to ragged spaced-delimited file of regions

    Returns:
        regions: A dictionary keyed by chromosome containing a list of region
            slices per chromosome.
    """
    regions = {}
    with open(path, 'rb') as inputs:
        chrom = 0
        for line in inputs:
            chrom += 1
            region_ids = np.fromstring(line, dtype='int', sep=sep)
            uniq = np.unique(region_ids)
            regions[chrom] = [slice(np.min(np.where(region_ids == i)[0]),
                                    np.max(np.where(region_ids == i)[0] + 1))
                              for i in uniq]
    return regions


def load_reference_sequence(path, fmt='fasta'):
    """Load reference sequence into dictionary-like object

    Args:
        path: String path to a reference sequence file
        fmt: The format of the reference sequence file

    Returns:
        A dictionary-like object containing the reference sequence
    """
    return SeqIO.to_dict(SeqIO.parse(path, fmt))


def get_num_reads_by_region(reads, regions):
    """Compute number of reads per region.

    Args:
        reads: Dictionary of vectors of read counts as produced by
            load_read_counts.
        regions: Dictionary of region slices as produced by load_regions.

    Returns:
        num_reads: Dictionary of pandas Series of read counts per region.
    """
    num_reads = {}
    for chrom, regions in regions.iteritems():
        num_reads[chrom] = pd.Series([reads[chrom][r].sum() for r in regions])
    return num_reads


def get_dinucleotides(ref, region, max_length=300):
    """Gets cut dinucleotide for reads with possible centers within the given
    region.

    Args:
        ref: Reference sequence as BioPython Sequence object
        region: Slice defining region
        max_length: Maximum read length

    Returns:
        dinucleotides: pandas DataFrame of cut dinucleotides for forward and
            reverse reads starting at each position that could lead to a center
            within the region. Has columns 'fwd' and 'rev'.
    """
    max_half = int(max_length) / 2 + 1
    first = max(0, region.start - max_half)
    last = min(len(ref), region.stop + max_half)
    dinucleotides = pd.DataFrame({'basepair': np.arange(first, last, dtype=int),
                                  'fwd': np.zeros(last - first, dtype='S2'),
                                  'rev': np.zeros(last - first, dtype='S2')})
    for i in xrange(last - first):
        basepair = i + first
        if basepair == 0:
            dinucleotide = ''
        else:
            dinucleotide = ref[(basepair - 1):(basepair + 1)]
        dinucleotides['fwd'][i] = str(dinucleotide)
        if i > 0:
            dinucleotides['rev'][i - 1] = str(
                Seq.reverse_complement(dinucleotide))
    dinucleotides = dinucleotides.set_index(['fwd', 'rev'])
    dinucleotides.sortlevel(inplace=True)
    return dinucleotides


def pmap_unordered(func, args, num_processes=None):
    """Map a function to a sequence of arguments using multiple processes. This
    exists because multiprocessing.Pool has poor support for closures, which are
    the cleanest way to implement parallel processing in this case.

    Args:
        func: Function to map to arguments.
        args: Iterable sequence of arguments for f. Must also have len().
        num_processes: Number of processes to use.

    Returns:
        values: Unordered list of results from evaluating f on each entry of
            args.
    """
    # Define a mapper function. This wraps f and provides communication with the
    # queue.
    def mapper(args, queue):
        """Mapper function for subsets of supplied arguments."""
        for arg in args:
            queue.put(func(arg))
    # Decide on number of processes and number of chunks.
    if num_processes is None:
        num_processes = mp.cpu_count()
        print 'Using %d processes by CPU count' % num_processes
    chunksize = len(args) / num_processes
    if len(args) % num_processes != 0:
        chunksize += 1
    procs = []
    queue = mp.JoinableQueue()
    num_map_shards = 0
    for i in xrange(num_processes):
        chunk = args[(i * chunksize):min((i + 1) * chunksize, len(args))]
        if len(chunk) < 1:
            continue
        num_map_shards += len(chunk)
        proc = mp.Process(target=mapper, args=(chunk, queue))
        proc.start()
        procs.append(proc)
    assert num_map_shards == len(args)
    values = []
    while len(values) < len(args):
        val = queue.get()
        values.append(val)
    return values


def draw_region(num_reads, region, ref, p_length, dn_dist):
    """Draw a set of reads from a region assuming uniform positioning
    conditioning on MNase dinucleotide preference and the observed distribution
    of fragment lengths.

    Args:
        num_reads: Integer number of reads to draw within the given region.
        ref: Reference sequence for the entire chromosome that region is on as a
            BioPython Seq object. NB: This is NOT a SeqRecord object.
        region: Slice defining the region of interest, as produced by
            load_regions.
        p_length: DiscreteRV for the marginal distribution of fragment lengths
            as produced by load_length_distribution.
        dn_dist: DinucleotideDistributions object containing all marginal and
            conditional cut dinucleotide distributions needed for sampling.

    Returns:
        reads: Integer vector of centers of the sampled reads. This is its
            center along its chromosome starting from 0, NOT within its region.
            This is randomly rounded if the center is fractional.
        cuts: DataFrame of cut statistics, indexed by (fwd, rev) cut
            dinucleotide with counts in a column named num_reads.
        lengths: Series of read length statistics indexed by length with values
            equal to the number of reads with each length.
    """
    # Build cut site to position index.
    dinucleotides = get_dinucleotides(ref, region,
                                      max_length=max(p_length.support()))
    # Build data structure for cut counts.
    cuts = pd.DataFrame({'num_reads': np.zeros(len(CUT_INDEX)),
                         'p': np.zeros(len(CUT_INDEX))},
                        index=CUT_INDEX)
    # Enumerate feasible draws for each dinucleotide pair and get probabilities
    # for multinomial draw of cut pair counts.
    feasible = {}
    min_length = p_length.support().min()
    max_length = p_length.support().max()
    for fwd, rev in cuts.index.to_series():
        possible = pd.DataFrame.from_records(
            itertools.ifilter(
                lambda pair: ((pair[1] - pair[0]) >= min_length and
                              (pair[1] - pair[0]) <= max_length),
                itertools.product(
                    dinucleotides['basepair'].loc[fwd],
                    dinucleotides['basepair'].loc[(slice(None), rev)])),
            columns=['fwd', 'rev'])
        possible['length'] = possible['rev'] - possible['fwd']
        possible['center'] = (possible['fwd'] + possible['rev']) / 2.
        possible = possible[possible['center'] > region.start]
        possible = possible[possible['center'] < region.stop]
        # Only include the cut dinucleotide pair and its feasible set if it has
        # members.
        if len(possible) > 0:
            possible['p'] = [p_length.pmf(l) for l in possible['length']]
            possible['p'] /= possible['p'].sum()
            dn = (fwd, rev)
            feasible[dn] = possible.copy()
            cuts['p'].loc[dn] = dn_dist.p_fwd.pmf(fwd) * \
                    dn_dist.p_rev_given_fwd[fwd].pmf(rev)
    # Normalize cut pair probabilities to account for omitted infeasible pairs
    # and execute the multinomial draw.
    cuts['p'] /= cuts['p'].sum()
    cuts['num_reads'] = np.random.multinomial(num_reads, cuts['p'])
    # Setup data structure for drawn fragment length histogram.
    lengths = pd.Series(np.zeros(len(p_length.support())),
                        index=p_length.support())
    # We iterate over the feasible cut dinucleotide pairs, which is a set of
    # constant size, instead of the individual reads we want to draw.
    reads = {}
    for cut_dn in feasible.keys():
        num_to_draw = cuts['num_reads'].loc[cut_dn]
        # Skip cut dinucleotide pairs for which we sampled 0 reads. This checks
        # the result of the multinomial draw, which is distinct from the
        # previous check for the existence of feasible fragments.
        if num_to_draw < 1:
            continue
        possible = feasible[cut_dn]
        # Sample the given number of reads from the set of feasible fragments,
        # allowing for duplication and weighting by the length distribution.
        # Note that this allows the realized length distribution to deviate from
        # the target one but guarantees an exact match to the target cut
        # dinucleotide pair distribution, as intended.
        read_ids = np.random.choice(possible.index, size=num_to_draw,
                                    replace=True, p=possible['p'])
        draws = possible.loc[read_ids]
        # Randomly round read centers that end up at fractional positions.
        rand_round = np.random.randint(1, size=len(draws))
        centers = (np.floor(draws['center']) * (1 - rand_round) +
                   np.ceil(draws['center']) * rand_round)
        reads[cut_dn] = centers
        # We update the read length histogram in one shot, relying on
        # index-based series alignment.
        lengths_hist = np.bincount(draws['length'],
                                   minlength=max(p_length.support())+1)
        lengths_hist = pd.Series(lengths_hist[min(p_length.support()):],
                                 index=p_length.support())
        lengths += lengths_hist
    if len(reads) > 0:
        reads = np.concatenate(reads.values()).astype(np.int64)
    else:
        reads = np.empty(0, dtype=np.int64)
    return reads, cuts, lengths


def draw_chromosome(num_reads, regions, ref, p_length, dn_dist,
                    num_processes=None):
    """Draw a set of reads from a chromosome assuming uniform positioning
    conditioning on MNase dinucleotide preference and the observed distribution
    of fragment lengths.

    Args:
        num_reads: Dictionary keyed by region ID containing the number of reads
            to sampler per region.
        ref: Reference sequence for the entire chromosome that region is on as a
            BioPython Seq object. NB: This is NOT a SeqRecord object.
        regions: Dictionary keyed by region ID of slices defining the region of
            interest, as produced by load_regions.
        p_length: DiscreteRV for the marginal distribution of fragment lengths
            as produced by load_length_distribution.
        dn_dist: DinucleotideDistributions object containing all marginal and
            conditional cut dinucleotide distributions needed for sampling.
        num_processes: Optional number of processes to use for parallel sampling
            at the region level.

    Returns:
        reads: Integer vector of the number of read centers sampled at each base
            pair along the chromosome.
        cuts: DataFrame of cut statistics, indexed by (fwd, rev) cut
            dinucleotide with counts in a column named num_reads.
        lengths: Series of read length statistics indexed by length with values
            equal to the number of reads with each length.
    """
    counts = np.zeros(len(ref), dtype=int)
    cuts = pd.DataFrame({'num_reads': np.zeros(len(CUT_INDEX))},
                        index=CUT_INDEX)
    lengths = pd.Series(np.zeros(len(p_length.support())),
                        index=p_length.support())
    draws = None
    # Branch for serial logic
    if num_processes is not None and num_processes <= 1:
        draws = [None] * len(regions)
        for i, region in enumerate(regions):
            draws[i] = draw_region(num_reads[i], region, ref, p_length,
                                   dn_dist)
    else:
        run_region = lambda i: draw_region(num_reads[i], regions[i], ref,
                                           p_length, dn_dist)
        draws = pmap_unordered(run_region, range(len(regions)),
                               num_processes=num_processes)
    # Compute aggregates from individual draws
    for centers, cutstats, lengthstats in draws:
        counts += np.bincount(centers, minlength=len(counts))
        cuts += cutstats
        lengths += lengthstats
    return counts, cuts, lengths


def draw_genome(num_reads, regions, ref, p_length, dn_dist, num_processes=None,
                ref_pattern='Scchr{:02d}'):
    """Draw a set of reads from genome assuming uniform positioning conditioning
    on MNase dinucleotide preference and the observed distribution of fragment
    lengths.

    Args:
        num_reads: Dictionary keyed by chromosome of dictionaries keyed by
            region ID containing the number of reads to sampler per region.
        ref: Reference sequence as a dictionary of SeqRecord objects keyed by
            chromsome name.
        regions: Dictionary keyed by chromosome of dictionaries keyed by region
            ID of slices defining the region of interest, as produced by
            load_regions.
        p_length: DiscreteRV for the marginal distribution of fragment lengths
            as produced by load_length_distribution.
        dn_dist: DinucleotideDistributions object containing all marginal and
            conditional cut dinucleotide distributions needed for sampling.
        num_processes: Optional number of processes to use for parallel sampling
            at the region level.
        ref_pattern: string.format pattern to convert chromosome numbers
            (starting from 1) to names in ref.

    Returns:
        reads: Dictionary keyed by chromosome of integer vectors of the number
            of read centers sampled at each base pair along the chromosome.
        cuts: DataFrame of cut statistics, indexed by (fwd, rev) cut
            dinucleotide with counts in a column named num_reads.
        lengths: Series of read length statistics indexed by length with values
            equal to the number of reads with each length.
    """
    reads = {}
    cuts = pd.DataFrame({'num_reads': np.zeros(len(CUT_INDEX))},
                        index=CUT_INDEX)
    lengths = pd.Series(np.zeros(len(p_length.support())),
                        index=p_length.support())
    for chrom, chrom_regions in regions.iteritems():
        print "Starting chromosome ", chrom
        ref_name = ref_pattern.format(chrom)
        counts, cutstats, lengthstats = draw_chromosome(
            num_reads[chrom], chrom_regions, ref[ref_name].seq, p_length,
            dn_dist, num_processes)
        reads[chrom] = counts
        cuts += cutstats
        lengths += lengthstats
        print "Finished chromosome ", chrom
    return reads, cuts, lengths


def save_counts(counts, path):
    """Saves counts from dictionary with keys containing ids to a ragged CSV
    file."""
    # Make directory if needed
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    with open(path, 'wb') as out:
        ids = sorted(counts.keys())
        for i in ids:
            np.savetxt(
                out, counts[i][np.newaxis, :], fmt='%.1f', delimiter=',')


def save_stats(stats, path):
    """Saves diagnostic statistics to CSVs."""
    # Make directory if needed
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    stats.to_csv(path)


def main():
    """Main function for script.
    Handles argument parsing and calls worker functions.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference_sequence', default='',
                        help='Path to reference sequence file (FASTA)',
                        required=True)
    parser.add_argument('--regions', default='',
                        help='Path to regions as computed by '
                        'cplate_segment_genome',
                        required=True)
    parser.add_argument('--read_counts', default='',
                        help='Path to ragged CSV of read center counts',
                        required=True)
    parser.add_argument('--lengths', default='',
                        help='Path to file of length frequencies',
                        required=True)
    parser.add_argument('--dinucleotide_statistics', default='',
                        help='Path to dinucleotide statistics CSV',
                        required=True)
    parser.add_argument('--output_dir', default='',
                        help='Destination directory for output',
                        required=True)
    parser.add_argument('-p', '--processes', default=None, type=int,
                        help='Number of processes to use')
    parser.add_argument('--rng_seed', default=20131109,
                        help='Seed for numpy RNG')
    args = parser.parse_args()

    # Load reference sequence, read summaries, and dinucleotide statistics
    reference_sequence = load_reference_sequence(args.reference_sequence)
    reads = load_read_counts(args.read_counts)
    regions = load_regions(args.regions)
    p_length = load_length_distribution(args.lengths)
    stats = load_dinucleotide_statistics(args.dinucleotide_statistics)

    # Compute region-level statistics and dinucleotide distributions
    num_reads = get_num_reads_by_region(reads, regions)
    dn_dist = DinucleotideDistributions(stats)

    # Run sampling
    np.random.seed(args.rng_seed)
    counts, cuts, lengths = draw_genome(num_reads, regions, reference_sequence,
                                        p_length, dn_dist,
                                        num_processes=args.processes)

    # Save count output
    save_counts(counts, os.path.join(args.output_dir, 'control.csv'))
    # Save diagnostic statistics
    save_stats(cuts, os.path.join(args.output_dir, 'cuts.csv'))
    save_stats(lengths, os.path.join(args.output_dir, 'lengths.csv'))

if __name__ == '__main__':
    main()
