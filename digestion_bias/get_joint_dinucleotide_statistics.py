#!/usr/bin/python
"""
Compute marginal cut site dinucleotide statistics from a SAM/BAM file of aligned
reads and a reference genome.
"""

import argparse
import itertools
import sys

import pysam
from Bio import SeqIO

DINUCLEOTIDES = ("AA", "AT", "AG", "AC",
                 "TA", "TT", "TG", "TC",
                 "GA", "GT", "GG", "GC",
                 "CA", "CT", "CG", "CC")
HEADER = ("dinucleotide1", "dinucleotide2", "cut_frequency",
          "control_frequency")
DELIMITER = ","

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--aligned_reads_source",
                        default="",
                        help="Path to aligned reads file (SAM or BAM).")
    parser.add_argument("--reference_sequence_source",
                        default="",
                        help="Path to file containing reference sequence"
                        " (FASTA).")
    parser.add_argument("--destination",
                        default="",
                        help="Destination for frequency results."
                        " Default is stdout.")
    parser.add_argument("--window_size",
                        default=10L,
                        type=int,
                        help="Window size around cut site for control"
                        " dinucleotides.")
    return parser.parse_args(argv)

def load_reference_sequence(path, fmt="fasta"):
    """Load reference sequence into dictionary-like object

    Args:
        path: String path to a reference sequence file
        format: The format of the reference sequence file

    Returns:
        A dictionary-like object containing the reference sequence
    """
    return SeqIO.to_dict(SeqIO.parse(path, fmt))

def get_cut_dinucleotides(aligned_read, reference_sequence, getrname):
    """
    aligned_read must be forward, not reverse complement. A cut site is
    extracted from both the beginning and end of the given read, and both
    dinucleotides are returned as (string, string). Returns None if either cut
    site is not actually a dinucleotide.
    """
    # Extract rname (chromosome) for aligned read
    rname = getrname(aligned_read.tid)
    # Extract the extended sequence for the read, including the full cut site
    first, last = (aligned_read.pos, aligned_read.pos + aligned_read.tlen)
    sequence = reference_sequence[rname].seq[first - 1 : last + 1]
    # Return only the dinucleotide containing the cut site.
    # Sorting the dinucleotides reduces the number of distinct pairs by
    # enforcing symmetry between fragment ends.
    cuts = [str(sequence[:2]), str(sequence[-2:].reverse_complement())]
    # Validity checks
    if not all([c in DINUCLEOTIDES for c in cuts]):
        return None
    return cuts

def get_control_dinucleotides(aligned_read, reference_sequence, getrname,
                              window_size=10L):
    """Gets list of control dinucleotide pairs from windows around an aligned
    read's cut site.
    """
    # Extract rname (chromosome) for aligned read
    rname = getrname(aligned_read.tid)
    half_window = int(window_size / 2)

    # Extract the extended sequence for the read, including the full cut site
    first, last = (min(aligned_read.positions),
                   max(aligned_read.positions) + 1L)
    if aligned_read.is_reverse:
        sequence = reference_sequence[rname].seq[first - half_window :
                                                 last + 1 + half_window]
        sequence = sequence.reverse_complement()
    else:
        sequence = reference_sequence[rname].seq[first - 1 - half_window :
                                                 last + half_window]

    # Return only the dinucleotide containing the cut site
    return [str(sequence[i : i + 2]) for i in range(2 * half_window + 1)]

class FakeMate(object):
    """Incomplete version of pysam.AlignedRead to avoid Samfile.mate()
    lookups."""
    def __init__(self, read):
        """Build fake mate from pysam.AlignedRead."""
        self.tid = read.tid
        self.aend = read.pos + read.tlen
        self.pos = self.aend - read.rlen
        self.positions = range(self.pos, self.aend)
        self.is_reverse = True

def get_dinucleotide_frequencies(aligned_reads_path, reference_sequence,
                                 window_size=10):
    """Gets observed and control cut frequencies for each aligned read in the
    given file using the given reference sequence.
    """
    cut_frequencies = dict(zip([DELIMITER.join(pair) for pair in
                                itertools.product(DINUCLEOTIDES, repeat=2)],
                               [0] * len(DINUCLEOTIDES)**2))
    control_frequencies = dict(zip([DELIMITER.join(pair) for pair in
                                    itertools.product(DINUCLEOTIDES, repeat=2)],
                                   [0] * len(DINUCLEOTIDES)**2))
    # Reading the reads one line at a time; minimal memory usage
    with pysam.Samfile(aligned_reads_path, 'rb') as sam:
        for read in sam:
            # Discarding unknown chromosomes
            if read.tid < 0:
                continue
            # Only using the forward read of each pair
            if not read.is_proper_pair or read.tlen < 0:
                continue
            # Filter as in Gaffney et al 2012
            if read.mapq < 10:
                continue
            # Make a fake mate, relying on duck typing
            mate = FakeMate(read)
            # Get cut dinucleotides from each end of the read
            cuts = get_cut_dinucleotides(read, reference_sequence, sam.getrname)
            if cuts is None:
                continue
            cut_frequencies[DELIMITER.join(cuts)] += 1
            # Count local controls within a small window
            controls = [
                get_control_dinucleotides(read, reference_sequence,
                                          sam.getrname, window_size),
                get_control_dinucleotides(mate, reference_sequence,
                                          sam.getrname, window_size),
            ]
            for pair in itertools.product(*controls):
                if any([x not in DINUCLEOTIDES for x in pair]):
                    continue
                control_frequencies[DELIMITER.join(pair)] += 1

    return cut_frequencies, control_frequencies


def main():
    """Driver function for the script. Handles argument parsing and IO.
    """
    args = parse_args(sys.argv[1:])

    if args.reference_sequence_source == "":
        raise ValueError("Reference sequence source not specified")
    if args.aligned_reads_source == "":
        raise ValueError("Aligned reads source not specified")

    reference_sequence = load_reference_sequence(args.reference_sequence_source)
    cut_frequencies, control_frequencies = get_dinucleotide_frequencies(
        args.aligned_reads_source, reference_sequence)

    # Write output
    if args.destination == "":
        output_file = sys.stdout
    else:
        output_file = open(args.destination, 'wb')

    print >> output_file, DELIMITER.join(HEADER)
    for pair in itertools.product(DINUCLEOTIDES, repeat=2):
        cuts = DELIMITER.join(pair)
        if cut_frequencies[cuts] < 1 and control_frequencies[cuts] < 1:
            continue
        print >> output_file, DELIMITER.join(
            (cuts, str(cut_frequencies[cuts]), str(control_frequencies[cuts])))

    if output_file is not sys.stdout:
        output_file.close()

if __name__ == "__main__":
    main()
