#!/usr/bin/python

import argparse
import sys

import numpy as np
import pandas as pd
from scipy import stats

# Constants

RNG_SEED = 20131109

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",
                        default="powerAnalysis/data/simChrom_y.txt")
    parser.add_argument("--output_prefix",
                        default="comparisons/puffin/data/simulated/")
    parser.add_argument("--digestion_error_dist",
                        default="powerAnalysis/data/digestion_error_dist.txt")
    return parser.parse_args(argv)

def draw_read(center, cdf):
    return center + (np.array([-1, 1]) *
        (np.searchsorted(cdf, np.random.uniform(size=2)) - 1))

def counts_to_bed(input_file, cdf, output_prefix):
    chrom = 0
    for line in open(input_file, "rb"):
        chrom += 1
        try:
            y = np.fromstring(line, dtype=int, sep=",")
        except:
            continue

        bed = []
        for i in xrange(y.size):
            for n in xrange(y[i]):
                read = draw_read(i, cdf)
                if np.any(read < 0):
                    continue
                bed.append("chr%d %d %d\n" % (chrom, read[0], read[1]))

        out = open(output_prefix + "reads_%02d.bed" % chrom, "wb")
        out.writelines(bed)
        out.close()
    return chrom

def load_digestion_error_cdf(fname):
    # Only keeping the CDF. The offsets are not needed because we are working
    # from simulated centers.
    return np.cumsum(np.loadtxt(fname, unpack=True)[1])

def main():
    np.random.seed(RNG_SEED)
    args = parse_args(sys.argv[1:])
    cdf = load_digestion_error_cdf(args.digestion_error_dist)
    counts_to_bed(args.input, cdf, args.output_prefix)

if __name__ == '__main__':
    main()

