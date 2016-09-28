#!/usr/bin/python

import gc
import getopt
import os
import sys
import tarfile

import pandas as pd
import yaml

SEP = ' '
HELP = '''
Usage: condense_summaries [options] CONFIG ARCHIVE [CONFIG ACHIVE ...]

Options:
  -h, --help            Show this help message and exit
  -c CHROM, --chrom=CHROM
                        Comma-separated indices of chromosomes to analyze;
                        defaults to 1
  --all                 Run all chromosomes

Details of the required format for the YAML CONFIG files can be found it further
documentation.
'''

# Variables to extract from summaries
VARS = ['theta', 'theta_med', 'se_theta', 'b', 'b_med', 'se_b', 'n_eff',
        'p_local_concentration_pm0', 'p_local_concentration_pm1',]
# Detections to include
PM = (0, 1)


def main():
    '''
    Main function for option-parsing and startup.
    '''
    # Set default values for options
    chrom_list  = None
    run_all     = False
    
    # Parse arguments and options
    opts, args = getopt.getopt(sys.argv[1:], "hc:",
                               ["help", "chrom=", "all"])
    for option, value in opts:
        if option in ('-h', "--help"):
            print >> sys.stderr, HELP
            sys.exit(2)
        elif option in ('-c', '--chrom'):
            chrom_list = [int(x) for x in value.split(',')]
        elif option == '--all':
            run_all = True
        else:
            print >> sys.stderr, "Error -- unknown option %s" % option
            sys.exit(1)

    # Check for logical consistency
    if run_all and (chrom_list is not None):
        print >> sys.stderr, "Error -- cannot have all with chrom"
        sys.exit(1)

    # Set default chrom value
    if chrom_list is None:
        chrom_list = [1]

    if len(args) > 0 and len(args) % 2 == 0:
        cfg_paths = args[::2]
        archive_paths = args[1::2]
    else:
        print >> sys.stderr, ("Error -- need paired paths to YAML "
                              "configurations and output archives")
        sys.exit(1)
    
    # Iterate over configurations
    for cfg_path, archive_path in zip(cfg_paths, archive_paths):
        # Parse YAML configuration
        cfg_file = open(cfg_path, 'rb')
        cfg = yaml.load(cfg_file)
        cfg_file.close()
        
        if run_all:
            chrom_list = range(1, cfg['data']['n_chrom']+1)
        
        # Detect compression for archive
        if archive_path[-2:] == 'gz':
            compression = 'gz'
        elif archive_path[-3:] == 'bz2':
            compression = 'bz2'
        else:
            compression = ''
        
        # Setup archive
        archive = tarfile.open(archive_path, 'w:' + compression)

        # Extract patterns for MCMC summaries
        pattern_summaries = cfg['mcmc_output']['summary_pattern'].format(**cfg)
        pattern_detections = cfg['mcmc_output']['detections_pattern']
        pattern_detections = pattern_detections.format(**cfg)

        # Extract path to scratch directory
        path_scratch = cfg['mcmc_summaries']['path_scratch']

        # Iterate over chromosomes
        for chrom in chrom_list:
            # Load summaries from chromosome
            path_summaries = pattern_summaries % chrom
            summaries = pd.read_table(path_summaries, sep=SEP)

            # Subset summaries to desired columns
            summaries = summaries[VARS]

            # Save to scratch directory 
            basename = os.path.basename(pattern_summaries % chrom)
            path_tmp = os.path.join(path_scratch, basename)
            summaries.to_csv(path_tmp, sep=SEP, index=False)

            # Add condensed summary file to archive
            archive.add(path_tmp, arcname=basename, recursive=False)

            # Add detection files
            for pm in PM:
                path_detections = pattern_detections % (chrom, pm)
                archive.add(path_detections,
                            arcname=os.path.basename(path_detections),
                            recursive=False)
            
            # Cleanup to avoid excess memory consumption
            del summaries
            gc.collect()
        
        # Finalize the archive
        archive.close()

if __name__ == '__main__':
    main()

