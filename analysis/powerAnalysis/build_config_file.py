import copy

import yaml

# Global constants
N_CHROM = 6

# Set default arguments

# Base config file to use
base_cfg_id = 'H_1-combined'
path_base_cfg = 'config/H_1-combined.yml'

# Pattern for output
pattern_output = 'powerAnalysis/data/simChrom_%s.%s'

# Number of replicates
n_replicates = 10

# Path for output configuration
path_cfg = 'powerAnalysis/powerAnalysis.yml'


# Function definitions

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


def main():
    # Load configuration
    with open(path_base_cfg, 'rb') as f:
        cfg = yaml.load(f)
    
    cfg = clean_config(cfg)

    # Change ID for power analysis
    cfg['id'] = 'powerAnalysis'

    # Change data paths
    cfg['data']['chrom_path'] = pattern_output % ('y', 'txt')
    cfg['data']['null_path'] = pattern_output % ('null', 'txt')
    cfg['data']['regions_path'] = pattern_output % ('regions', 'txt')
    cfg['data']['template_path'] = cfg['data']['template_path'].format(
        id=base_cfg_id)
    cfg['data']['n_chrom'] = n_replicates
    
    # Save revised config
    with open(path_cfg, 'wb') as f:
        yaml.dump(cfg, f, default_flow_style=False)


if __name__=='__main__':
    cfg = main()
