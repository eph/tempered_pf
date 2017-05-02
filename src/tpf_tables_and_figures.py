import argparse
import json

import numpy as np
import pandas as p

parser = argparse.ArgumentParser(description="Create Figures and "
                                 "Tables for 'Tempered Particle Filtering'")

parser.add_argument('simulations', metavar='sims', type=str, nargs='+',
                    help='The JSON files with simulation output from tpf_driver.')

parser.add_argument('--output', action='store', type=str, default='',
                    help='The filename for figure output (default: no filename, plot interactively)')
# parser.add_argument('--figure-details', action='store', default=None,
#                     help='Flag to call additional (hardcoded) figure details')

args = parser.parse_args()
print(args.simulations)


sims = []

for sim_file in args.simulations:
    with open(sim_file) as sim_file:
        sim = json.load(sim_file)

    bias_series = p.Series(sim['output']['likhat']) - sim['output']['truth']
    sim['output']['bias_series'] = bias_series
    sim['output']['bias1'] = bias_series.mean()
    sim['output']['std1'] = bias_series.std()
    sim['output']['bias2'] = (np.exp(bias_series)-1).mean()

    if sim['inputs']['filter'] == 'bootstrap':
        sim['inputs']['name'] = 'BSPF' 
    elif sim['inputs']['filter'] == 'resample':
        sim['inputs']['name'] = 'resample' 
    elif sim['inputs']['filter'] == 'opt':
        sim['inputs']['name'] = 'opt' 
    else:
        tpf_string = 'TPF$(r^*={:1.0f})$'
        sim['inputs']['name'] = tpf_string.format(sim['inputs']['rstar'])
                                 
                                 
    sims.append(sim)
    

fl = '{: 7.3f}'.format
inl = '{: 7d}'.format

rows = [['                                 '] + [sim['inputs']['name'] for sim in sims],
        ['Number of Particles              '] + [inl(sim['inputs']['npart']) for sim in sims],
        ['Number of Repetitions            '] + [inl(sim['inputs']['nsim']) for sim in sims],
        ['Bias $\hat{\Delta}_1$            '] + [fl(sim['output']['bias1']) for sim in sims],
        ['StdD $\hat{\Delta}_1$            '] + [fl(sim['output']['std1']) for sim in sims],
        ['Bias $\hat{\Delta}_2$            '] + [fl(sim['output']['bias2']) for sim in sims],
        [r'$T^{-1}\sum_{t=1}^{T}N_{\phi,t}$ '] + [fl(np.mean(sim['output']['avg_iterations'])) for sim in sims],
        ['Average Run Time (s)             '] + [fl(sim['output']['average_time']) for sim in sims],
]


table = ' \\\\ \n'.join([' & '.join(r) for r in rows]) + ' \\\\'
print(table)





