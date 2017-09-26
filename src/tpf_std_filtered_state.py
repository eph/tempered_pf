import argparse 
import json

import numpy as np 
import pandas as p

import tqdm
from figures import saved_figure, figure_defaults, despine


parser = argparse.ArgumentParser(description='Std of filtered state')

parser.add_argument('simulations', metavar='sims', type=str, nargs="+",
                    help='The JSON files with simulation output from tpf_driver using --save-states')
parser.add_argument('--output', action='store', type=str, default='',
                    help='The filename for figure output (default: no filename, plot interactively)')


var_names = ['c', 'ppi', 'R', 'z', 'y', 'g', 'ylag', 'Eppi', 'Ec', 'Ey', 'Ez']

stds = []
means = []
args = parser.parse_args()

for sim_file in tqdm.tqdm(args.simulations):
    with open(sim_file) as sim:
        results = json.load(sim)

    if results['inputs']['sample'] == 'great_moderation':
        index = p.period_range(freq='Q', start='1983Q1', periods=80)
    else:
        index = p.period_range(start='2003Q1', freq='Q', periods=44)


    reps = filter(lambda x: x.isdigit(), results['output'].keys())

    T = len(index)
    sim_results = p.DataFrame()
    for r in reps:
        time_series = [results['output'][r]['mean_filtered_states']['{:03d}'.format(n+1)]
                       for n in range(T)]

        res_r = p.DataFrame(np.array(time_series), index=index, columns=var_names)
        res_r['sim'] = int(r)
        sim_results = sim_results.append(res_r)

    stds.append(sim_results.groupby(sim_results.index).std())
    means.append(sim_results.groupby(sim_results.index).mean())

cp = figure_defaults(n=5)
with saved_figure(args.output) as (fig,ax):
    for i, std in enumerate(stds):
        std.g.plot(ax=ax, linewidth=5, color=cp[i])

    ax.legend(['BSPF, M=40k', 'TPF(r*=2), M=40k', 'TPF(r*=3), M=40k',
               'TPF(r*=2), M=4k', 'TPF(r*=3), M=4k'], loc='lower right', ncol=2)

    despine()

