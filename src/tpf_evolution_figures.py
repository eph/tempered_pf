import argparse 
import json

import numpy as np 
import pandas as p

from figures import saved_figure, figure_defaults, despine
parser = argparse.ArgumentParser(description='BSPF in the Great Recession Plot')

parser.add_argument('simulation', metavar='sims', type=str, 
                    help='The JSON files with simulation output from tpf_driver using --save-states')
parser.add_argument('--output', action='store', type=str, default='',
                    help='The filename for figure output (default: no filename, plot interactively)')
args = parser.parse_args()

with open(args.simulation) as sim_file:
    sim = json.load(sim_file)


gamQ = sim['inputs']['para'][9]
var_names = ['c', 'ppi', 'R', 'z', 'y', 'g', 'ylag', 'Eppi', 'Ec', 'Ey', 'Ez']

forecast = p.DataFrame(sim['output']['001']['fcst_states']['024'])
forecast.columns = var_names + ['weights']
forecast['ygr'] = forecast.y - forecast.ylag + forecast.z + gamQ

update = p.DataFrame(sim['output']['001']['update_states']['024'])
update.columns = var_names + ['weights']
update['ygr'] = update.y - update.ylag + update.z + gamQ

from scipy.stats import norm
true_filtered_density = norm(loc=-2.977+gamQ, scale=np.sqrt(0.013452427))
xgrid = np.linspace(-6,6,1000)
cp = figure_defaults(grayscale=False)
with saved_figure(args.output) as (fig,ax):

    forecast.ygr.plot(kind='kde', ax=ax, linewidth=5, alpha=0.7, linestyle='dashed')
    despine()
    ax.axvline(x=update.ygr.mean(), linewidth=7, color=cp[1], alpha=0.7)
    ax.plot(xgrid, true_filtered_density.pdf(xgrid), linewidth=5, color=cp[2], alpha=0.7)
    ax.set_xlim(-4,2)
    ax.set_ylabel('')
    ax.legend(['Forecast Density', 'BSPF Filtered Density', 'True Filtered Density'],
              fontsize=18)

