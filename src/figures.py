from contextlib import contextmanager
import matplotlib.pyplot as plt

from matplotlib import rc

import seaborn as sns

import numpy as np

despine = sns.despine

preamble = r"""
\usepackage[helvet]{sfmath}
\usepackage[scaled]{helvet}
\renewcommand\familydefault{\sfdefault}
\usepackage[T1]{fontenc}
"""
def figure_defaults(n=4, grayscale=False):
    #rc('text', usetex=True)
    #rc('text.latex', preamble=preamble)
    #rc('lines', linewidth=3)

    if grayscale:
        cp = sns.color_palette("gray")
    else:
        cp = sns.color_palette("cubehelix", n)

    sns.set_palette(cp)
    sns.set_style('white')
    sns.set_context(rc={'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': True, 
                            'font.family':'serif', 'font.sans-serif':[u'Helvetica']})

    return sns.color_palette()

@contextmanager
def saved_figure(fname='', **kwds):
    """
    Saves a figure in `fname`.
    """
    fig, ax = plt.subplots(**kwds)
    try:
        yield (fig, ax)
    finally:
        if fname=='':
            plt.show()
        else:
            fig.savefig(fname, bbox_inches='tight')
            plt.close(fig)

        
        
def add_arrows(ax, scaling=0.008, xarrow=None, yarrow=None, xscaling=None, yscaling=None):
    """
    Adds arrows to x and y axes. 
    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    xspan = xmax - xmin
    yspan = ymax - ymin

    if xscaling is None:
        xscaling = scaling
    if yscaling is None:
        yscaling = scaling

    if xarrow is None:
        xarrow = xscaling*xspan, xscaling*yspan

    if yarrow is None:
        yarrow = yscaling*yspan, yscaling*yspan

    ax.arrow(xmin, ymin, xmax-xmin, 0, color='k', lw=0,
             head_length=xarrow[0], head_width=4*xarrow[1], clip_on=False)

    ax.arrow(xmin, ymin, 0, ymax-ymin, color='k', lw=0,
             head_length=4*yarrow[0], head_width=yarrow[1], clip_on=False)
