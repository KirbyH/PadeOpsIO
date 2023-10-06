#!/usr/bin/python

# Hannah Johlas, 2019
# Revised Kirby Heck, 2023

import matplotlib
import matplotlib.pyplot as plt

# Plot formatting
# matplotlib.rcParams['figure.figsize']     = [16.0/2.54, 6.0/2.54]               # figure dimensions, in inches
matplotlib.rcParams['savefig.format']     = 'pdf'                               # image file type (pdf, png)
matplotlib.rcParams['savefig.dpi']        = 600                                 # figure resolution (makes a little better .png)
matplotlib.rcParams['savefig.pad_inches'] = 0.10                                # remove extra whitespace
matplotlib.rcParams['savefig.bbox']       = 'tight'                           # remove extra whitespace = 'tight'
matplotlib.rcParams['lines.linewidth']    = 1.0
matplotlib.rcParams['legend.numpoints']   = 1
matplotlib.rcParams['legend.handlelength']  = 1.0
matplotlib.rcParams['legend.handletextpad'] = 0.2
matplotlib.rcParams['legend.columnspacing'] = 1.0
matplotlib.rcParams['legend.borderpad']   = 0.3
matplotlib.rcParams['lines.markersize']   = 2
matplotlib.rcParams['font.family']        = 'serif'
# matplotlib.rcParams['font.serif']         = 'cmr10'
matplotlib.rcParams['xtick.major.pad']    = 2
matplotlib.rcParams['ytick.major.pad']    = 2
matplotlib.rcParams['axes.labelpad']      = 1
matplotlib.rcParams['xtick.labelsize']    = 12
matplotlib.rcParams['ytick.labelsize']    = 12
matplotlib.rcParams['axes.labelsize']     = 12
matplotlib.rcParams['legend.fontsize']    = 12
# matplotlib.rcParams['axes.titlesize']     = 'small'
matplotlib.rcParams['figure.subplot.bottom'] = 0.0
matplotlib.rcParams['figure.subplot.top']    = 1.0
matplotlib.rcParams['figure.subplot.right']  = 1.0
matplotlib.rcParams['figure.subplot.left']   = 0.0
matplotlib.rcParams['text.usetex']           = True
matplotlib.rcParams.update({'figure.autolayout': False})
