import numpy as np
import pyPLUTO as pp
from astropy.io import ascii
import os
import sys
from ipywidgets import interactive, widgets,fixed
from IPython.display import Audio, display
%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.animation import FuncAnimation,FFMpegWriter
from matplotlib import rc,rcParams
from scipy.integrate import quad

rc('text', usetex=True)
rcParams['figure.figsize'] = (15., 6.0)
rcParams['ytick.labelsize'],rcParams['xtick.labelsize'] = 17.,17.
rcParams['axes.labelsize']=19.
rcParams['legend.fontsize']=17.
rcParams['text.latex.preamble'] = ['\\usepackage{siunitx}']
# import seaborn
# seaborn.despine()
# seaborn.set_style('white', {'axes.linewidth': 0.5, 'axes.edgecolor':'black'})
# seaborn.despine(left=True)
import importlib
cd Notebooks/
import f
test1=np.load('../test1/test1.npz')
mkdir Images

importlib.reload(f)
f.quadruple(test1,np.log10(test1['TMP']),rows=1,cols=3,tdk='yrs',Save_Figure='test1')
f.quadruple(test1,np.log10(test1['RHO']),rows=1,cols=3,tdk='yrs',Save_Figure='test1_rho')

np.log10(test1['TMP']).max()
np.log10(test1['TMP']).max()
np.log10(test1['TMP'][:,:,-1]).sum()
np.log10(test1['TMP'][:,:,-1]).sum()
np.log10(test1['TMP'][:,:,-1]).sum()
