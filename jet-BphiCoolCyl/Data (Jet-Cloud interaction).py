# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:hydrogen
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.2'
#       jupytext_version: 1.0.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import numpy as np
#import pyPLUTO as pp
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
rcParams['figure.figsize'] = (17., 14.0)
rcParams['ytick.labelsize'],rcParams['xtick.labelsize'] = 17.,17.
rcParams['axes.labelsize']=19.
rcParams['legend.fontsize']=17.
rcParams['text.latex.preamble'] = ['\\usepackage{siunitx}']
# import seaborn
# seaborn.despine()
# seaborn.set_style('white', {'axes.linewidth': 0.5, 'axes.edgecolor':'black'})
# seaborn.despine(left=True)
import importlib

# %%
#test0=np.load('../test0/test0.npz')
test5=np.load('test5.npz')

# %%
rcParams['figure.figsize'] = (22,22)

# %%
test5['RHO'][-1].shape

# %%
plt.plot(test1['Y'],np.log10(test1['RHO'][-1,:,400]))

# %%
plt.contourf(test5['X'],test5['Y'],np.log10(test5['RHO'][-1,:,:]),cmap='jet',levels=np.linspace(-7,5,13))
plt.colorbar()
plt.xlim(0,7)
plt.ylim(20,40)

# %%
plt.contourf(test5['X'],test5['Y'],np.log10(test5['TMP'][-1,:,:]),cmap='jet',levels=np.linspace(0,16,17))
plt.colorbar()
plt.xlim(0,7)
plt.ylim(20,40)

# %%
PU=test5['PU']

# %%
PU

# %%
gamma=1./(1.-(test5['Vz']**2+test5['Vr']**2+test5['Vphi']**2))

# %%
V=np.sqrt(test5['Vz']**2+test5['Vr']**2+test5['Vphi']**2)

# %%
V.max()

# %%
plt.contourf(test5['X'],test5['Y'],np.sqrt(test5['Vz'][70,:,:]**2+test5['Vr'][70,:,:]**2+test5['Vphi'][70,:,:]**2),cmap='jet')#,levels=np.linspace(0,16,17))
plt.colorbar()
plt.xlim(0,4)
plt.ylim(20,40)

# %%
gamma.shape

# %%
gamma[:,:,:].max()

# %%
plt.contourf(test5['X'],test5['Y'],gamma[-1,:,:],cmap='jet',levels=np.arange(0,17,0.25))
plt.colorbar()
plt.xlim(0,2.)
#plt.ylim(20,40)

# %%
plt.contourf(test5['X'],test5['Y'],np.log10(np.abs(test5['Bphi'][-1,:,:])),cmap='jet',levels=np.arange(-4,-2,0.25))
plt.colorbar()
plt.contour(test5['X'],test5['Y'],gamma[-1,:,:],cmap='jet',levels=np.arange(0,17,0.5))

plt.xlim(0,2.)
#plt.ylim(20,40)

# %%
sigma=test1['Bphi']**2/(gamma**2*test1['RHO'])

# %%
plt.contourf(test1['X'],test1['Y'],2*np.log10(sigma[-1,:,:]),cmap='jet',levels=np.linspace(-50,0,17))
plt.colorbar()
plt.xlim(0,4)
plt.ylim(20,40)

# %%
plt.contourf(test1['X'],test1['Y'],2*np.log10(test1['Bphi'][-1,:,:]),cmap='jet',levels=np.linspace(-10,0,17))
plt.colorbar()
plt.xlim(0,4)
plt.ylim(20,40)

# %%
importlib.reload(f)
f.quadruple(test1,np.log10(test1['RHO']),rows=2,cols=3,tdk='yrs',Save_Figure='',scale=(6,10),color='jet',zlim=[-1,3])

# %%
test1['T']

# %%
importlib.reload(f)
f.quadruple(test1,np.log10(test1['RHO']),rows=2,cols=3,tdk='yrs',color='jet',scale=(6,10),zlim=[-1,3])

# %%
beta=test1['PRS']/test1['Bz']**2

# %%
beta.min(),beta.max()

# %%
rcParams['figure.figsize'] = (12,22)
T=-1
#plt.contourf(np.log10(test1['RHO'][:,:,T]).T,origin='lower',cmap='Greys',levels=np.linspace(-3,5,50))
plt.contourf(np.log10(test1['TMP'][:,:,T]).T,origin='lower',cmap='jet',levels=np.linspace(2,15,100))
plt.contour(test1['Bz'][:,:,T].T**2,origin='lower',linewidth=7)
plt.contour(test1['PRS'][:,:,T].T**2,origin='lower',linestyles='--',linewidth=7)
plt.xlim(0,150)
plt.ylim(100,300)

# %%
rcParams['figure.figsize'] = (12,22)
plt.contourf(test1['Bz'][:,:,-1].T**2,origin='lower',cmap='Greys')
plt.contour(np.log10(test1['RHO'][:,:,-1]).T,origin='lower',linewidth=5,levels=np.linspace(-2,3,3))

# %%
f.quadruple(test1,test1['Bz']**2,rows=2,cols=3,color='Greys',scale=(6,10),tdk='yrs',Save_Figure='')

# %%
test1['Vy'].shape

# %%
f.quadruple(test1,np.log10(np.sqrt(test1['Vy']**2+test1['Vx']**2)),rows=2,cols=3,color='jet',scale=(6,15),tdk='kyrs',zlim=[-10,0],Save_Figure='beta')

# %%
f.quadruple(test1,np.log10(beta),rows=2,cols=3,color='jet',scale=(6,10),tdk='yrs',zlim=[-5,5],Save_Figure='beta')

# %%
f.quadruple(test1,np.log10(test1['PRS']),rows=2,cols=3,tdk='yrs')

# %%
test1['Bz'].shape

# %%
rcParams['figure.figsize'] = (12,12)
plt.plot(test1['Bz'][:,0,1]**2)
plt.plot(test1['PRS'][:,0,1])
#plt.plot(test1['X'][:],test1['Bz'][:,40,1]**2)
#plt.plot(test1['X'][:],test1['Bz'][:,50,1]**2)
plt.yscale('log')
plt.ylim(1e-10,1e-1)
plt.xlim(0,100)

# %%
rcParams['figure.figsize'] = (12,12)
plt.plot(test1['X'][:],test1['Bz'][:,0,1]**2)
plt.plot(test1['X'][:],test1['PRS'][:,0,1])
#plt.plot(test1['X'][:],test1['Bz'][:,40,1]**2)
#plt.plot(test1['X'][:],test1['Bz'][:,50,1]**2)
plt.yscale('log')
plt.ylim(1e-10,1e-1)

# %%
plt.plot(test1['X'][:],np.log10(beta[:,1,2]),'--')
plt.ylim(-2,2)

# %%
rcParams['figure.figsize'] = (12,12)
plt.plot(test1['X'][:],test1['Bz'][:,0,1]**2)
plt.plot(test1['X'][:],test1['PRS'][:,0,1])
#plt.plot(test1['X'][:],test1['Bz'][:,50,1]**2)
plt.yscale('log')
plt.ylim(1e-10,1e-1)

# %%
plt.plot(test1['X'][:],test1['Bz'][:,0,1])
# plt.plot(test1['X'][:],test1['Bz'][:,40,1])
# plt.plot(test1['X'][:],test1['Bz'][:,50,1])
#plt.yscale('log')

# %%
np.log10(test1['TMP'].max())

# %%
plt.plot(test1['X'][:],test1['TMP'][:,0,0])
plt.plot(test1['X'][:],test1['TMP'][:,100,3])
#plt.plot(test1['X'][:],test1['TMP'][:,0,-1])
plt.yscale('log')


# %%
plt.plot(test1['X'][:],test1['TMP'][:,0,0])
plt.plot(test1['X'][:],test1['TMP'][:,300,3])
#plt.plot(test1['X'][:],test1['TMP'][:,0,-1])
plt.yscale('log')


# %%
plt.plot(test1['X'][:],test1['PRS'][:,0,0])
plt.plot(test1['X'][:],test1['PRS'][:,0,1])
#plt.plot(test1['X'][:],test1['TMP'][:,0,-1])
plt.yscale('log')


# %%
plt.plot(test1['X'][:],test1['RHO'][:,0,0])
plt.plot(test1['X'][:],test1['RHO'][:,0,1])
#plt.plot(test1['X'][:],test1['TMP'][:,0,-1])
plt.yscale('log')

# %%
importlib.reload(f)
f.pprofile(test1,'TMP',tdk='yrs',steps=3,yprop=0,ix=0)

# %%
