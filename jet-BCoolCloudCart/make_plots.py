import numpy as np
from astropy.io import ascii
import astropy.units as u
import os
import sys
#from ipywidgets import interactive, widgets,fixed
#from IPython.display import Audio, display
#%matplotlib inline

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.animation import FuncAnimation,FFMpegWriter, PillowWriter
from matplotlib import rc,rcParams

#from celluloid import Camera
from scipy.integrate import quad
rc('text', usetex=True)
rcParams['figure.figsize'] = (17., 14.0)
rcParams['ytick.labelsize'],rcParams['xtick.labelsize'] = 17.,17.
rcParams['axes.labelsize']=19.
rcParams['legend.fontsize']=17.
rcParams['text.latex.preamble'] = ['\\usepackage{siunitx}']
# rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
# rc('animation', html='html5')

cwd=os.getcwd()
wdir=cwd+'/'
outdir=cwd[cwd.rfind('/')+1:]
sim=np.load('{}.npz'.format(outdir))
print('Loaded Simulation of size {}'.format(sim['RHO'].shape))
PU=sim['PU']
print('Loaded simulation units')

for var in ['RHO','TMP','Vz','Bphi']:
    print('{}: min={:.2e} mean={:.2e} max={:.2e}'.format(var,sim[var].min(),sim[var].mean(),sim[var].max()))
V=np.sqrt(sim['Vz']**2+sim['Vr']**2+sim['Vphi']**2)
gamma=1./np.sqrt(1.-V**2); print('Gamma min:{:.2f} mean:{:.2f} max:{:.2f}'.format(gamma.min(),gamma.mean(),gamma.max()))
sigma=sim['Bphi']**2/(gamma**2*sim['RHO']) ; print('Sigma min:{:.2e} mean:{:.2e} max:{:.2e}'.format(sigma.min(),sigma.mean(),sigma.max()))

def make_plot(TT,save=True,animate=False,snap=False,x0=0,y0=0,dx=2,dy=2):
    plt.clf()
    fig,ax=plt.subplots(ncols=2,nrows=2,figsize=(30,30))
    [[ax0,ax1],[ax2,ax3]] = ax
    if snap: camera = Camera(fig)
    for T in TT:
        print('Time: {} [{}])'.format(sim['T'][T],T))
        # for var in ['RHO','TMP','Vz','Bphi']:
        #     print('{}: min={:.2f} mean={:.2f} max={:.2f}'.format(var,sim[var][T,:,:].min(),sim[var][T,:,:].mean(),sim[var][T,:,:].max()))
        # print('gamma: min={:.2f} mean={:.2f} max={:.2f}'.format(gamma[T,:,:].min(),gamma[T,:,:].mean(),gamma[T,:,:].max()))
        # print('sigma: min={:.2f} mean={:.2f} max={:.2f}'.format(sigma[T,:,:].min(),sigma[T,:,:].mean(),sigma[T,:,:].max()))

        im0 = ax0.contourf(sim['X'],sim['Y'],np.log10(sim['RHO'][T,:,:]),origin='lower',cmap='jet',levels=np.linspace(-5,5,50))
        ax0.set_title('Density')
        ax0.set_xlim(x0,x0+dx)
        ax0.set_ylim(y0,y0+dy)
        divider = make_axes_locatable(ax0)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im0, cax=cax, orientation='vertical')

        im1=ax1.contourf(sim['X'],sim['Y'],np.log10(sim['TMP'][T,:,:]),origin='lower',cmap='jet',levels=np.linspace(0,17,50))
        ax1.set_title('Temperature')
        ax1.set_xlim(x0,x0+dx)
        ax1.set_ylim(y0,y0+dy)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im1, cax=cax, orientation='vertical');

        im2=ax2.contourf(sim['X'],sim['Y'],gamma[T,:,:],origin='lower',cmap='jet',levels=np.linspace(0,5,50))
        ax2.set_title('Gamma')
        ax2.set_xlim(x0,x0+dx)
        ax2.set_ylim(y0,y0+dy)
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im2, cax=cax, orientation='vertical');

        # im3=ax1.contourf(sim['X'],sim['Y'],np.log10(sigma[T,:,:]),origin='lower',cmap='Greys',levels=np.linspace(0,17,50))
        # ax3.set_title('Gamma')
        # divider = make_axes_locatable(ax3)
        # cax = divider.append_axes('right', size='5%', pad=0.05)
        # fig.colorbar(im3, cax=cax, orientation='vertical');

        im3=ax3.contourf(sim['X'],sim['Y'],sim['Vz'][T,:,:],origin='lower',cmap='RdBu',levels=np.linspace(-1,1,52))
        ax3.set_title('Vz')
        ax3.set_xlim(x0,x0+dx)
        ax3.set_ylim(y0,y0+dy)
        divider = make_axes_locatable(ax3)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im3, cax=cax, orientation='vertical');
        if snap: camera.snap()
        if save: plt.savefig('simplot{:03d}'.format(T), bbox_inches='tight')
        #plt.cla()
    if snap:
        animation = camera.animate()
        #animation.save('animation.gif')#,writer=PillowWriter(fps=24)) doesn work
    if animate:
        os.system('ffmpeg -r 60 -s 1920x1080 -i simplot%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p simanim.mp4')

make_plot([60,70,75,80],x0=0,y0=37,dx=5,dy=8,snap=False,animate=False,save=True)
