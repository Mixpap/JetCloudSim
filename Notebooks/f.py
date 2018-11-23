import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
#import seaborn
#import yt
import pyPLUTO as pp
from astropy.io import ascii
import os
import sys
from ipywidgets import interactive, widgets,fixed
from IPython.display import Audio, display
import matplotlib.ticker as ticker
from matplotlib.animation import FuncAnimation,FFMpegWriter
from matplotlib import rc,rcParams
rc('text', usetex=True)
rcParams['text.latex.preamble'] = ['\\usepackage{siunitx}']

global mp
mp=1.67e-24
global kb
kb=1.38e-16
global Gcgs
Gcgs = 6.6743e-8
global RHO0
RHO0=1.67e-24
global V0
V0=3e10
global L0
L0=3e19
global PRS0
PRS0=RHO0*V0**2
global t0
t0=L0/V0
global Temp0
Temp0=V0**2*(mp/kb)
global B0
B0=np.sqrt(4.*np.pi*RHO0*V0**2)
global G0
G0=Gcgs*(RHO0*t0**2)
global M0
M0=RHO0*L0**3
global dt


def MakeGif(VAR,X,Y,Vx,Vy,Bx,By,VARc,vfile,Contours=False,qV=False,qB=False,step=1,v1label='',v2label='',nn=10):
        T=np.arange(0,VAR.shape[2],1)
        t=T[::step]
        fig=plt.figure(figsize=(10,10))
        fig.set_tight_layout(True)
        ext=[X.min(),X.max(),Y.min(),Y.max()]
        ax1=plt.subplot()
        ax1.get_yaxis().get_major_formatter().set_useOffset(False)
        #ax1.add_artist(plt.Circle((0, 0), 1.0, color='r',fill=False))
        pc = ax1.imshow(VAR[:,:,t[0]].T,cmap='viridis',origin='lower',aspect='equal',extent=ext,
                            vmin=VAR.min(axis=2).min(),vmax=VAR.max(axis=2).max())
        if VARc!='':
                pcc=pc.axes.contour(X,Y,VARc[:,:,t[0]].T,origin='lower',aspect='equal',extent=ext,
                                    levels=np.linspace(VARc.min(axis=2).min(),VARc.max(axis=2).max(),8))
                cbc=plt.colorbar(pcc,ax=ax1,fraction=0.046, pad=0.04)
                cbc.formatter.set_powerlimits((0, 0))
                cbc.update_ticks()
        cb=plt.colorbar(pc,ax=ax1,fraction=0.046, pad=0.04)
        cb.formatter.set_powerlimits((0, 0))
        cb.update_ticks()
        def update(i):
            ax1.cla()
            pc = ax1.imshow(VAR[:,:,i].T,cmap='viridis',origin='lower',aspect='equal',extent=ext,
                            vmin=VAR.min(axis=2).min(),vmax=VAR.max(axis=2).max())
            if Contours:
                pcc=pc.axes.contour(X,Y,VARc[:,:,i].T,origin='lower',aspect='equal',extent=ext,
                                    levels=np.linspace(VARc.min(axis=2).min(),VARc.max(axis=2).max(),8))
            if qV:
                pc.axes.quiver(X[::nn],Y[::nn],Vx[:,:,i][::nn,::nn].T,Vy[:,:,i][::nn,::nn].T)
            if qB:
                pc.axes.quiver(X[::nn],Y[::nn],Bx[:,:,i][::nn,::nn].T,By[:,:,i][::nn,::nn].T,color='r')
            label = 'Time = {0:.1f} kyrs'.format(i*dt/1000.)
            ax1.set_title(label)
            return ax1
        anim = FuncAnimation(fig, update, frames=range(VAR.shape[2]), interval=200)
        anim.save(vfile+'.gif',writer='imagemagic',bitrate=-1,dpi=72,codec='libx264',extra_args=['-pix_fmt', 'yuv420p'])
        print("ffmpeg -f gif -i "+vfile+".gif "+vfile+".mp4")
        os.system("ffmpeg -f gif -i "+vfile+".gif "+vfile+".mp4")
def div(F):
    """ compute the divergence of n-D scalar field `F` """
    return reduce(np.add,np.gradient(F))
def grad(F):
    return np.sqrt(np.gradient(F)[0]**2+np.gradient(F)[1]**2)

def quadruple(d,VAR,tdk='Myrs',Save_Figure='',cl='',nn=0,mspeed='km',rows=2,cols=2,scale=(6,6),xlim=[None,None],
              ylim=[None,None],zlim=[None,None],color='viridis',datafolder='Images/'):
    """
    Plot a rows(=2) x cols(=2) Variable
    """
    X,Y=d['X'],d['Y']
    Vx=d['Vx'] if nn>0 else 0
    Vy=d['Vy'] if nn>0 else 0
    T=np.linspace(0,d['T'].shape[0]-1,rows*cols,dtype=int)
    fig, axes = plt.subplots(nrows=rows, ncols=cols, sharex=True, sharey=True,
                            figsize=(cols*scale[0],rows*scale[1]))
    td=1e3 if tdk=='kyrs' else 1e6
    for i,ax in enumerate(axes.flat):
        ext=[X.min(),X.max(),Y.min(),Y.max()]
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        #ax.add_artist(plt.Circle((0, 0), 1.0, color='r',fill=False,linestyle='--'))
        label = '{:.1f} {}'.format(d['T'][T[i]]/td,tdk)
        ax.set_title(label,fontsize=20)
        ax.grid(False)
        zmin=VAR.min() if zlim[0]==None else zlim[0]
        zmax=VAR.max() if zlim[1]==None else zlim[1]
        pc = ax.imshow(VAR[:,:,T[i]].T,cmap=color,origin='lower',aspect='equal',
                       extent=ext,vmin=zmin,vmax=zmax)
        if nn>0:
            k=nn #distance from boundaries for first/last arrows
            sc=2. if mspeed =='max' else 5. if mspeed == 'c' else 1e-4
            q=pc.axes.quiver(X[k:-k:nn],Y[k:-k:nn],
                            Vx[:,:,T[i]][k:-k:nn,k:-k:nn].T,
                            Vy[:,:,T[i]][k:-k:nn,k:-k:nn].T,
                             scale=sc,alpha=0.5,width=0.002)
            if mspeed == 'c':
                pc.axes.quiverkey(q,0.05,1.02,1.,r'$1\si{c}$',labelpos='E',fontproperties={'weight': 'bold'})
            elif mspeed == 'max':
                mV=np.max(np.sqrt(Vx[np.argmin((d['Y']-ylim[0])**2):np.argmin((d['Y']-ylim[1])**2),
                                     np.argmin((d['X']-xlim[0])**2):np.argmin((d['X']-xlim[1])**2),T[i]]**2+
                                  Vy[np.argmin((d['Y']-ylim[0])**2):np.argmin((d['Y']-ylim[1])**2),
                                     np.argmin((d['X']-xlim[0])**2):np.argmin((d['X']-xlim[1])**2),T[i]]**2))
                pc.axes.quiverkey(q,0.05,1.02,mV,'{:.2f} c'.format(mV),labelpos='E',
                                  fontproperties={'weight': 'bold'})
            else:
                pc.axes.quiverkey(q,0.02,1.02,3.36e-6,r'$1\si{km.s^{-1}}$',labelpos='E',fontproperties={'weight': 'bold'})

    ax.set_xlim(xlim[0],xlim[1])
    ax.set_ylim(ylim[0],ylim[1])
    plt.tight_layout()
    cbar_ax = fig.add_axes([0., 1.015, 1., 0.025*(np.float(cols)/rows)])#*(np.float(cols)/rows)
    cb=fig.colorbar(pc, cax=cbar_ax,orientation="horizontal",label=cl)
    cb.ax.tick_params(labelsize=17)
    cb.ax.xaxis.offsetText.set(size=20)
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
    if (Save_Figure!=''):
        plt.savefig(datafolder+Save_Figure,bbox_inches='tight')

def RadiusPlot(c,tdk='Myrs',radscale=80./256.,mid=256/2,Save_Figure='',datafolder='../Document/DataImages/'):
    td=1e6 if tdk=='Myrs' else 1e3
    plt.figure(figsize=(7.5,6.5))
    Radius=np.array([])
    VAR=c['RHO']
    #T=np.linspace(0,c['T'].shape[0]-1,1,dtype=int)
    for t in range(VAR.shape[2]):
        Radius=np.append(Radius,mid-np.abs(np.diff(VAR[mid,:,t])).argmax())
    plt.ylabel('Radius (pc)',fontsize=20)
    plt.xlabel('Time {}'.format(tdk),fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.plot(c['T']/td,Radius*radscale)
    plt.tight_layout()
    if (Save_Figure!=''):
        plt.savefig(datafolder+Save_Figure,bbox_inches='tight')


def plutoplot2D(w,T,
                Velocity=False,Magnetic_Field=False,log10=False,Div=False,
                show=['Density'],Contours=['Density'],
                VelocityD='0',x=str(128),y=str(128),Save_Figure='',MakeVideo=''):
    X=w['X']
    Y=w['Y']
    x=int(x)
    y=int(y)
    VAR=np.log10(w[show]) if log10 else np.abs(div(w[show])) if Div else w[show]
    VARt=VAR[x,y,:]
    dt=w['T'][1]-w['T'][0]
    td=1.
    if Contours!='None':
        VARc=np.log10(w[Contours]) if log10 else np.abs(div(w[Contours])) if Div else w[Contours]
        VARct=VARc[x,y,:]
    else:
        VARc=0.0
        VARct=0.0
    fig =plt.figure(figsize=(12,10))
    fig.set_tight_layout(True)
    gs = gridspec.GridSpec(12, 10)
    ax,pc = pimshow(gs[:10,:],VAR,X,Y,T,dt=dt)
    if Contours!='None':
        axc,pcc = pimcontour(gs[:10,:],VARc,X,Y,T,pc)
        #plt.colorbar(pcc,ax=ax,fraction=0.046, pad=0.04)
    axp1=plt.subplot(gs[10:,:])

    x,y=int(x),int(y)
    #print x,X[x],y,Y[y],VAR[x,y,T],VARc[x,y,T]
    ax.plot(X[x],Y[y],'o',color='r')
    TT=np.arange(0,VARt.shape[0],1)*dt/td
    axp1.plot(TT,VARt)
    axp1.vlines(T*dt/td,VARt.min(),VARt.max(),linestyle='--')
    axp1.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    axp1.set_ylabel(show)
    if Contours!='None':
        axp2=axp1.twinx()
        axp2.plot(TT,VARct,linestyle='--')
        axp2.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
        axp2.set_ylabel(Contours)
    if Velocity:
        Vx,Vy=w['Vx'],w['Vy']
        nn=int(VelocityD)
        if (nn>0):
            q=pc.axes.quiver(X[::nn],Y[::nn],Vx[:,:,T][::nn,::nn].T,Vy[:,:,T][::nn,::nn].T,
                    scale=1e-4,alpha=0.5,width=0.002)
            pc.axes.quiverkey(q,0.02,1.02,3.36e-6,r'$1\si{km.s^{-1}}$',labelpos='E',
                    fontproperties={'weight': 'bold'})
    if Magnetic_Field:
        Bx,By=w['Bx'],w['By']
        nn=int(VelocityD)
        if (nn>0):
            pc.axes.quiver(X[::nn],Y[::nn],Bx[:,:,T][::nn,::nn].T,By[:,:,T][::nn,::nn].T,color='r')
    if (Save_Figure!=''):
        plt.tight_layout()
        plt.savefig(Save_Figure)
    if (MakeVideo!=''):
        MakeGif(VAR,X,Y,Vx,Vy,Bx,By,VARc,MakeVideo,Contours=Contours,
                qV=Velocity,qB=Magnetic_Field,v1label=show,v2label=Contours,step=20)

def pprofile(c,VAR,VAR2=None,steps=5,itlim=-1,ix=128,yprop=128,tdk='Myrs',yl='n $(\si{cm^{-3}})$',yl2='Pressure',
              alpha=1.0,alpha2=1.0,sc1='log',sc2='linear',secopt='left',Save_Figure='',xlim=[None,None],ylim=[None,None],):
    
    if (VAR2!= None and secopt=='left'): 
        plt.figure(figsize=(12,5.2))
    else: #ax.set_aspect('equal')
        plt.figure(figsize=(7.5,6.5))
    ax=plt.subplot()
    datafolder='../Document/DataImages/'
    T=np.linspace(0,c['T'][:itlim].shape[0]-1,steps,dtype=int)
    d=c[VAR]
    if VAR2!=None:
        ax2=ax.twinx()
        ax2.set_ylabel(yl2)
        if( secopt=='left'):
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
            ax2.yaxis.tick_left()
            ax2.yaxis.set_label_position("left")
            ax.vlines(0,0.8*d.min(),1.4*d.max(),linewidth=0.5)
            ls='-'
        else:
            ls='--'
        d2=c[VAR2] 
        ax2.set_yscale(sc2)
    td= 1e6 if tdk == 'Myrs' else 1e3 if tdk == 'kyrs' else 1.
    for t in T:
        ax.plot(c['X'][ix:]*10.,d[yprop,ix:,t],label='{:.1f} {}'.format(c['T'][t]/td,tdk),alpha=alpha)
        if VAR2!=None:
            if secopt=='left': 
                ax2.plot(c['X'][:ix]*10,d2[yprop,:ix,t],linestyle=ls,label='{:.1f} {}'.format(c['T'][t]/td,tdk),
                                alpha=alpha2)
            else:
                ax2.plot(c['X'][ix:]*10,d2[yprop,ix:,t],linestyle=ls,label='{:.1f} {}'.format(c['T'][t]/td,tdk),
                                alpha=alpha2)
    if steps<8: 
        ax.legend(loc='best')
    ax.set_ylim(0.8*d.min(),1.4*d.max()) 
    if VAR2!=None:
        ax2.set_xlim(xlim[0],xlim[1])
        ax2.set_ylim(ylim[0],ylim[1])
        if secopt=='left':
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
    ax.set_xlabel('X $(\si{pc})$')
    ax.set_ylabel(yl)
    ax.set_yscale(sc1)
    ax.set_xlim(xlim[0],xlim[1])
    ax.set_ylim(ylim[0],ylim[1])
    plt.tight_layout()
    if Save_Figure != '': plt.savefig(datafolder+Save_Figure,bbox_inches='tight')