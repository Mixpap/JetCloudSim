import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import configparser
import re
from astropy import units as u
from astropy.constants import c,m_p,k_B
from astropy.visualization import quantity_support
import os


config = configparser.ConfigParser(allow_no_value=True,delimiters=' ')
config.read('pluto.ini'); print('loaded pluto.ini')
Xstring=re.findall(r"[-+]?\d*\.\d+|\d+", config['Grid']['X1-grid'])
Ystring=re.findall(r"[-+]?\d*\.\d+|\d+", config['Grid']['X2-grid'])

XX=np.array([])
j=1
for i in range(int(Xstring[0])):
	XX=np.append(XX,np.linspace(float(Xstring[j]),float(Xstring[j+2]),int(Xstring[j+1])))
	j=j+2

#print(XX,XX.shape)
YY=np.array([])
j=1
for i in range(int(Ystring[0])):
	YY=np.append(YY,np.linspace(float(Ystring[j]),float(Ystring[j+2]),int(Ystring[j+1])))
	j=j+2

#print(YY,YY.shape)
#Xmin=float(Xstring[1]);Xmax=float(Xstring[3]);Nx=int(Xstring[2])
#Ymin=float(Ystring[1]);Ymax=float(Ystring[3]);Ny=int(Ystring[2])
#print(Xstring)
#XX= np.linspace(Xmin,Xmax,Nx)
#YY= np.linspace(Ymin,Ymax,Ny)

dblout=np.loadtxt('dbl.out',str); print('loaded dbl.out')
Nbl=dblout.shape[0]
simvars=dblout[0,6:]
TTall=dblout[:,1].astype(float)

L0=10*u.pc
rho0=1.67e-24*(u.g*u.cm**(-3))
v0=c.cgs
print('Making Pluto Units (PU) for L0 = {} / rho0= {} / v0 ={}'.format(L0,rho0,v0))

gauss_B = (u.g/u.cm)**(0.5)/u.s
equiv_B = [(u.G, gauss_B, lambda x: x, lambda x: x)]
PU={'rho':rho0,
    'vx1':v0,
    'vx2':v0,
    'vx3':v0,
    'Bx1':(v0*np.sqrt(4.*np.pi*rho0)).to(u.G,equivalencies=equiv_B),
    'Bx2':(v0*np.sqrt(4.*np.pi*rho0)).to(u.G,equivalencies=equiv_B),
    'Bx3':(v0*np.sqrt(4.*np.pi*rho0)).to(u.G,equivalencies=equiv_B),
    'prs':(rho0*v0**2).to(u.Ba),
    'tmp':(v0**2*m_p/k_B).cgs,
    'T':(L0/v0).cgs,
    'L':(L0).cgs}

cwd=os.getcwd()
wdir=cwd+'/'
outdir=cwd[cwd.rfind('/')+1:]
d=5
start=0
end=Nbl
A=np.fromfile("data.{:04d}.dbl".format(start))
B=A.reshape((simvars.shape[0],YY.shape[0],XX.shape[0]))
M={var:B[start,:,:] for start,var in enumerate(simvars)}
TT=np.array(TTall[start])
for idbl in range(start+1,end,d):
    A=np.fromfile("data.{:04d}.dbl".format(idbl))
    B=A.reshape((simvars.shape[0],YY.shape[0],XX.shape[0]))
    #M=np.dstack((M,{var:B[i,:,:]*PU[var] for i,var in enumerate(simvars)}))
    M=np.append(M,{var:B[i,:,:] for i,var in enumerate(simvars)})
    TT=np.append(TT,TTall[idbl])
    
RHO= np.array([M[i]['rho'] for i,t in enumerate(TT)])#*PU['rho']
PRS= np.array([M[i]['prs'] for i,t in enumerate(TT)])#*PU['prs']
Br= np.array([M[i]['Bx1'] for i,t in enumerate(TT)])#*PU['Bx1']
Bz= np.array([M[i]['Bx2'] for i,t in enumerate(TT)])#*PU['Bx2']
Bphi= np.array([M[i]['Bx3'] for i,t in enumerate(TT)])#*PU['Bx3']
TMP= np.array([M[i]['tmp'] for i,t in enumerate(TT)])#*PU['tmp']
Vr=np.array([M[i]['vx1'] for i,t in enumerate(TT)])#*PU['vx1']
Vz=np.array([M[i]['vx2'] for i,t in enumerate(TT)])#*PU['vx2']
Vphi=np.array([M[i]['vx3'] for i,t in enumerate(TT)])#*PU['vx3']
np.savez_compressed(outdir,RHO=RHO,PRS=PRS,Br=Br,Bz=Bz,Bphi=Bphi,TMP=TMP,Vr=Vr,Vz=Vz,Vphi=Vphi,T=TT,X=XX,Y=YY,PU=PU)
