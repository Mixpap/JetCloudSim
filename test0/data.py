import pyPLUTO as pp
import numpy as np
import os
t0=1000000000.0
second_to_yrs=3.17e-8
cwd=os.getcwd()
wdir=cwd+'/'
outdir=cwd[cwd.rfind('/')+1:]
print wdir,outdir
nlinf = pp.nlast_info(w_dir=wdir)
D=pp.pload(0,w_dir=wdir)
PRS=D.prs
RHO=D.rho
X=D.x1
Y=D.x2
Vx=D.vx1
Vy=D.vx2 if 'vx2' in D.vars else 0.
Bx=D.bx1 if 'bx1' in D.vars else 0.
By=D.bx2 if 'bx2' in D.vars else 0.
x_HI=D.X_HI if 'X_HI' in D.vars else 0.
x_HII=D.X_HII if 'X_HII' in D.vars else 0 
x_H2=D.X_H2 if 'X_H2' in D.vars else 0 
T=D.SimTime

for t in range(1,nlinf['nlast']):
	D=pp.pload(t,w_dir=wdir)
	PRS=np.dstack((PRS,D.prs))
	RHO=np.dstack((RHO,D.rho))
	Vx=np.dstack((Vx,D.vx1))
	Vy=np.dstack((Vy,D.vx2))
	Bx=np.dstack((Bx,D.bx1)) if 'bx1' in D.vars else 0.
	By=np.dstack((By,D.bx2)) if 'bx2' in D.vars else 0.
	x_HI=np.dstack((x_HI,D.X_HI)) if 'X_HI' in D.vars else 0.
	x_HII=np.dstack((x_HII,D.X_HII)) if 'X_HII' in D.vars else 0.
	x_H2=np.dstack((x_H2,D.X_H2)) if 'X_H2' in D.vars else 0.
	T=np.append(T,D.SimTime)
	
np.savez_compressed(outdir,RHO=RHO,PRS=PRS,Vx=Vx,Vy=Vy,Bx=Bx,By=By,X=X,Y=Y,
x_HI=x_HI,x_HII=x_HII,x_H2=x_H2,T=T*t0*second_to_yrs)
print nlinf['time'],nlinf['nlast'],T[-1]
#T=np.linspace(0.,nlinf['time'],nlinf['nlast'])*t0*second_to_yrs
