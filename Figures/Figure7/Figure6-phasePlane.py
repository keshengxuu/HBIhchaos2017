# -*- coding: utf-8 -*-
"""

"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker

ISI222=np.loadtxt("data/gsd0222_ISIs.txt")
Vtraces=np.loadtxt("data/gsd0222_Traces.txt.gz")
Vthresh=-6.30248577e+01  #only for gsd=0.222


def plotBranch(ax,sTable,farg='k-',width=2,maxgap=0.05,hilo=False):
    distances=np.sqrt(np.diff(sTable[:,3])**2 + np.diff(sTable[:,6])**2)
    limits=np.where(distances>maxgap)[0] +1 
    limits=np.append(limits,len(sTable))
    limits=np.insert(limits,0,0)
    for inic,fin in zip(limits[:-1],limits[1:]):
#        print inic,fin
        ax.plot(sTable[inic:fin,3],sTable[inic:fin,6],farg,lw=width)
        if hilo:
            ax.plot(sTable[inic:fin,3],sTable[inic:fin,10],farg,lw=width)
    return limits



#%%

plt.figure(1,figsize=(13,4))
plt.clf()

#ax1=plt.subplot2grid((3,2),(0,0))
ax1=plt.subplot2grid((1,7),(0,0),colspan=2)
#ax1=plt.subplot(131)

ax1.plot(Vtraces[64000:74000,-1]/1000,Vtraces[64000:74000,0],'k')
ax1.plot((72,77),(Vthresh,Vthresh),'r-')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Voltage (mV)')


#ax2=plt.subplot2grid((3,2),(0,1))
#ax2=plt.subplot(132)
ax2=plt.subplot2grid((1,7),(0,2),colspan=2)

ax2.plot(ISI222[:,0]/1000,ISI222[:,1],'k.',ms=3)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('intervals (ms)')

#ax3=plt.subplot2grid((3,2),(1,0),rowspan=2,colspan=2,projection='3d')
#ax3=plt.subplot(133,projection='3d')
ax3=plt.subplot2grid((1,7),(0,4),colspan=3,projection='3d')

ax3.plot(Vtraces[:,1],np.zeros_like(Vtraces[:,2]),Vtraces[:,3],lw=0.5,alpha=0.5)
ax3.plot(np.zeros_like(Vtraces[:,1]),Vtraces[:,2],Vtraces[:,3],lw=0.5,alpha=0.5)
ax3.plot(Vtraces[:,1],Vtraces[:,2],np.zeros_like(Vtraces[:,3]),lw=0.5,alpha=0.5)
ax3.plot(Vtraces[:,1],Vtraces[:,2],Vtraces[:,3],'k')

ax3.azim=50
ax3.zaxis.set_rotate_label(False) 
ax3.set_xlabel('$\mathsf{a_{sd}}$',fontsize='x-large')
ax3.set_ylabel('$\mathsf{a_{sr}}$',fontsize='x-large')
ax3.set_zlabel('$\mathsf{a_{h}}$',fontsize='x-large',rotation=90)
ax3.set_xticks((0,0.2,0.4,0.6))
ax3.set_yticks((0,0.2,0.4,0.6))
ax3.set_zticks((0,0.05,0.1,0.15))


plt.figtext(0.01,0.92,'A',fontsize='xx-large')
plt.figtext(0.29,0.92,'B',fontsize='xx-large')
plt.figtext(0.6,0.92,'C',fontsize='xx-large')

plt.tight_layout()