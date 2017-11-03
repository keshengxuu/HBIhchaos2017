# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 09:22:20 2016

@author: porio
"""
from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

isis = np.loadtxt("data/ISIs_gsdT36gh02.txt.gz")

Vtraces = [np.loadtxt("data/Vtrace-gsd%04.0f.txt"%(g*1e4)) for g in (0.21,0.23,0.245,0.263,0.28,0.315)]
time=np.arange(0,4,0.00005)
                                          
#%%
fig=plt.figure(1,figsize=(8,8))
#fig=plt.figure(1,figsize=(8,6))
fig.clf()
#
#ax1=plt.subplot2grid((5,3),(3,0))
#ax2=plt.subplot2grid((5,3),(3,1))
#ax3=plt.subplot2grid((5,3),(3,2))
#ax4=plt.subplot2grid((5,3),(4,0))
#ax5=plt.subplot2grid((5,3),(4,1))
#ax6=plt.subplot2grid((5,3),(4,2))
#ax0=plt.subplot2grid((5,3),(0,0),colspan=3,rowspan=3)
#
ax1=plt.subplot2grid((7,3),(3,0),colspan=1,rowspan=2)
ax2=plt.subplot2grid((7,3),(3,1),colspan=1,rowspan=2)
ax3=plt.subplot2grid((7,3),(3,2),colspan=1,rowspan=2)
ax4=plt.subplot2grid((7,3),(5,0),colspan=1,rowspan=2)
ax5=plt.subplot2grid((7,3),(5,1),colspan=1,rowspan=2)
ax6=plt.subplot2grid((7,3),(5,2),colspan=1,rowspan=2)
ax0=plt.subplot2grid((7,3),(0,0),colspan=3,rowspan=3)
    
im0=ax0.scatter(isis[:,0],isis[:,1],s=2,c=isis[:,2],marker='.',
           edgecolors='none',cmap='jet',vmin=0,vmax=1.2,)  
ax0.set_yscale('log')
ax0.set_xlim(0.18,0.33)
ax0.set_ylim((5,5000))
ax0.set_xlabel('$\mathsf{g_{sd} (mS/cm^2)}$',fontsize='large')
ax0.set_ylabel('ISI (ms)')

pos=plt.axes([0.12,0.65,0.014,0.14])
cb=plt.colorbar(im0,cax=pos,ticks=(0,0.6,1.2),extend='max')
cb.set_label('Lyapunov Exponent',fontsize='x-small')
cb.ax.tick_params(labelsize='x-small')

for label,pos in zip(range(1,7),(0.21,0.23,0.245,0.265,0.28,0.32)):
    ax0.text(pos,260,'(%i)'%label,va='bottom',ha='center',fontsize='small')

    
i=1

for trace,ax in zip(Vtraces,(ax1,ax2,ax3,ax4,ax5,ax6)):
    ax.plot(time,trace,'k',lw=1)
    ax.axis('off')
    if ax in (ax1,ax2,ax3):
        ax.set_xlim((1,2.5))
    else:
        ax.set_xlim((1,2))
    ax.plot((1,1.25),(-83,-83),'k',lw=2, clip_on = False)        
    ax.set_ylim((-80,10))
    ax.text(1,5,"(%d)"%i,va='bottom',ha='right')
    i+=1
    
plt.figtext(0.02,0.99,'A',fontsize='xx-large',va='top',ha='center')
plt.figtext(0.02,0.6,'B',fontsize='xx-large',va='top',ha='center')    
    
    
plt.tight_layout()
