# -*- coding: utf-8 -*-
"""
Firing rate statistics
FIXED DISTRIBUTION OF gsd, gh

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors 

#Table=np.loadtxt("data/ghgsdT36_finalTable.txt",skiprows=1)

BifTable=np.loadtxt("data/bifurc2D_oscil_ghgsd_ToPlot.dat",delimiter=',')
#Table=np.loadtxt("bifurc_oscil_gsd/bifurc2Dghgsd_oscil_allinfo.dat")

PLpoints=[]

HB1=(BifTable[:,4]==2)
LP1=(BifTable[:,4]==4)
LP2=(BifTable[:,4]==33)+(BifTable[:,4]==34)
LP3=(BifTable[:,4]==48)+(BifTable[:,4]==49)

PDs=[(BifTable[:,4]==i)+(BifTable[:,4]==(i+1)) for i in (6,8,10,12,14,16,25,27,29,31,44,46)]



Table=np.loadtxt('data/mlefinalTable.txt')
lya0=Table[:,2]
ghrange=np.array([np.min(Table[:,1]),np.max(Table[:,1])])
gsdrange=np.array([np.min(Table[:,0]),np.max(Table[:,0])])

ghN=len(np.unique(Table[:,1]))
gsdN=len(np.unique(Table[:,0]))
ratio=np.diff(gsdrange)/np.diff(ghrange)
lya00=np.reshape(lya0,(gsdN,ghN)).transpose()
extent = (gsdrange[0],gsdrange[1],ghrange[0],ghrange[1])

cmap=plt.get_cmap('jet')
cmap.set_under('k',1)
norm=colors.Normalize(vmin=-0.00005,vmax=0.006)

plt.figure(1,figsize=(12,5))
plt.clf()

#plt.subplot(131)
#
#plt.imshow(ImgLyap,origin='lower',extent=(gsdrange[0],gsdrange[1],ghrange[0],ghrange[1]),
#           aspect=ratio,interpolation='none')
##plt.colorbar()           
#           
#plt.xlabel("$g_{sd}$",fontsize='x-large')
#plt.ylabel("$g_{h}$",fontsize='x-large')

plt.subplot2grid((1,9),(0,0),colspan=4)

plt.plot(BifTable[HB1,0],BifTable[HB1,1],'g-',lw=1)
for cond in PDs:
    plt.plot(BifTable[cond,0],BifTable[cond,1],'r-',lw=1)
for cond in LP1,LP2,LP3:
    plt.plot(BifTable[cond,0],BifTable[cond,1],'m-',lw=2)   

plt.gca().set_aspect(ratio)
   
plt.text(0.165,0.5,'HB',color='g',size='x-large')
plt.text(0.198,0.48,'PD',color='r',size='x-large')
plt.text(0.255,0.42,'PD',color='r',size='x-large')
plt.text(0.3,0.3,'LP',color='m',size='x-large')
plt.text(0.233,0.43,'LP',color='m',size='x-large')
    
plt.xlim((0.14,0.35))
plt.ylim((0,0.6))
plt.xlabel("$\mathsf{g_{sd}}$",fontsize='x-large')
plt.ylabel("$\mathsf{g_{h}}$",fontsize='x-large')


plt.subplot2grid((1,9),(0,4),colspan=5)
gci=plt.imshow(lya00,origin='lower',extent=extent,norm=norm,aspect='auto',interpolation='none',cmap=cmap)
#plt.imshow(ImgLyap,origin='lower',extent=(gsdrange[0],gsdrange[1],ghrange[0],ghrange[1]),
#          aspect=ratio,interpolation='none')
plt.plot(BifTable[HB1,0],BifTable[HB1,1],'g-',lw=1)
for cond in PDs:
    plt.plot(BifTable[cond,0],BifTable[cond,1],'r-',lw=1)
for cond in LP1,LP2,LP3:
    plt.plot(BifTable[cond,0],BifTable[cond,1],'m-',lw=2)   

plt.xlim([0.14,0.3495])
plt.ylim([0,0.5975])
plt.xticks([0.15,0.2,0.25,0.3,0.3495],['0.15','0.20','0.25','0.3','0.35'], fontsize = 12)
plt.yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.5975],['0.0','0.1','0.2','0.3','0.4','0.5','0.6'], fontsize = 12)
plt.xlabel("$\mathsf{g_{sd}}$",fontsize='x-large')
plt.ylabel("$\mathsf{g_{h}}$",fontsize='x-large')

#pos=plt.subplot2grid((1,9),(0,8))
cbar = plt.colorbar(gci,extend='both', ticks=np.linspace(0,0.006,4),pad = 0.01)
cbar.set_label(u'MLE',fontsize='large')

plt.figtext(0.01,0.92,'A',fontsize='xx-large')
plt.figtext(0.44,0.92,'B',fontsize='xx-large')

plt.tight_layout()
plt.savefig('Figure8.png')