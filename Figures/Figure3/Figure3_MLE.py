# -*- coding: utf-8 -*-
"""
Firing rate statistics
FIXED DISTRIBUTION OF gsd, gh

"""
import numpy as np
#import contextlib
import matplotlib.pyplot as plt
#from matplotlib import colors
#import matplotlib.gridspec as gridspec
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.transforms as mtransforms
#import spikesclassify_function as sc
#import matplotlib as mpl
#from matplotlib import rc, cm
#from matplotlib.colors import Normalize

#loading the data
#Table0 = np.loadtxt("finalTable/finalTable_gsdgsrzoomT36.txt",skiprows=1)
#Table1 = np.loadtxt("finalTable/finalTable_ghgsd36.txt",skiprows=1)

#loading the MLE data
Table0= np.loadtxt("finalTable/mlefinalTable_gsdgsr36.txt",skiprows=1)
Table1= np.loadtxt("finalTable/mlefinalTable_gsdgh36.txt",skiprows=1)
Table1=Table1[:,[1,0,2]]
Table1=Table1[np.lexsort((Table1[:,1],Table1[:,0]))]

#Fonts!
plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

#changing the xticks and  yticks fontsize for all sunplots
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 
plt.rc('font',size=18)


# setting the fontsize of the plots
fsize = 18
ffsize= 28
im=[]

cmap=plt.get_cmap('jet')
cmap.set_under('k',1.)
#%%

ghmin1, ghmax1, gsdmin1, gsdmax1= 0.35, 0.4, 0.197, 0.207 #range chaos 1 gh/gsd T36
ghmin2, ghmax2, gsdmin2, gsdmax2= 0.1, 0.15, 0.197, 0.207 #range nonchaos 1 gh/gsd T36

def plotLE(ax,dataTable,index,axlabels):
    y_range=np.array([np.min(dataTable[:,0]),np.max(dataTable[:,0])])#*1000
    x_range=np.array([np.min(dataTable[:,1]),np.max(dataTable[:,1])])#*1000
    yN=len(np.unique(dataTable[:,0]))
    xN=len(np.unique(dataTable[:,1]))

    ratio=np.diff(x_range)/np.diff(y_range)
    print ratio

    Lyap=np.reshape(dataTable[:,-1],(yN,xN))
    
    im=ax.imshow(Lyap,origin='lower',extent=(x_range[0],x_range[1],y_range[0],y_range[1]),
           aspect=ratio,interpolation='none',cmap=cmap,vmin=0,vmax=0.005)
    ax.set_xlabel(axlabels[0],fontsize='large')
    ax.set_ylabel(axlabels[1],fontsize='large')
           
    return im
           
plt.figure(1,figsize=(12,5))
plt.clf()

ax1=plt.subplot(121)
im1=plotLE(ax1,Table1,0,(u"$g_{sd} [mS/cm^2]$",u"$g_{h} [mS/cm^2]$"))
xL=ax1.get_xlim()
yL=ax1.get_ylim()
ax1.plot((gsdmin1,gsdmin1,gsdmax1,gsdmax1,gsdmin1),(ghmin1,ghmax1,ghmax1,ghmin1,ghmin1),'r-',lw=3)
ax1.plot((gsdmin2,gsdmin2,gsdmax2,gsdmax2,gsdmin2),(ghmin2,ghmax2,ghmax2,ghmin2,ghmin2),'g-',lw=3)
ax1.set_xlim(xL)
ax1.set_ylim(yL)


ax2=plt.subplot(122)
im2=plotLE(ax2,Table0,0,('$g_{sr} [mS/cm^2]$','$g_{sd} [mS/cm^2]$'))

cbar=plt.colorbar(im2,extend='max')#,ticks=(0,0.3,0.6,0.9,1.2,1.5))
cbar.set_label('Lyapunov Exponent')

plt.figtext(0.01,0.9,'A',fontsize='x-large')
plt.figtext(0.48,0.9,'B',fontsize='x-large')

plt.tight_layout()