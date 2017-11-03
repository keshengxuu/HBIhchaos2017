# -*- coding: utf-8 -*-
"""
Firing rate statistics
FIXED DISTRIBUTION OF gsd, gh

"""
import numpy as np
import contextlib
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as mtransforms
import spikesclassify_function as sc
import matplotlib as mpl
from matplotlib import rc, cm
from matplotlib.colors import Normalize

#loading the data
Table0 = np.loadtxt("Data/finalTable_gsdgsrZoomT36gh00.txt",skiprows=1)
Table1 = np.loadtxt("Data/mlegsdgsrfinalTable_HB5.txt",skiprows=1)
#loading complexity data
TableC0 = np.loadtxt("Data/gsdgsrT36gh00_ComplexfinalTable.txt",skiprows=1)

#Fonts!
plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

#changing the xticks and  yticks fontsize for all sunplots
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.rc('font',size=15)


# setting the fontsize of the plots
fsize = 18
ffsize= 18
im=[]


alf = 0.8
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

fig,axs=plt.subplots(nrows=4, ncols=3,figsize=(8,6))
plt.clf()
#axs = axs.ravel()
lyap_min = 0
lyap_max = 1.6
cmap=plt.get_cmap('jet')
cmap.set_under('k',1)
norm=colors.Normalize(vmin=0,vmax=0.006)
aux =sc.lyap_filter(Table0)
auxC_0 =sc.Complex_filter(TableC0)
va,vb,vc=sc.spikes_classfy(Table0)
mle1=sc.mle_filter(Table1)


#assign the subplots to corresponding the variables such as ax1,ax2,ax3,...
#ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9=[plt.subplot(3,3,i) for i in range(1,10)]
ax1,ax2,ax3,ax4=[plt.subplot(2,2,i) for i in range(1,5)]

#data=[aux[0],aux[1],aux[2]]
#data1=[aux1[0],aux1[1],aux1[2]]
#data2=[aux2[0],aux2[1],aux2[2]]
axes=[ax1,ax2,ax3,ax4]
axsum=[plt.subplot(2,2,i) for i in range(1,5)]
#%%
#figure subplots
# plot of gsr against gsd
im1=ax1.imshow(aux[0],origin='lower', interpolation='none',vmin =0, vmax = 20, extent=aux[3],aspect='auto')
im2=ax2.imshow(mle1[0],origin='lower',extent=mle1[2],norm=norm,aspect='auto',interpolation='none',cmap=cmap)


im3=ax3.imshow(auxC_0[0],origin='lower', interpolation='none',vmin =0, vmax = 1, extent=auxC_0[1],aspect='auto')
#plots of spikesclassfy
cmaps = mpl.colors.ListedColormap(['black','blue','palegreen','green','darkorange','brown','darkred'])
cmaps1 = mpl.colors.ListedColormap(['brown'])

im4a= ax4.imshow(vb, origin = 'lower', vmin=1, vmax=10, extent=aux[3], interpolation='nearest', cmap=cm.Oranges, norm = MidpointNormalize(midpoint=2), aspect='auto')#, alpha = 0.5)
#pb= ax.imshow(v1b, origin = 'lower', vmin=1, vmax=10, extent=[0,0.5,0, 0.4], interpolation='nearest', cmap=cm.Oranges, norm = LogNorm(vmin=0.01, vmax=10),aspect='auto')#, alpha = 0.5)
im4b = ax4.imshow(va, origin = 'lower', vmin=0, vmax=6, extent=aux[3], interpolation='nearest',cmap = cmaps,  aspect='auto')#, alpha = 0.5)

im4c = ax4.imshow(vc, origin = 'lower',  extent=aux[3], interpolation='nearest', cmap=cmaps1, aspect='auto')#,alpha = 0.5)

plt.autoscale(False)
#%%
#setting labels and tickes of first colomn,namely,plot of gsr against gsd
#ax1.set_ylabel(r'$\mathsf{ g_{sd}}$', fontsize = ffsize)
for ax in axes:
#    ax.set_xlabel(u'$g_{sr}[mS/mc^2]$')
#    ax.set_ylabel(r'$\mathsf{ g_{sd} \; (mS/mc^2)}$')
    ax.set_xlim([0.2,0.345])
    ax.set_ylim([0.1,0.4])
    ax.set_yticks(np.linspace(0.1,0.4,4))
    ax.set_xticks(np.linspace(0.2,0.345,4))
    ax.set_xticklabels(['%0.2f'%(0.2+i*0.05) for i in range(4)])
    ax.tick_params(axis='x', pad=8)
#setting the xticksets
ax1.set_ylabel(r'$\mathsf{g_{sd} \; (mS/cm^2)}$')
ax3.set_ylabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$')
ax3.set_xlabel(r'$\mathsf{g_{sr} \; (mS/cm^2)}$')
ax4.set_xlabel(r'$\mathsf{g_{sr}\; (mS/cm^2)}$')


# setting the text for the subplots
index_list=['A','B','C','D']
axindex=[ax1,ax2,ax3,ax4]
for ax,plotor in zip(axindex,index_list):
    ax.text(0.16,0.41,'%s'%plotor,fontsize =18)


for ax in axsum:
    divider1 = make_axes_locatable(ax)
    cax1= divider1.append_axes("right","35%", pad="3%")
    fig.delaxes(cax1)

##To create an axes,n axes at position rect [left, bottom, width, height]
cax = fig.add_axes([0.42, 0.56, 0.015, 0.38])
cbar=fig.colorbar(im1, extend='max',cax=cax)
cbar.set_label(r'$\mathsf{Firing Rate}$ ',fontsize=15)
cbar.set_ticks(np.linspace(0,20,5))
cbar.ax.tick_params(labelsize=10)
#
cax1 = fig.add_axes([0.90, 0.56, 0.015, 0.38])
cbar1=fig.colorbar(im2,extend='both',cax=cax1)
cbar1.set_label(r'$\mathsf{MLE}$',fontsize=15)
cbar1.set_ticks(np.linspace(0,0.006,4))
#change the appearance of ticks anf tick labbel
cbar1.ax.tick_params(labelsize=10)

#colorbar of comlexity
cax2 = fig.add_axes([0.42, 0.10, 0.015, 0.38])
cbar2=fig.colorbar(im3, extend='max', cax=cax2)
cbar2.set_label(r'$\mathsf{Lempel-Ziv Complexity}$',fontsize=15)
cbar2.set_ticks(np.linspace(0,1,6))
cbar2.ax.tick_params(labelsize=12)
#colorbar of firing pattern
cax4 = fig.add_axes([0.9, 0.10, 0.015, 0.38])
cbar4 = fig.colorbar(im4b, ticks=[ 0.4, 1.3, 2.1, 3, 3.8, 4.7, 5.5],cax=cax4)
cbar4.ax.set_yticklabels(['0', '1', '2', '3', '4', '5', '6'])
cbar4.set_label( r'$\mathsf{Firing  Pattern}$', fontsize=15,labelpad=10)
cbar4.ax.tick_params(labelsize=10)



##To create an axes,n axes at position rect [left, bottom, width, height]
#cbaxes = fig.add_axes([0.17, 0.26, 0.07, 0.01]) 
#cb = plt.colorbar(pb, ) 
#cbb = plt.colorbar(im4a,orientation="horizontal",ticks=[1,3,5,7,9], cax = cbaxes )#,ticks=[1,3,5,7,9], shrink = 1.1)

#subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
#left  = 0.125  # the left side of the subplots of the figure
#right = 0.9    # the right side of the subplots of the figure
#bottom = 0.1   # the bottom of the subplots of the figure
#top = 0.9      # the top of the subplots of the figure
#wspace = 0.2   # the amount of width reserved for blank space between subplots
#hspace = 0.2   # the amount of height reserved for white space between subplots

plt.subplots_adjust(bottom=0.1,left=0.1,wspace = 0.15,hspace = 0.2,right=1, top=0.94)
plt.savefig('Figure4_withoutgh.png')
