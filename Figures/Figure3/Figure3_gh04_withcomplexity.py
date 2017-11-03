# -*- coding: utf-8 -*-
"""
Firing rate statistics
FIXED DISTRIBUTION OF gsd, gh

"""
import numpy as np
#import contextlib
import matplotlib.pyplot as plt
from matplotlib import colors
#import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.transforms as mtransforms
import spikesclassify_function as sc
import matplotlib as mpl
from matplotlib import rc, cm
from matplotlib.colors import Normalize

#loading the data
Table0 = np.loadtxt("finalTable/finalTable_gsdgsrzoomT36.txt",skiprows=1)
Table1 = np.loadtxt("finalTable/finalTable_ghgsd36.txt",skiprows=1)
Table2 = np.loadtxt("finalTable/finalTable_ghgsr36.txt",skiprows=1)#T
#loading complexity data
TableC0 = np.loadtxt("finalTable/gsdgsrT36ZoomComplexfinalTable_v2.txt",skiprows=1)
TableC1 = np.loadtxt("finalTable/ghgsdT36ComplexfinalTable_v2.txt",skiprows=1)
TableC2 = np.loadtxt("finalTable/ghgsrT36ComplexfinalTable_v2.txt",skiprows=1)#T
#loading the MLE data
Table3= np.loadtxt("finalTable/mlefinalTable_gsdgsr36.txt",skiprows=1)
Table4= np.loadtxt("finalTable/mlefinalTable_gsdgh36.txt",skiprows=1)
Table5= np.loadtxt("finalTable/mlefinalTable_gsrgh36.txt",skiprows=1)

#Fonts!
plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

#changing the xticks and  yticks fontsize for all sunplots
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('font',size=20)


# setting the fontsize of the plots
fsize = 18
ffsize= 28
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


lyap_min = 0
lyap_max = 1.6
cmap=plt.get_cmap('jet')
cmap.set_under('k',1)
norm=colors.Normalize(vmin=0,vmax=0.006)

aux =sc.lyap_filter(Table0)
aux1 =sc.lyap_filter(Table1)
aux2 =sc.lyap_filter(Table2)
auxC_0 =sc.Complex_filter(TableC0)
auxC_1 =sc.Complex_filter(TableC1)
auxC_2 =sc.Complex_filter(TableC2)
mle1=sc.mle_filter(Table3)
mle2=sc.mle_filter(Table4)
mle3=sc.mle_filter(Table5)

va,vb,vc=sc.spikes_classfy(Table0)
v1a,v1b,v1c=sc.spikes_classfy(Table1)
v2a,v2b,v2c=sc.spikes_classfy(Table2)


fig,axs=plt.subplots(nrows=4, ncols=3,figsize=(13,12))
plt.clf()
#axs = axs.ravel()


#assign the subplots to corresponding the variables such as ax1,ax2,ax3,...
#ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9=[plt.subplot(3,3,i) for i in range(1,10)]
ax1,ax5,ax9,ax2,ax6,ax10,ax3,ax7,ax11,ax4,ax8,ax12=[plt.subplot(4,3,i) for i in range(1,13)]

#data=[aux[0],aux[1],aux[2]]
#data1=[aux1[0],aux1[1],aux1[2]]
#data2=[aux2[0],aux2[1],aux2[2]]
axes=[ax1,ax2,ax3,ax4]
axes1=[ax5,ax6,ax7,ax8]
axes2=[ax9,ax10,ax11,ax12]
axsum=[plt.subplot(4,3,i) for i in range(1,13)]
#%%
#figure subplots
# plot of gsr against gsd
im1=ax1.imshow(aux[0],origin='lower', interpolation='none',vmin =0, vmax = 20, extent=aux[3],aspect='auto')
#plotting the MLE
im2=ax2.imshow(mle1[0],origin='lower',extent=mle1[2],norm=norm,aspect='auto',interpolation='none',cmap=cmap)
#im2=ax2.imshow(aux[1],origin='lower', interpolation='none',vmin =lyap_min, vmax = lyap_max, extent=aux[3],aspect='auto')
im3=ax3.imshow(auxC_0[0],origin='lower', interpolation='none',vmin =0, vmax = 1, extent=auxC_0[1],aspect='auto')
#plots of spikesclassfy
cmaps = mpl.colors.ListedColormap(['black','blue','palegreen','green','darkorange','brown','darkred'])
cmaps1 = mpl.colors.ListedColormap(['brown'])
im4a= ax4.imshow(vb, origin = 'lower', vmin=1, vmax=10, extent=aux[3], interpolation='nearest', cmap=cm.Oranges, norm = MidpointNormalize(midpoint=2), aspect='auto')#, alpha = 0.5)
#pb= ax.imshow(v1b, origin = 'lower', vmin=1, vmax=10, extent=[0,0.5,0, 0.4], interpolation='nearest', cmap=cm.Oranges, norm = LogNorm(vmin=0.01, vmax=10),aspect='auto')#, alpha = 0.5)
im4b = ax4.imshow(va, origin = 'lower', vmin=0, vmax=6, extent=aux[3], interpolation='nearest',cmap = cmaps,  aspect='auto')#, alpha = 0.5)

im4c = ax4.imshow(vc, origin = 'lower',  extent=aux[3], interpolation='nearest', cmap=cmaps1, aspect='auto')#,alpha = 0.5)

#%%
# plots  of second column,namely,# plot of gsd against gh
im5=ax5.imshow(aux1[0],origin='lower', interpolation='none', extent=aux1[3],aspect='auto')
#plots of the mle
im6=ax6.imshow(mle2[1],origin='lower',extent=mle2[3],norm=norm,aspect='auto',interpolation='none',cmap=cmap)
#plots of complecity
#im6=ax6.imshow(aux1[1],origin='lower', interpolation='none',vmin =lyap_min, vmax = lyap_max, extent=aux1[3],aspect='auto')
im7=ax7.imshow(auxC_1[0],origin='lower', interpolation='none',vmin =0, vmax = 1, extent=auxC_1[1],aspect='auto')
#plots of spikesclassfy
cmaps = mpl.colors.ListedColormap(['black','blue','palegreen','green','darkorange','brown','darkred'])
im8a= ax8.imshow(v1b, origin = 'lower', vmin=1, vmax=10, extent=aux1[3], interpolation='nearest', cmap=cm.Oranges, norm = MidpointNormalize(midpoint=2), aspect='auto')#, alpha = 0.5)
#pb= ax.imshow(v1b, origin = 'lower', vmin=1, vmax=10, extent=[0,0.5,0, 0.4], interpolation='nearest', cmap=cm.Oranges, norm = LogNorm(vmin=0.01, vmax=10),aspect='auto')#, alpha = 0.5)
im8b = ax8.imshow(v1a, origin = 'lower', vmin=0, vmax=6, extent=aux1[3], interpolation='nearest',cmap = cmaps,  aspect='auto')#, alpha = 0.5)

ax6.plot((0.17,0.33),(0.2,0.2),'w--',lw=3)
ax8.plot((0.17,0.33),(0.2,0.2),'w--',lw=3)


#im8c = ax8.imshow(v1c, origin = 'lower',  extent=aux1[3], interpolation='nearest', cmap=cmaps1, aspect='auto')#,alpha = 0.5)

#%%
# plot of gsr against gh
im9=ax9.imshow(aux2[0],origin='lower', interpolation='none', extent=aux2[3],aspect='auto')
im10=ax10.imshow(mle3[1],origin='lower',extent=mle3[3],norm=norm,aspect='auto',interpolation='none',cmap=cmap)
#im10=ax10.imshow(aux2[1],origin='lower', interpolation='none',vmin =lyap_min, vmax = lyap_max, extent=aux2[3],aspect='auto')

#plot of  complexcity
im11=ax11.imshow(auxC_2[0],origin='lower', interpolation='none',vmin =0, vmax = 1, extent=auxC_2[1],aspect='auto')

#plots of spikesclassfy
cmaps =mpl.colors.ListedColormap(['black','blue','palegreen','green','darkorange','brown','darkred'])
im12a= ax12.imshow(v2b, origin = 'lower', vmin=1, vmax=10, extent=aux2[3], interpolation='nearest', cmap=cm.Oranges, norm = MidpointNormalize(midpoint=2), aspect='auto')#, alpha = 0.5)
#pb= ax.imshow(v1b, origin = 'lower', vmin=1, vmax=10, extent=[0,0.5,0, 0.4], interpolation='nearest', cmap=cm.Oranges, norm = LogNorm(vmin=0.01, vmax=10),aspect='auto')#, alpha = 0.5)
im12b = ax12.imshow(v2a, origin = 'lower', vmin=0, vmax=6, extent=aux2[3], interpolation='nearest',cmap = cmaps,  aspect='auto')#, alpha = 0.5)
#im12c = ax12.imshow(v2c, origin = 'lower',  extent=aux2[3], interpolation='nearest', cmap=cmaps1, aspect='auto')#,alpha = 0.5)
plt.autoscale(False)
#%%
#setting labels and tickes of first colomn,namely,plot of gsr against gsd
#ax1.set_ylabel(r'$\mathsf{ g_{sd}}$', fontsize = ffsize)
for ax in axes:
#    ax.set_xlabel(r'$\mathsf{ g_{sr}}$', fontsize = ffsize)
    ax.set_ylabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$',fontsize='large')
    ax.set_xlim([0.2,0.345])
    ax.set_ylim([0.1,0.4])
    ax.set_yticks(np.linspace(0.1,0.4,4))
    ax.set_xticks(np.linspace(0.2,0.345,4))
    ax.set_xticklabels(['%0.2f'%(0.2+i*0.05) for i in range(4)])
    ax.tick_params(axis='x', pad=8)
#setting the xticksets
for ax in (ax1,ax2,ax3):
    ax.set_xticks([])
ax4.set_xlabel(r'$\mathsf{ g_{sr} \; (mS/cm^2)}$',fontsize='large')


#setting labels and tickes of second column,namely,plot of gsd against gh

for ax in axes1:
#    ax.set_xlabel(r'$\mathsf{ g_{sd}}$', fontsize = ffsize)
    ax.set_ylabel(r'$\mathsf{ g_{h} \; (mS/cm^2)}$',fontsize='large')
    ax.set_xlim([0.17,0.33])
    ax.set_ylim([0,0.595])
    ax.set_xticks(np.linspace(0.17,0.33,5))
    ax.set_yticks(np.linspace(0.0,0.595,4))
    ax.set_yticklabels(['%0.1f'%(i*0.2) for i in range(4)])
    ax.tick_params(axis='x', pad=8)
#setting the xticksets
for ax in (ax5,ax6,ax7):
    ax.set_xticks([])
ax8.set_xlabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)} $',fontsize='large')

#etting labels and tickes of thirds column,namely,plot of gsr against gh.
for ax in axes2:
#    ax.set_xlabel(r'$\mathsf{ g_{sr}}$', fontsize = ffsize)
    ax.set_ylabel(r'$\mathsf{ g_{h} \; (mS/cm^2)}$',fontsize='large')
    ax.set_xlim([0.2,0.32])
    ax.set_ylim([0,0.6])
    ax.set_xticks(np.linspace(0.2,0.32,4))
    ax.set_yticks(np.linspace(0.0,0.6,4))
    ax.tick_params(axis='x', pad=8)
#setting the xticksets
for ax in (ax9,ax10,ax11):
    ax.set_xticks([])
ax12.set_xlabel(r'$\mathsf{ g_{sr} \; (mS/cm^2)}$',fontsize='large')
#setting labels and tickes of firing patterns



# setting the text for the subplots
index_list=['A','B','C','D']
axindex=[ax1,ax2,ax3,ax4]
for ax,plotor in zip(axindex,index_list):
    ax.text(0.15,0.38,'%s'%plotor,fontsize =24)


for ax in axsum:
    divider1 = make_axes_locatable(ax)
    cax1= divider1.append_axes("right","18%", pad="3%")
    fig.delaxes(cax1)

#colorbar of  Firing rate
##To create an axes,n axes at position rect [left, bottom, width, height]
cax = fig.add_axes([0.94, 0.77, 0.01, 0.21])
cbar=fig.colorbar(im1, extend='max',cax=cax)
cbar.set_label(u"Firing Rate ",fontsize=15)
cbar.set_ticks(np.linspace(0,20,5))
cbar.ax.tick_params(labelsize=12)

#colorbar of MLE
cax1 = fig.add_axes([0.935, 0.54, 0.01, 0.21])
#plt.colorbar(gci,extend='both', ticks=np.linspace(0,0.005,6))
cbar1=fig.colorbar(im2,extend='both',cax=cax1)
cbar1.set_label(u"MLE",fontsize=15)
cbar1.set_ticks(np.linspace(0,0.006,4))
#change the appearance of ticks anf tick labbel
cbar1.ax.tick_params(labelsize=12)
#cax1 = fig.add_axes([0.94, 0.54, 0.01, 0.21])
#cbar1=fig.colorbar(im2, extend='max',cax=cax1)
#cbar1.set_label('Lyapunov  Exponent ',fontsize=18)
#cbar1.set_ticks(np.linspace(0,1.6,5))
#change the appearance of ticks anf tick labbel
#cbar1.ax.tick_params(labelsize=14)

#colorbar of comlexity
cax2 = fig.add_axes([0.94, 0.305, 0.01, 0.21])
cbar2=fig.colorbar(im3, extend='max', cax=cax2)
cbar2.set_label(u"Lempel-Ziv Complexity",fontsize=15)
cbar2.set_ticks(np.linspace(0,1,6))
cbar2.ax.tick_params(labelsize=12)
#colorbar of firing rate
cax3 = fig.add_axes([0.94, 0.07, 0.012, 0.21])
cbar3 = fig.colorbar(im8b, ticks=[ 0.4, 1.3, 2.1, 3, 3.8, 4.7, 5.5],cax=cax3)
cbar3.ax.set_yticklabels(['0', '1', '2', '3', '4', '5', '6'])
cbar3.set_label(u"Firing  pattern", fontsize=15,labelpad=10)
cbar3.ax.tick_params(labelsize=12)



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

plt.subplots_adjust(bottom=0.07,left=0.09,wspace = 0.1,hspace = 0.1,right=0.98, top=0.98)
plt.savefig('Figure3.png',dpi=300)
plt.savefig('Figure3.eps',dpi=300)
