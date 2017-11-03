# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 13:40:38 2016

@author: keshengxu
"""
import numpy as np
import contextlib
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as mtransforms
import csv

#Fonts!
plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

#change the font size of xtick and  ytick
plt.rc('xtick',labelsize=12)
plt.rc('ytick', labelsize=12)


ffsize= 18
fffsize=12
lyap_min = 0
lyap_max = 1.6


#loading the data
Table0 = np.loadtxt("Data/finalTable/finalTable_tsrth36.txt",skiprows=1)
Table1 = np.loadtxt("Data/finalTable/finalTable_tsdthT36.txt",skiprows=1)
ISItable=np.loadtxt("Data/ISIs_th.txt.gz")
MLE=np.loadtxt('Data/mlethfinalTable_HB5.txt')
MLETable=np.loadtxt('Data/mlethtsdfinalTable_HB5.txt',skiprows=1)
MLETable1=np.loadtxt('Data/mlethtsrfinalTable_HB5.txt',skiprows=1)


def lyap_filter(Table):

    Lyap=Table[:,-1]
    Freq=Table[:,2]

    ghrange=np.array([np.min(Table[:,0]),np.max(Table[:,0])])
    gsdrange=np.array([np.min(Table[:,1]),np.max(Table[:,1])])
    ghN=len(np.unique(Table[:,0]))
    gsdN=len(np.unique(Table[:,1]))

    ratio=np.diff(gsdrange)/np.diff(ghrange)

    Lyap2=np.reshape(Lyap,(ghN,gsdN))
    Freq2=np.reshape(Freq,(ghN,gsdN))

    minFR=2
    maxFR=4.5

    maxFR=6


    Mask1=0.6 + 0.4* ((Freq2>minFR)*(Freq2<maxFR))

    cmap=plt.get_cmap('jet')
    cmap.set_under('k',1.)

    FreqN=np.minimum(Freq2,20*np.ones_like(Freq2))/20.  #Maxima freq será 20.
    FreqN[FreqN==0]=-10

    ImgFreq=cmap(FreqN)
    ImgFreq2=ImgFreq*Mask1[:,:,None]

#ImgFreq3=ImgFreq*nonchaosMask[:,:,None]

    ImgLyap=cmap(Lyap2)
    ImgLyap2=ImgLyap*Mask1[:,:,None]
    ImgLyap2[Lyap2==-10,3]=1
    extent = (gsdrange[0],gsdrange[1],ghrange[0],ghrange[1])
    output = [ImgFreq, ImgLyap, ImgLyap2, extent, ratio]
    return output


def mle_filter(Table):
    lya0=Table[:,2]
    ghrange=np.array([np.min(Table[:,1]),np.max(Table[:,1])])
    gsdrange=np.array([np.min(Table[:,0]),np.max(Table[:,0])])
    ghN=len(np.unique(Table[:,1]))
    gsdN=len(np.unique(Table[:,0]))
    ratio=np.diff(gsdrange)/np.diff(ghrange)
    lya00=np.reshape(lya0,(gsdN,ghN)).transpose()
    extent = (gsdrange[0],gsdrange[1],ghrange[0],ghrange[1])
    output =[lya00,extent]
    return output

# plot  the  figures
aux = lyap_filter(Table0)
aux1 = lyap_filter(Table1)
mle1=mle_filter(MLETable)
mle2=mle_filter(MLETable1)

cmap=plt.get_cmap('jet')
cmap.set_under('k',1)
norm=colors.Normalize(vmin=0.000,vmax=0.006)

#%%
fig = plt.figure(4, figsize=(10, 8))
fig.clf()
gs = gridspec.GridSpec(100, 100, wspace=0, hspace=0.1)


#The  First subplot
ax1 = plt.subplot(gs[:45,:92])

im1=ax1.scatter(ISItable[:,0],ISItable[:,1],s=1.6,linewidth=0,c=ISItable[:,2],
                         vmin=0 ,vmax=1.6, alpha =0.9)

plt.xlabel(r'$\mathsf{\tau_{h} \; (ms)}$', fontsize = ffsize)
plt.ylabel(r'$\mathsf{ISI \; (ms) }$',position=(85,0.4), fontsize = ffsize)
#plt.grid(True)
plt.yscale('log')
plt.xlim([94,220])
plt.ylim([10**1.8,10**5.2])
plt.xticks([100,130,160,190,220])
plt.yticks([10**2,10**3,10**4],[r'$\mathsf{10^2}$',r'$\mathsf{10^3}$',r'$\mathsf{10^4}$'],fontsize=14)

ax11 = ax1.twinx()
#s2 = np.sin(2*np.pi*t)
#ax11.plot(t, s2, 'r.')
#ax11.set_ylabel('sin', color='r')
plt.grid(True)
plt.plot(MLE[:,0],MLE[:,1],'b')
plt.xlim([94,220])
plt.ylim([-0.016,0.00801])
plt.ylabel(r'$\mathsf{MLE}$',position=(222,0.8),fontsize = 15,color='b')
plt.yticks([0,0.005],['0.000','0.005'])
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.tick_params(axis="y", labelcolor="b")#,pad=8)
#ax11.set_xticks([1000,1300,1600,1900,2200],['100','130','160','190','220'])
for tl in ax11.get_yticklabels():
    tl.set_color('b')

cax5 = fig.add_axes([0.84, 0.62, 0.01, 0.16])
cbar=fig.colorbar(im1, extend='max',cax=cax5,ticks=[ 0.0,0.8,1.6])
cbar.set_label(r'$\mathsf{Lyapunov \;Exponent}$',fontsize=10)
cbar.ax.set_yticklabels(['0.0','0.8','1.6'],fontsize =8)
#%%
#The second subplot
ax2 = plt.subplot(gs[54:,:40])
im2=ax2.imshow(mle1[0],origin='lower',extent=mle1[1],norm=norm,aspect='auto',interpolation='none',cmap=cmap)
ax2.set_xlabel(r'$\mathsf{\tau_{h} \; (ms)}$', fontsize = ffsize)
ax2.set_ylabel(r'$\mathsf{ \tau _{sd} \; (ms)}$', fontsize = ffsize)
ax2.set_ylim([5,40])
ax2.set_xticks(np.linspace(50,500,6))
ax2.set_yticks(np.linspace(5,40,6))
ax2.plot((50,500),(10,10),'w--',lw=3)



#divider = make_axes_locatable(plt.gca())
#cax = divider.append_axes("right", "5%", pad="3%")
#cbar=plt.colorbar(im1,extend='max', cax=cax, ticks=np.linspace(0,20,5))
#cbar.set_label(r'$\mathsf{Firing\; Rate }$',fontsize=fffsize)

#%%
#the third subplot
ax3 = plt.subplot(gs[54:,52:92])
im3 =ax3.imshow(mle2[0],origin='lower',extent=mle2[1],norm=norm,aspect='auto',interpolation='none',cmap=cmap)
plt.xlabel(r'$\mathsf{\tau_{h} \; (ms)}$', fontsize = ffsize)
plt.ylabel(r'$\mathsf{ \tau _{sr} \; (ms)}$', fontsize = ffsize)
ax3.set_xticks(np.linspace(50,500,6))
ax3.set_yticks(np.linspace(10,100,6))
ax3.plot((50,500),(35,35),'w--',lw=3)

ax1.text(-0.115,0.9, 'A',  transform=ax1.transAxes, fontsize = 20)
ax2.text(-0.26,0.9, 'B',  transform=ax2.transAxes, fontsize = 20)
ax3.text(-0.2,0.9, 'C',  transform=ax3.transAxes, fontsize = 20)


ax8.plot((0.17,0.33),(0.2,0.2),'w--',lw=3)




### setting the colorbar
#pltgca=[ax2,ax3]
#for plc in pltgca:
#    divider = make_axes_locatable(plc)
#    cax = divider.append_axes("right", "5%", pad="3%")
#    fig.delaxes(cax)
#
#divider3 = make_axes_locatable(ax1)
#cax3 = divider3.append_axes("right", "0.7%", pad="3%")
#fig.delaxes(cax3)
##
###To create an axes,n axes at position rect [left, bottom, width, height]
cax5 = fig.add_axes([0.915, 0.08, 0.01, 0.42])
#plt.colorbar(gci,extend='both', ticks=np.linspace(0,0.005,6))
cbar5=fig.colorbar(im2,extend='both',cax=cax5)
cbar5.set_label(r'$\mathsf{MLE}$',fontsize=15)
cbar5.set_ticks(np.linspace(0.00,0.006,4))

plt.subplots_adjust(left = 0.1,bottom=0.08, right=0.98, top=0.98, wspace=0.05, hspace=0.45)
plt.draw()
plt.savefig('Figure6_timeconstant.png',dpi=300)



