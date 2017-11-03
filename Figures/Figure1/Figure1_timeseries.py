# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 15:03:37 2016
keshengxuu@gmail.com
@author: keshengxu
"""
import numpy as np
from scipy.stats import norm
from scipy import integrate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os.path
from matplotlib import colors
import matplotlib.gridspec as gridspec


#Fonts!
plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)

font= { 'family' : 'sans-serif',
#        'color' : 'darkred',
	'weight' : 'normal',
	'size'   : 12,
	}
#
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward',0))  # outward by 10 points
#            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


#loading the data
time = np.loadtxt('Data/time.txt')
v1,v2,v3,v4,v5=[np.loadtxt('Data/var_t%g.txt'%i) for i in range(1,6)]
ISI1,ISI2,ISI3,ISI4,ISI5=[np.loadtxt('Data/ISI%g.txt'%i)for i in range(1,6)]
spikes6=np.loadtxt('Data/spikes-T36.300.txt')
ISI6 = np.diff(spikes6)

#ISIs=np.diff(spkf)
binss = 40#20
#linf = 20  Antiguo
linf = 0
lsup = 800

fig,axs=plt.subplots(figsize=(8,6))
plt.clf()

ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15=[plt.subplot(5,3,i)for i in range(1,16)]

axes_v=[ax1,ax4,ax7,ax10,ax13]
axes_ISI=[ax2,ax5,ax8,ax11,ax14]
aexs_cou=[ax3,ax6,ax9,ax12,ax15]
color_set=['chocolate','darkorange','orange','green','palegreen']
#saving the Membrane potential
Vdata=[v1,v2,v3,v4,v5]
ISIdata=[ISI1,ISI2,ISI3,ISI4,ISI5]
ISIdata1=[ISI1,ISI2,ISI3,ISI4,ISI6]
simu_time=[np.cumsum(ISI1),np.cumsum(ISI2),np.cumsum(ISI3),np.cumsum(ISI4),np.cumsum(ISI5)]

# set the order of the plots
plot_order=['A','B','C','D','E']
#figure out the membrane action potential
for ax,V,cl,plotor in zip(axes_v,Vdata,color_set,plot_order):
    ax.plot(time,V,color='black', linewidth = 1)
    ax.set_ylabel(u'Voltage (mV)',fontsize=11)
    ax.set_ylim([-90,30])
    ax.set_xlim([30000,34000])
    adjust_spines(ax, ['left'])
    ax.set_yticks(np.linspace(-90,30,3))
    ax.text(28500,32,'%s'%plotor,fontsize=12)
    ax.axes.tick_params(direction="out")

#adjust_spines(ax13, ['left', 'bottom'])
#for ax in axes_v[0:5]:
#    adjust_spines(ax, ['left'])
#    ax.tick_params(labelsize=10)
#adjust_spines(ax13, ['left','bottom'])
#ax13.set_xticklabels(['30','31','32','33','34','35'],fontsize = 8)
#ax13.spines['bottom'].set_position(('outward',5))
##Change the appearance of ticks and tick labels.
#ax13.tick_params(axis='x',direction='out',labelsize=10)

#figure out the ISI
for ax,st,ISI in zip(axes_ISI,simu_time,ISIdata):
    ax.plot(st/1000,ISI, 'bo', ms=1)
    adjust_spines(ax, ['left','bottom'])
    ax.set_ylim([10**1,10**3])
    ax.set_xlim([0,150])
    ax.set_yscale('log')
    ax.set_ylabel('ISI (ms)',fontdict=font)
    ax.tick_params(axis='y',direction='out', length=4, width=1)


for ax in axes_ISI[0:4]:
    ax.axes.tick_params(axis='x',direction='out', length=0, width=2)
    ax.set_xticks([])
ax14.axes.tick_params(axis='x',direction='out')
ax14.set_xticks(np.linspace(0,150,4))
ax14.set_xlabel('Time (s)',fontdict=font)


# list of temperaure
tem_text=['20','24.76','26','33','36.3']
# Histgram of ISI
for ax,ISI,temt in zip(aexs_cou,ISIdata1,tem_text):
    hist, bins = np.histogram(ISI, bins=binss,range = (linf, lsup))
    widths = np.diff(bins)
    ax.bar(bins[:-1],np.sqrt(hist), widths)
#    ax.bar(bins[:-1],hist, widths)
    adjust_spines(ax, ['left','bottom'])
    ax.set_ylabel('Event Count \n (sq.root)',fontdict=font,fontsize=10)
#    ax.set_xlabel('ISI (ms)',fontdict=font)
    ax.set_ylim([0,32])
    ax.set_yticks(np.linspace(0,30,3))
    ax.tick_params(axis='y',direction='out', length=4, width=1)
    ax.text(500,25,r'$\mathsf{ %s^{\circ} C}$'%temt,fontdict=font,fontsize=14)

#set the label and ticks
for ax in aexs_cou[0:4]:
    ax.axes.tick_params(axis='x',direction='out', length=0, width=2)
    ax.set_xticks([])
ax15.axes.tick_params(axis='x',direction='out')
ax15.set_xticks(np.linspace(0,800,5))
ax15.set_xlabel('ISI (ms)',fontdict=font)


#define a scalebar under the menberane action potentional
def figura_scalebar(ax):
    ax.set_ylim([-90,30])
    ax.set_xlim([30000,34000])
    ax.hlines(-85,32000,33000,color='black',linewidth=2)
    ax.text(32300,-120,r'$\mathsf{ 1\;s}$',fontsize=12)
#    ax.axis('off')
#
#ax0 = fig.add_axes([0.1, -0.1, 0.01, 0.271])
##figura_scalebar(ax0)
ax0 = fig.add_subplot(5,3,13)
figura_scalebar(ax0)

plt.subplots_adjust(bottom=0.08,left=0.1,wspace = 0.4,hspace = 0.18,right=0.98, top=0.97)
plt.draw()
plt.savefig('Figure1_timeseries.png',dpi=600)
