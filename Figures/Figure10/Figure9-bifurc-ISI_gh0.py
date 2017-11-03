# -*- coding: utf-8 -*-
"""

"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker

#Fonts!
plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

bifTable=np.loadtxt('data/bifurc1D_oscil_gh0_ToPlot.dat')
ISIdata=np.loadtxt('data/gsdT36gh00_ISI-plotData.txt.gz')
Table=np.loadtxt("data/gsdT36gh00_finalTable.txt", delimiter='\t')

br1S=(bifTable[:,1]==1)*(bifTable[:,0]==1)
br1U=(bifTable[:,1]==1)*(bifTable[:,0]==2)

PerStBranch = (bifTable[:,1]==-2)*(bifTable[:,0]==3)
PerUnstBranch = (bifTable[:,1]==-2)*(bifTable[:,0]==4)

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

ISIdata=np.c_[ISIdata,np.zeros(len(ISIdata))]
for line in ISIdata:
    line[-1]=Table[Table[:,0]==line[0]/1e7,-1]

#%%

#plt.rc('text', usetex=False)
#plt.rc('font', family='Helvetica')

plt.figure(1,figsize=(8,8))
plt.clf()

ax1=plt.subplot2grid((5,2),(2,0),rowspan=3,colspan=2)
plotBranch(ax1,bifTable[br1S],'r-',2,maxgap=0.5)
plotBranch(ax1,bifTable[br1U],'k-',1,maxgap=1)

print plotBranch(ax1,bifTable[PerStBranch],'g-',hilo=True,maxgap=5,width=1)
print plotBranch(ax1,bifTable[PerUnstBranch],'b-',hilo=True,maxgap=5,width=1)    

#plt.plot((0.2,0.2,0.24,0.24,0.2),(-55,-45,-45,-55,-55),'k:')
#plt.plot((0.244,0.244,0.25,0.25,0.244),(-40,-20,-20,-40,-40),'k:')

HBpoint=((0.211772,-69.7123))
LPpoint=((0.318099,-31.6646))

plt.annotate('HB',xy=HBpoint,xytext=(HBpoint[0]-0.01,HBpoint[1]+5),
                 arrowprops=dict(facecolor='black', arrowstyle='->')) 
plt.annotate('LP',xy=LPpoint,xytext=(LPpoint[0]-0.01,LPpoint[1]+5),
                 arrowprops=dict(facecolor='black', arrowstyle='->'))                  

#
#plt.plot((0.2,0.2,0.24,0.24,0.2),(-55,-45,-45,-55,-55),'k:')
#plt.plot((0.244,0.244,0.25,0.25,0.244),(-40,-20,-20,-40,-40),'k:')

ax1.set_ylabel("Voltage (mV)")
ax1.set_xlabel(u"$\mathsf{g_{sd}(mS/cm^2)}$",fontsize='large')

ax1.set_xlim((0.15,0.33))
ax1.set_ylim((-90,0))
#ax1.text(0.2,-45,'C1',ha='center',va='bottom')
#ax1.text(0.244,-20,'C2',ha='center',va='bottom')

ax1i=plt.axes([0.17,0.4,0.2,0.15], axisbg='w')
print(plotBranch(ax1i,bifTable[br1S],'r-',2,maxgap=0.5))
print(plotBranch(ax1i,bifTable[br1U],'k-',1,maxgap=1))
ax1i.plot((0.15,0.34,0.34,0.15,0.15),(-90,-90,0,0,-90),'k:',lw=0.5)

ax1i.set_xlim((0.1,1))
ax1i.set_ylim((-100,20))
ax1i.set_xticks([0,0.5,1])
ax1i.set_yticks([-80,-40,0])

ax2=plt.subplot2grid((5,2),(0,0),rowspan=2,colspan=2)
plt.scatter(ISIdata[:,0]/1e4,ISIdata[:,1],c=ISIdata[:,2],
            marker='.',edgecolors='none',cmap='jet',vmin=0,vmax=1)
ax2.set_ylim((200,1000))
ax2.set_yscale('log')
#ax2.tick_params(axis='y',which='minor',width=0)


ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#ax2.get_yaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
ax2.set_yticks([200,500,1000])
ax2.set_ylabel("interval (ms)")

ax2.set_xlim((0.15,0.33))
#ax2.set_xticklabels([])
ax2.set_xlabel(u"$\mathsf{g_{sd}(mS/cm^2)}$",fontsize='large')
pos=plt.axes([0.15,0.7,0.02,0.25])
plt.colorbar(cax=pos,ticks=(0,0.5,1))

"""

ax3=plt.subplot2grid((7,2),(5,0),rowspan=2)
for cond in PerStBranches:
    plotBranch(ax3,bifTable[cond],'g-',hilo=True,maxgap=0.5,width=1)
    
for cond in PerUnstBranches:
    plotBranch(ax3,bifTable[cond],'b-',hilo=True,maxgap=0.5,width=1)    
#
PDs=[(0.209072,-54.0793),(0.215991,-51.2519),(0.217399,-50.5465),(0.217653,-50.389),
     (0.214,-47.78),(0.214332,-47.413),(0.214351,-47.375)]

for pd in PDs:
    plt.annotate('',xy=pd,xytext=(pd[0]-0.002,pd[1]+0.5),
                 arrowprops=dict(facecolor='black', arrowstyle='->'))    
    
    
ax3.set_xlim((0.2,0.24))
ax3.set_ylim((-55,-45))
ax3.set_xticks((0.2,0.22,0.24))
ax3.set_ylabel("Voltage (mv)")

ax4=plt.subplot2grid((7,2),(5,1),rowspan=2)
for cond in PerStBranches:
    plotBranch(ax4,bifTable[cond],'g-',hilo=True,maxgap=0.5,width=1)
    
for cond in PerUnstBranches:
    plotBranch(ax4,bifTable[cond],'b-',hilo=True,maxgap=0.5,width=1)    

plt.xlim((0.244,0.25))
plt.ylim((-40,-20))
ax4.set_xticks((0.244,0.247,0.25))

plt.figtext(0.02,0.97,'A',fontsize='xx-large')
plt.figtext(0.02,0.7,'B',fontsize='xx-large')
plt.figtext(0.02,0.25,'C',fontsize='xx-large')
"""
plt.tight_layout()