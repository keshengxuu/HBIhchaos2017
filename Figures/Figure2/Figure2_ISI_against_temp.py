"""
Created on Mon Mar  7 21:26:41 2016
keshengxuu@gmail.com 
@author: ksxuu
"""
import numpy as np
import contextlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rc, cm
import matplotlib as mlp
from numpy.ma import masked_array
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec
from matplotlib import rc, cm
import csv
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1 import make_axes_locatable


plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'



fsize = 10
line_width = 0.7
#change the font size of xtick and  ytick
plt.rc('xtick',labelsize=14)
plt.rc('ytick', labelsize=14)

# define the fucntion to read spike train using CSV_reader
def read_csv_file(filename):
    """Reads a CSV file and return it as a list of rows."""
    data = []
    for row in csv.reader(open(filename),delimiter=','):
# transfer the string type to the float type
        row=list((float(x) for x in row))
        data.append(row)
    return data


Xvals=read_csv_file('Data/Xvals.txt' )
Xvals2=read_csv_file('Data/Xvals2.txt' )
X1vals=np.array([x[0] for x in Xvals])


Yvals=read_csv_file('Data/Yvals.txt' )
Yvals2=read_csv_file('Data/Yvals2.txt' )



Cvals=np.loadtxt('Data/Cvals.txt')
Cvals2=np.loadtxt('Data/Cvals2.txt')


Cmax=np.max(Cvals)

Cmax2=np.max(Cvals2)


fsize = 14
ffsize=18
line_width = 0.7

fig=plt.figure(1, figsize=(10,8))
plt.clf()
#Move the edge of an axes to make room for tick labels
#fig.subplots_adjust(left = 0.1,right = 0.95,bottom=0.1,top = 0.95,wspace = 0.1,hspace = 0.15)

ax1 = plt.subplot(311)
for X,Y,C in zip(Xvals,Yvals,Cvals):
    ax1.scatter(X,Y,s=1,linewidth=0,c=C*np.ones_like(Y),vmin=0,vmax=1.6,alpha=0.9)
plt.xlabel(r'$\mathsf{ Temperature \; (^{\circ} C)}$',fontsize = ffsize)
plt.ylabel(r'$\mathsf{ISI \; (ms) }$',fontsize = ffsize)
plt.xticks([10,15,20,25,30, 35],['10','15','20','25','30','35'])
plt.yscale('log')
plt.grid(True)
plt.ylim([10**1,3*10**3])
plt.xlim(5,37)
plt.text(2.5,10000,'A', fontsize = 20)
plt.text(28,1000,r'$\mathsf{g_{h}=0.4}$', fontsize = 20)

# create an axes on the right side of ax1. The width of cax will be 1%
# of ax1 and the padding between cax and ax will be fixed at 0.05 inch.
divider2 = make_axes_locatable(ax1)
cax2= divider2.append_axes("right","5%", pad="3%")
fig.delaxes(cax2)

zoom1=28.5
zoom2=29.5
i1=min(((X1vals>=zoom1)*(X1vals<=zoom2)).nonzero()[0])
i2=max(((X1vals>=zoom1)*(X1vals<=zoom2)).nonzero()[0])
Xvals1=Xvals[i1:i2]
Yvals1=Yvals[i1:i2]
Cvals1=Cvals[i1:i2]

axins = plt.subplot(323)
for X,Y,C in zip(Xvals1,Yvals1,Cvals1):
    axins.scatter(X,Y,s=3,linewidth=0,c=C*np.ones_like(Y),vmin=0,vmax=1.6,alpha=0.9)
#plt.xlabel(r'$\mathsf{ Temperature \; (^{\circ} C)}$',  fontdict=font)
plt.ylabel(r'$\mathsf{ISI \; (ms) }$',fontsize = ffsize)
plt.xlim(zoom1,zoom2)
plt.xticks([28.5,29,29.5])
#plt.xticks([28,29,30],['28','29','30'])
plt.grid(True)
plt.yscale('log')
#plt.yticks([100,600],['100','600'], fontsize = fsize)
plt.ylim([2*10**1,0.5*10**3])
#plt.text(28.1,800,'C', fontsize = 20)

# create an axes on the right side of ax1. The width of cax will be 1%
# of ax1 and the padding between cax and ax will be fixed at 0.05 inch.
divider3 = make_axes_locatable(axins)
cax3= divider3.append_axes("right","16.5%", pad="3%")
fig.delaxes(cax3)

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax1, axins, loc1=2, loc2=2, fc="none", ec="0.5")
plt.draw()

zoom1=33.7
zoom2=35.3
i1=min(((X1vals>=zoom1)*(X1vals<=zoom2)).nonzero()[0])
i2=max(((X1vals>=zoom1)*(X1vals<=zoom2)).nonzero()[0])
Xvals1=Xvals[i1:i2]
Yvals1=Yvals[i1:i2]
Cvals1=Cvals[i1:i2]

axins1 = plt.subplot(324)
for X,Y,C in zip(Xvals1,Yvals1,Cvals1):
    im=axins1.scatter(X,Y,s=3,linewidth=0,c=C*np.ones_like(Y),vmin=0,vmax=1.6,alpha=0.9)

#plt.xlabel(r'$\mathsf{ Temperature \; (^{\circ} C)}$', fontdict=font)
plt.ylabel(r'$\mathsf{ISI \; (ms) }$',fontsize = ffsize)
plt.xlim(zoom1,zoom2)
plt.xticks([34,35],['34','35'])
#plt.xticks([34,35,36],['34','35','36'])
plt.grid(True)
plt.yscale('log')
#plt.yticks([100,1000],['100','1000'], fontsize = fsize)
plt.ylim([10**2,10**3])

# create an axes on the right side of ax1. The width of cax will be 1%
# of ax1 and the padding between cax and ax will be fixed at 0.05 inch.
divider4 = make_axes_locatable(axins1)
cax4= divider4.append_axes("right","16.5%", pad="3%")
fig.delaxes(cax4)

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes arear'$\mathsf{
mark_inset(ax1, axins1, loc1=2, loc2=2, fc="none", ec="0.5")
plt.draw()


# subplots of  bifurcation of the model without ih
ax4 = plt.subplot(313)
for X,Y,C in zip(Xvals2,Yvals2,Cvals2):
    im0=plt.scatter(X,Y,s=1,linewidth=0,c=C*np.ones_like(Y),vmin=0,vmax=1.6,alpha=0.9)

plt.xlabel(r'$\mathsf{ Temperature \; (^{\circ} C)}$',fontsize = ffsize)
plt.ylabel(r'$\mathsf{ISI \; (ms) }$',fontsize = ffsize)
plt.xlim(5,37)
plt.xticks([10,15,20,25,30, 35],['10','15','20','25','30','35'])
plt.grid(True)
plt.yscale('log')
plt.ylim([10**1,10**4])
plt.text(2.5,10000,'B', fontsize = 20)
plt.text(28,3000,r'$\mathsf{g_{h}=0}$', fontsize = 20)

# create an axes on the right side of ax1. The width of cax will be 1%
# of ax1 and the padding between cax and ax will be fixed at 0.05 inch.
divider1 = make_axes_locatable(ax4)
cax1= divider1.append_axes("right","5%", pad="3%")
fig.delaxes(cax1)




#To create an axes,n axes at position rect [left, bottom, width, height]
cax5 = fig.add_axes([0.922, 0.12, 0.015, 0.8])
cbar=fig.colorbar(im0, extend='max',cax=cax5,ticks=[ 0.0, 0.4, .8,1.2, 1.6])
cbar.set_label(r'$\mathsf{Lyapunov \; Exponent }$',fontsize=15)
cbar.ax.set_yticklabels(['0.0', '0.4', '0.8', '1.2','1.6'],fontsize = fsize)

#subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
#left  = 0.125  # the left side of the subplots of the figure
#right = 0.9    # the right side of the subplots of the figure
#bottom = 0.1   # the bottom of the subplots of the figure
#top = 0.9      # the top of the subplots of the figure
#wspace = 0.2   # the amount of width reserved for blank space between subplots
#hspace = 0.2   # the amount of height reserved for white space between subplots
#plt.subplots_adjust(bottom=0.08,left=0.1,wspace = 0.05,hspace = 0.3,right=0.98, top=0.97)
plt.tight_layout()
plt.savefig('Figure2_ISI_against_temp.png',dpi=600)

