# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 13:28:33 2015

@author: porio
"""

from __future__ import division
import numpy as np
import NPE
from mpi4py import MPI

comm=MPI.COMM_WORLD
rank=comm.Get_rank()
threads=comm.Get_size()


folder="spikes_gsdT36"


Table=np.loadtxt(folder + "/dataTable.txt", delimiter='\t')
#Table=np.insert(Table,1,np.zeros(len(Table)),axis=1)

Table=Table[np.lexsort((Table[:,1],Table[:,0]))]
#Table=Table[np.lexsort((Table[:,0]))]
LyapTable=[]
L=len(Table)

for i in range(L):
    row=Table[i]
    if row[2]>0.25 and row[2]<50 and i%threads==rank:
#        spikes=np.loadtxt(folder + "/spikes-th%04d.txt"%(row[0]*10))        
#        spikes=np.loadtxt(folder + "/spikes-T%06.3f.txt"%row[0])        
#        spikes=np.loadtxt(folder + "/spikes-gsd%04.0f-gsr%04.0f.txt"%(row[0]*1e7,row[1]*1e7))
#        spikes=np.loadtxt(folder + "/spikes-gh%04.0f-gsd%04.0f.txt"%(row[0]*1e7,row[1]*1e7))
#        spikes=np.loadtxt(folder + "/spikes-th%04.0f-gsd%04.0f.txt"%(row[0]*10,row[1]*1e7))
        spikes=np.loadtxt(folder + "/spikes-gsd%04.0f.txt"%(row[0]*1e7))
#        spikes=np.loadtxt(folder + "/spikes-tsd%04.0f-th%04.0f.txt"%(row[0]*100,row[1]*10))
#        spikes=np.loadtxt(folder + "/spikes-tsr%04.0f-th%04.0f.txt"%(row[0]*10,row[1]*10))
        
        ISI = np.diff(spikes)
        
        if len(ISI)>5000:
            ISIs2=ISI[:5000]
        else:
            ISIs2=ISI
        L1,_,p1,N1,_=NPE.Lyap2(ISIs2,d=7,eps=0.1,maxpercent=0.01,plot=False,npoints=6)
        L2,_,p2,N2,_=NPE.Lyap2(ISIs2,d=9,eps=0.1,maxpercent=0.01,plot=False,npoints=6)
        L3,_,p3,N3,_=NPE.Lyap2(ISIs2,d=11,eps=0.1,maxpercent=0.01,plot=False,npoints=6)
        with open(folder + '/LyapTable-%02d.txt'%rank,'a') as dfile:
            dfile.write('%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n'%(row[0],row[1],L1,p1,L2,p2,L3,p3))
    
