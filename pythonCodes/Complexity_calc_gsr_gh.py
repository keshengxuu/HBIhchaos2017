# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 15:33:22 2016

@author: jmaidana
"""


from __future__ import division
import numpy as np
import complexity_test
from mpi4py import MPI
import time as time_comp

comm=MPI.COMM_WORLD
rank=comm.Get_rank()
threads=comm.Get_size()


#folder="spikes_gsdgsrT36gh06"
#folder="spikes_ghgsdT36"
folder="spikes_ghgsrT36"

Table=np.loadtxt(folder + "/dataTable.txt", delimiter='\t')

Table=Table[np.lexsort((Table[:,1],Table[:,0]))]
#LyapTable=[]
L=len(Table)

for i in range(L):
    row=Table[i]

    if row[2]>0.25 and row[2]<50 and i%threads==rank:
#        spikes=np.loadtxt(folder + "/spikes-gsd%04.0f-gsr%04.0f.txt"%(row[0]*1e7,row[1]*1e7))
#        spikes=np.loadtxt(folder + "/spikes-gh%04.0f-gsd%04.0f.txt"%(row[0]*1e6,row[1]*1e6))
        spikes=np.loadtxt(folder + "/spikes-gsd%04.0f-gsr%04.0f.txt"%(row[0]*1e7,row[1]*1e7))
#        spikes=np.loadtxt(folder + "/spikes-tsd%04.0f-th%04.0f.txt"%(row[0]*100,row[1]*10))
#        spikes=np.loadtxt(folder + "/spikes-tsr%04.0f-th%04.0f.txt"%(row[0]*10,row[1]*10))
        
        
        
        ISI = np.diff(spikes)
        
        if len(ISI)>3500: # for L.E. calculations was 5000
            ISIs2=ISI[:3500]
        else:
            ISIs2=ISI
        tic = time_comp.clock()
        try:    
            cval = complexity_test.complex_spk(np.cumsum(ISIs2))
        except:
            cval = np.array([0,0,0,0,0,0,0,0,0,0])
#        print 'algunos_valores:',cval, np.shape(cval)
        with open(folder + '/ComplexTable-%02d.txt'%rank,'a') as dfile:
            dfile.write('%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n'%(row[0],row[1],cval[0],cval[1],cval[2],cval[3],cval[4],cval[5],cval[6],cval[7],cval[8],cval[9]))
        toc = time_comp.clock()
        print 'time spent: ', (toc- tic) 
