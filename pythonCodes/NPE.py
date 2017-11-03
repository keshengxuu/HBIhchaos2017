# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 21:13:18 2015

@author: porio
"""

import numpy as np
import matplotlib.pyplot as plt
#import tentMap
import itertools as it
from scipy import stats

def neighbors(x, serieRm, N):
    """
    Calculate N neighbors of x within serie
    
        x:      vector
        serie:  array of vectors with same lenghth as x
       N:      number of neighbors to return
    """
    distances = np.sqrt(np.sum((x-serieRm)**2,1))
    S=np.argsort(distances)
    ind_neigh=S[1:N+1]    
    return ind_neigh, serieRm[ind_neigh]

def Lyap(serie,d=4,eps=10,npoints=10,plot=False,maxpercent=0.5,verbose=True):
    Rm = np.array([serie[i:i+d] for i in range(len(serie)-d+1)])
    N=float(len(Rm)-npoints-2)
    Ndists=N**2 - N   #Number of effective distances (point pairs)
    alldist=np.sqrt(np.sum((Rm[None,:,:] - Rm[:,None,:])**2,axis=-1))  
    pairs1=np.array(((alldist<eps)*(alldist>0)).nonzero()).T
    pairs=pairs1[(np.all(pairs1<N,1)*(np.abs(pairs1[:,0] - pairs1[:,1])>10)).nonzero()]

    #Reduce epsilon until only a 'maxpercent' fraction of pairs is considered    
    while len(pairs)/Ndists > maxpercent/100:
        eps*=0.9
        pairs1=np.array(((alldist<eps)*(alldist>0)).nonzero()).T
        pairs=pairs1[(np.all(pairs1<N,1)*(np.abs(pairs1[:,0] - pairs1[:,1])>10)).nonzero()]
    #Increase epsilon to have at least 10 pairs    
    while len(pairs)<10:
        eps*=1.1
        pairs1=np.array(((alldist<eps)*(alldist>0)).nonzero()).T
        pairs=pairs1[(np.all(pairs1<N,1)*(np.abs(pairs1[:,0] - pairs1[:,1])>10)).nonzero()]
    if verbose:
        print len(pairs),"point pairs out of",Ndists,"%.2g percent"%(100*float(len(pairs))/Ndists),"epsilon=",eps

    idx=idx=pairs[:,None,:] + np.arange(npoints)[:,None]
    dd=np.sqrt(np.sum((Rm[idx][:,:,0,:]-Rm[idx][:,:,1,:])**2,-1))
    ddm=np.mean(dd,0)
    slope,_,rval,pval,_=stats.linregress(range(1,npoints),np.log(ddm[1:]))
    if plot:
        plt.figure(10)
        plt.semilogy(ddm,'.')
    return slope,rval,pval,ddm    
    
def Lyap2(serie,d=4,eps=10,npoints=10,plot=False,maxpercent=0.5,min_eps=1e-9,verbose=False):
    """Output: [L,r,p,N,d]
    L      Lyapunov exponent
    r      r coefficient of linear regression
    p      p-value associated to linear regression
    N      final number of pairs used
    d      array containing the mean distances calculated
    """
    
    #Calculate reconstructed serie in dimension d    
    Rm = np.array([serie[i:i+d] for i in range(len(serie)-d+1)])
    #Leave the last 'npoints' out of the analysis    
    N=float(len(Rm)-npoints-2)
    Ndists=N**2 - N   #Number of effective distances (point pairs)
    Ns=len(Rm)-1
    #Calculate all the distances
#    alldist=np.sqrt(np.sum((Rm[None,:,:] - Rm[:,None,:])**2,axis=-1))  
    interm_mat=np.array([serie[None,i:i+Ns]-serie[i:i+Ns,None] for i in range(d)])    
    alldist=np.sqrt(np.sum(interm_mat**2,0))
    
    #Find the indices of the distances lower than epsilon and that are not zero    
    pairs1=np.array(((alldist<eps)*(alldist>0)).nonzero()).T
    #Eliminate the pairs that are close in time and those that involve points larger than N
    pairs=pairs1[(np.all(pairs1<N,1)*(np.abs(pairs1[:,0] - pairs1[:,1])>10)).nonzero()]
    
    #Reduce epsilon until only a 'maxpercent' fraction of pairs is considered    
    ii=0
    while len(pairs)/Ndists > maxpercent/100 and eps>min_eps:
        eps*=0.75
        pairs1=np.array(((alldist<eps)*(alldist>0)).nonzero()).T
        pairs=pairs1[(np.all(pairs1<N,1)*(np.abs(pairs1[:,0] - pairs1[:,1])>10)).nonzero()]
        if verbose:
            print len(pairs),"point pairs out of",Ndists,"%.2g percent"%(100*float(len(pairs))/Ndists),"epsilon=",eps
    #Increase epsilon to have at least 10 pairs    
    while len(pairs)<10:
        eps*=1.5
        pairs1=np.array(((alldist<eps)*(alldist>0)).nonzero()).T
        pairs=pairs1[(np.all(pairs1<N,1)*(np.abs(pairs1[:,0] - pairs1[:,1])>10)).nonzero()]
        if verbose:
            print len(pairs),"point pairs out of",Ndists,"%.2g percent"%(100*float(len(pairs))/Ndists),"epsilon=",eps
    
    print len(pairs),"point pairs out of",Ndists,"%.2g percent"%(100*float(len(pairs))/Ndists),"epsilon=",eps

    #idx contains the sequences whose distance will be calculated. 
    #It is a (Npairs x npoints x 2) array
    #For each pair, each of its components is extended npoints into the future
    idx=pairs[:,None,:] + np.arange(npoints)[:,None]
    #Then calculate the distances, only of the last component
    dd=np.abs(Rm[idx][:,:,0,-1]-Rm[idx][:,:,1,-1])
    # and the mean across the pairs
    ddm=np.mean(dd,0)
    #fit the log(distance) to a line - discard the first distance
    slope,_,rval,pval,_=stats.linregress(range(1,npoints),np.log(ddm[1:]))
    #and plot if choosed    
    if plot:
        plt.figure(10)
        plt.semilogy(ddm,'.')
    return slope,rval,pval,len(pairs),ddm      

def dists(serie,d=4):
    Rm = np.array([serie[i:i+d] for i in range(len(serie)-d+1)])
    return Rm,np.sqrt(np.sum((Rm[None,:,:] - Rm[:,None,:])**2,axis=-1))  
    
    
def Error(serieRm,n_frac=0.01,index=0,horizon=4):
    horizon=np.array(horizon)
    neighI,neigh = neighbors(serieRm[index],serieRm,len(serieRm)*n_frac)  # (NNeigh x d)

    indNc = neighI[neighI+max(horizon) < len(serieRm)-1]

    NLP=np.mean(serieRm[indNc[:,None] + horizon],0)      #  (h x d)
    error = np.sqrt(np.sum((serieRm[index+horizon] - NLP)**2,-1))
    Nerror = np.sqrt(np.sum((serieRm[index+horizon] - np.mean(serieRm,0))**2,-1))    
    return error/Nerror

def Error2(serieRm,n_frac=0.01,index=0,horizon=4):
    horizon=np.array(horizon)
    neighI,neigh = neighbors(serieRm[index],serieRm,len(serieRm)*n_frac)  # (NNeigh x d)

    indNc = neighI[neighI+max(horizon) < len(serieRm)-1]

    NLP=np.mean(serieRm[indNc[:,None] + horizon],0)      #  (h x d)
    error = np.abs(serieRm[index+horizon,-1] - NLP[:,-1])
    Nerror = np.abs(serieRm[index+horizon,-1] - np.mean(serieRm[:,-1]))
    return error/Nerror

def NPE(serie,hor=range(1,11),d=4):
    Rm = np.array([serie[i:i+d] for i in range(len(serie)-d+1)])
    NPEs = np.array([Error(Rm,index=i,horizon=hor) for i in range(len(Rm)-max(hor))])

    return np.mean(NPEs,0)#/np.mean(alldist)
    
def NPE2(serie,hor=range(1,11),d=4):
    Rm = np.array([serie[i:i+d] for i in range(len(serie)-d+1)])
    
    NPEs = np.array([Error2(Rm,index=i,horizon=hor,n_frac=frac) for i in range(len(Rm)-max(hor))])

    return np.mean(NPEs,0)#/np.mean(alldist)

#    return np.mean(NPEs,0), Rm#/np.mean(alldist)
        
    
filtered_ISIs =[]

if __name__=='__main__':
    
    import os
    
    spikes=np.loadtxt('Neuron/spikes_T/spikes-T36.520.txt')
    ISIs=np.diff(spikes)
    Lyap2(ISIs,d=5,plot=True,eps=0.005)

    h=range(1,20)
    d=5
    frac=0.005

#    """    
#    plt.figure(1)
#    plt.clf()
#    plt.figure(2)
#    plt.clf()
    pli=it.count(start=1)
    filelist=os.listdir('spikes/')
    filelist.sort()
    spikes=[np.loadtxt('spikes/'+filen) for filen in filelist]
#    NPEs=np.array([NPE(np.diff(spk),hor=h,d=d) for spk in spikes])
#        
#    plt.figure(1)
#    plt.plot(h,NPEs.T)
#    plt.figure(2)
#    for spk in spikes:
#        plt.subplot(1,4,pli.next())
#        plt.plot(spk[1:],np.diff(spk),'.')
#        
#        
#    plt.figure(1)
#    plt.legend(filelist,loc=0,fontsize='small')

#    """
        
#    serie4 = np.random.normal(scale= 0.01, size=2000)
#    serie3 = np.sin(np.linspace(0,100,2000))
#    serie = tentMap.Tent(L=2000)
#    serieR=np.random.permutation(serie)
#    serie2 = np.random.normal(size=2000)
#    h=range(1,11)
#    d=3
#    aaa = Lyap2(serie,d=5,plot=True)

    ISIs = [np.diff(spk) for spk in spikes]
    filtered_ISIs = [spk[(spk>160)*(spk<270)] for spk in spikes]
    CVs=[]
    LyapN=[]
    pvals=[]
    Npairs=[]

    
    for filISI in filtered_ISIs:
        if len(filISI)>10:
            CVs.append(np.std(filISI)/np.mean(filISI))
        else:
            CVs.append(0)
        
    for ISIk in ISIs:
        L,_,p,Np,_=Lyap2(ISIk,d=12,eps=0.5,npoints=7,verbose=False)
        print L,p,Np
        LyapN.append(L)
        pvals.append(p)
        Npairs.append(Np)

    data=np.array([re.findall(r'[\d.]+(?:e[+-]?\d+)?', fileN.replace('.txt','')) for fileN in filelist]).astype('float')
    
    Table=np.hstack((data,np.array(CVs)[:,None],np.array(LyapN)[:,None],np.array(pvals)[:,None],np.array(Npairs)[:,None]))
    
    np.savetxt('Result.txt',Table,fmt='%.6g',delimiter=',')
    
    
#    NPEs = NPE2(serie[100:],hor=h,d=d)
#    NPEs2 = NPE2(serie2[100:],hor=h,d=d)
#    NPEs3 = NPE2(serie3[100:],hor=h,d=d)
#    NPEsR = NPE2(serieR[100:],hor=h,d=d)
#    NPEs4 = NPE2(serie4[100:],hor=h,d=d)

#a=NPE(serie,h=6,n_frac=0.01)
#result=[np.sqrt(np.mean(NPE(serie,h=hor,n_frac=0.01)**2)) for hor in h]
#
#    plt.figure(3)    
#    plt.clf()    
##    plt.plot(h,NPEs)
#    plt.plot(h,np.squeeze(NPEs))
#    plt.plot(h,NPEs2)
#    plt.plot(h,NPEs3)
#    plt.plot(h,NPEsR)
#    plt.plot(h,NPEs4)
#    plt.legend(['Tent','rand','sine','Tent-R', 'yami'])

    
