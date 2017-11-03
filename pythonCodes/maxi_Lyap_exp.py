# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 14:56:13 2015
The Huber_braun neuronal model function
@author: ksxuu
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os.path
from matplotlib import colors

from matplotlib import rc, cm
import matplotlib.gridspec as gridspec
#for  MPI
from mpi4py import MPI  
comm=MPI.COMM_WORLD
numproc=comm.size
rank = comm.Get_rank()

def run_kut4(F,t,y,dt,args1):
        K0 = dt*F(y,t,args1)
        K1 = dt*F(y + K0/2.0,t + dt/2.0, args1)
        K2 = dt*F( y + K1/2.0,t + dt/2.0,args1)
        K3 = dt*F( y + K2,t + dt,args1)
        return (K0 + 2.0*K1 + 2.0*K2 + K3)/6.0


def HyB(Var,t,tempF):
    [rho,phi]=tempF
    [v,asd,asr,ah]=Var
    isd = rho*gsd*asd*(v - Ed)                                                                          # ecuacion II para I
    Imemb=isd + rho*gsr*(asr**2)/(asr**2+0.4**2)*(v-Er) \
                + rho*gh*ah*(v - Eh)+ rho*gl*(v - El)    # ecuacion I de la membrana 
    asdinf = 1/(1+np.exp(-zsd*(v-V0sd)))                                                                # ecuacion IV para asd
    ahinf= 1/(1+np.exp(-zh*(v-V0h)))                                                                    # ecuacion IV para ah
#setting the equation
    Det=np.array([-Imemb,
                phi*(asdinf - asd)/tsd,
                phi*(-eta*isd - kappa*asr)/tsr,
                phi*(ahinf-ah)/th])
    return Det


"""
The main function is starting from here          
"""
#Parámetros de este Modelo 
gsd = 0.21; gsr = 0.28;
gl = 0.06; gh = 0.4;
V0sd = -40; zsd = 0.11; tsd = 10;
eta = 0.014; kappa = 0.18; tsr = 35;
V0h= -85; zh = -0.14; th=125;

Ed = 50; Er = -90; El = -80; Eh = -30;


v=-60 #Fijamos arbitrariamente un voltaje inicial
t_trans=5000
delta_t=0.01

temp=36
phi = 3**((temp-25.)/10)
rho = 1.3**((temp-25.)/10)
tempF=[rho,phi]
#Luego calulamos el valor de las variables a ese voltaje
asd = 1/(1+np.exp(-zsd*(v-V0sd)));
ah= 1/(1+np.exp(-zh*(v-V0h)));

asr = -eta*rho*gsd*asd*(v - Ed)/kappa;                  # ecuacion V para asr



testep=30
for i in range(0,testep):
    if i%numproc == comm.rank:
        gsd=0.14+i*0.001
        #Luego calulamos el valor de las variables a ese voltaje
#Luego calulamos el valor de las variables a ese voltaje
        asd = 1/(1+np.exp(-zsd*(v-V0sd)));                                                                      # ecuacion IV para asd                                   
        ah= 1/(1+np.exp(-zh*(v-V0h)))  
        
        rho = 1.3**((temp-25.)/10)
        phi = 3**((temp-25.)/10)
        asr = -eta*rho*gsd*asd*(v - Ed)/kappa;
        tempF=[rho,phi]
        
        #initial conditions
        y0=[v,asd,asr,ah]

#    
#     Function to calculate the maximal Lyapunov exponent of
#     an attractor of a system dy/dt=f(y, t) (it is implicitly 
#     assumed that f(y, t) is independent of t)
#     Inputs: dydt - handle to a function that calculates dy/dt
#             y0 - initial condition in the basin of attraction
#             t_trans - time for transients to disappear
#             d0 - initial separation of two orbits
#             delta_t - time step
#             t_max - length of time (after t_trans) to integrate for
#                     (steps 3-5 are repeated until this time is reached)
#     Outputs: mle - running average of the maximal Lyapunov
#                    exponent at each time step
# integrate to get rid of transient behaviour:
        time = np.arange(0,t_trans,delta_t)
        Y=integrate.odeint(HyB, y0,time,args = ((rho,phi),))


        d0=0.0001
        t_max=10000
        y1 = Y[-1,:] #'#;    % final value of solution is y1
        y2=y1+np.append(d0,np.zeros(len(y1)-1));   # perturb by d0 to get y2
        N_steps =int(t_max/delta_t);   # number of steps required
        sl = np.zeros(N_steps);
        sl0=0
        sum0 = 0;
        t=0
        #   # integrate both orbits by a time delta_t:
        for I in range(N_steps):
        #    y = y + run_kut4(HyB,t,y,dt,tempF)
            y1 =y1+run_kut4(HyB,t,y1,delta_t,tempF);
            y2 =y2+run_kut4(HyB,t,y2,delta_t,tempF);
            t=t+delta_t
            d1 =np.linalg.norm(y2-y1);              # new separation
            lambda0 =np.log(d1/d0)/delta_t;   # Lyapunov exponent
            sum0 = sum0+lambda0;        # running sum of Lyapunov exponents
#            sl[I] = sum0/I;
            y2 = y1+(y2-y1)*d0/d1;   # renormalise y2 so separation is d0
        #end
        ## divide running sum by number of iterations to get running average:
#        mle = sl[0:I];
        mle=sum0/N_steps
        outfile = open("gsd%0.3f_vs_gh.txt" % gsd, "a")
#        heads='gsd\tgsr\tfiringRate'
        for k0 in range(2):
            if k0==0:
                outfile.write(str(gsd) + " ")
            else:
                outfile.write(str(mle) + "\n")
                
        #            print le,rank
        outfile.close()






#plt.plot(Y[:,0],Y[:,2])
