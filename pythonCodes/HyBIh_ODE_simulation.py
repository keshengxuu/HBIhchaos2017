# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:46:06 2016
keshengxuu@gmail.com & patricio.orio@uv.cl 
@author: keshengxu
"""


import numpy as np
from scipy.stats import norm
from scipy import integrate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os.path
from matplotlib import colors

from matplotlib import rc, cm
import matplotlib.gridspec as gridspec

#Forma clásica de matplotlib, pero con la ventaja de poder hacer subscripts

#plt.rc('font', size = 14)
#plt.rcParams['mathtext.fontset'] = 'custom'
#plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
#plt.rcParams['svg.fonttype'] = 'none'

#Fonts de Latex

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


font= { 'family' : 'serif',
#        'color' : 'darkred',
	'weight' : 'normal',
	'size'   : 18,
	}
"""
Aquí definimos una funcion que toma las variables a tiempo t
y devuelve las derivadas
"""
def HyB(Var,t,tempF):
    [rrho,pphi] = tempF
    [v,ar,asd,ca,ah]=Var
    ad = 1/(1+np.exp(-zd*(v-V0d)))
    isd = rrho*gsd*asd*(v - Ed)
    #Imemb=isd + rho*gd*ad*(v - Ed) + rho*(gr*ar + gsr*asr)*(v-Er) + gl*(v - El)
    Imemb = isd + rrho*gd*ad*(v - Ed) + rrho*(gr*ar + gsr*(ca**2)/(ca**2+0.4**2))*(v-Er) + rrho*gl*(v - El) \
                + rrho*gh*ah*(v - Eh)
    arinf = 1/(1+np.exp(-zr*(v-V0r)))
    asdinf = 1/(1+np.exp(-zsd*(v-V0sd)))
    ahinf = 1/(1+np.exp(-zh*(v-V0h)));

    return np.array([-Imemb,
                pphi*(arinf - ar)/tr,
                pphi*(asdinf - asd)/tsd,
                pphi*(-eta*isd - kappa*ca)/tsr,
                pphi*(ahinf-ah)/th])

def ISIS(spkf):
    ISIs = np.diff(spkf)
    return ISIs

#Parámetros del modelo
gd = 2.5; gr = 2.8; gsd = 0.21; gsr = 0.28;
gl = 0.06; gh = 0.4; gtrek=0.06;  #Ojo con gh, el prof alfredo tiene 0.6
V0d = -25; V0r = -25; zd = 0.25; zr = 0.25;tr = 2;
V0sd = -40; zsd = 0.11; tsd = 10;
eta = 0.014; kappa = 0.18; tsr = 35;
V0h= -85; zh = -0.14; th=125;


Ed = 50; Er = -90; El = -80; Eh = -30;

temp=14
rho=1.3**((temp-25.)/10)
phi = 3**((temp-25.)/10)

temp2 = 20
rho2=1.3**((temp2-25.)/10)
phi2 = 3**((temp2-25.)/10)
#Y luego el comando 'ode# -*- coding: utf-8 -*-

temp3 = 24.76
rho3=1.3**((temp3-25.)/10)
phi3 = 3**((temp3-25.)/10)

temp4 = 26
rho4=1.3**((temp4-25.)/10)
phi4 = 3**((temp4-25.)/10)

temp5 = 33
rho5=1.3**((temp5-25.)/10)
phi5 = 3**((temp5-25.)/10)

temp6 = 36.3
rho6=1.3**((temp6-25.)/10)
phi6 = 3**((temp6-25.)/10)


Tstop = 200000
#dt=0.1#ms OJO!!
dt=0.1
#voltaje inicial e inicialización de variables
v = -60
v2 = -60 #Fijamos arbitrariamente un voltaje inicial
#Luego calulamos el valor de las variables a ese voltaje
ad = 1/(1+np.exp(-zd*(v-V0d)));
ar = 1/(1+np.exp(-zr*(v-V0r)));
asd = 1/(1+np.exp(-zsd*(v-V0sd)));
ca = -eta*rho*gsd*asd*(v - Ed)/kappa;
ah = 1/(1+np.exp(-zh*(v-V0h)));

ad2 = 1/(1+np.exp(-zd*(v-V0d)));
ar2 = 1/(1+np.exp(-zr*(v-V0r)));
asd2 = 1/(1+np.exp(-zsd*(v-V0sd)));
ca2 = -eta*rho*gsd*asd*(v - Ed)/kappa;
ah2 = 1/(1+np.exp(-zh*(v-V0h)));
#Ahora viene la simulacion misma
#Creamos un vector con los valores iniciales
X=np.array([v,ar,asd,ca,ah])
X2=np.array([v2,ar2,asd2,ca2,ah2])
time = np.arange(0,Tstop,dt)
time2 = np.arange(0,Tstop,dt)
#Y luego el comando 'ode' se encarga del resto

Var_t = integrate.odeint(HyB, X, time, args = ((rho,phi),))
Var_t_2 = integrate.odeint(HyB, X, time, args = ((rho2,phi2),))
Var_t_3 = integrate.odeint(HyB, X, time, args = ((rho3,phi3),))
Var_t_4 = integrate.odeint(HyB, X, time, args = ((rho4,phi4),))
Var_t_5 = integrate.odeint(HyB, X, time, args = ((rho5,phi5),))
Var_t_6 = integrate.odeint(HyB, X, time, args = ((rho6,phi6),))


#Var_t contiene el curso temporal de todas las variables,
#en los tiempos especificados por el vector [dt:dt:Tstop]
#
aa = 400000
#bb = Tstop/dt-1 Antiguo
bb = 600000
cc = 680000


def SPKF(serie):
    #deteción de spikes
    spk_ind= np.where(np.diff((serie[aa:-1,0]>0)*1)==1)[0]+aa+1
    #interpolación de tiempo de spikes
    spkf=time[spk_ind]-serie[spk_ind,0]*dt/(serie[spk_ind+1,0]-serie[spk_ind,0])
    return spkf

def ISIS(spkf):
    ISIs = np.diff(spkf)
    return ISIs

spkf1 = SPKF(Var_t)
ISI1 = ISIS(spkf1)

spkf2 = SPKF(Var_t_2)
ISI2 = ISIS(spkf2)

spkf3 = SPKF(Var_t_3)
ISI3 = ISIS(spkf3)

spkf4 = SPKF(Var_t_4)
ISI4 = ISIS(spkf4)

spkf5 = SPKF(Var_t_5)
ISI5 = ISIS(spkf5)

spkf6 = SPKF(Var_t_6)
ISI6 = ISIS(spkf6)

np.savetxt('system_data/time_all_simulation.txt', time)

np.savetxt('system_data/5_variables_series/V_1.txt', Var_t)
np.savetxt('system_data/spkf1.txt', spkf1)
np.savetxt('system_data/ISI1.txt', ISI1)

np.savetxt('system_data/5_variables_series/V_2.txt', Var_t_2)
np.savetxt('system_data/spkf2.txt', spkf2)
np.savetxt('system_data/ISI2.txt', ISI2)

np.savetxt('system_data/5_variables_series/V_3.txt', Var_t_3)
np.savetxt('system_data/spkf3.txt', spkf3)
np.savetxt('system_data/ISI3.txt', ISI3)

np.savetxt('system_data/5_variables_series/V_4.txt', Var_t_4)
np.savetxt('system_data/spkf4.txt', spkf4)
np.savetxt('system_data/ISI4.txt', ISI4)

np.savetxt('system_data/5_variables_series/V_5.txt', Var_t_5)
np.savetxt('system_data/spkf5.txt', spkf5)
np.savetxt('system_data/ISI5.txt', ISI5)

np.savetxt('system_data/5_variables_series/V_5.txt', Var_t_6)
np.savetxt('system_data/spkf6.txt', spkf6)
np.savetxt('system_data/ISI6.txt', ISI6)

#binss = 30
##linf = 20  Antiguo
#linf = 0
#lsup = 800
#linf2 = 0
#lsup2 = 350
#lsup3 = 260
#fsize = 10
#
#
#plt.figure(7, figsize=(8,6))
#
#gs1 = gridspec.GridSpec(5, 3)
##v17
#gs1.update(wspace=0.37, hspace=0.5, left = 0.08, right = 0.87,bottom = 0.1, top = 0.9)

#ax1 = plt.subplot(gs1[0])
#plt.plot(time[bb:cc],Var_t[bb:cc,0], color= 'chocolate', linewidth = 1) #OJO modiificado para dt = 0.05
##plt.plot(time[800000:bb],Var_t[800000:bb,0], color= 'sienna', linewidth = 1) #OJO modiificado para dt = 0.05
##plt.xlim(40000,45000)
#plt.ylim(-80,28)
#plt.xticks([])
##plt.yticks([])
#plt.yticks([-70,-35,0,20],['-70','-35','0','20'], fontsize = fsize)
#
#ax12 = plt.subplot(gs1[1])
#plt.plot(ISIS(spkf = SPKF(Var_t)), 'bo', ms=1)
##plt.yticks([0,300,600],['0','300','600'], fontsize = fsize)
##plt.xlim(linf2,lsup2)
##plt.ylim(0,700)
##plt.xticks([])
#plt.yscale('log')
#
#ax13 = plt.subplot(gs1[2])
#plt.hist(ISIS(spkf = SPKF(Var_t)), bins=binss, range = (linf, lsup))
#plt.yticks([0,200,400],['0','200','400'], fontsize = fsize)
#plt.xlim(linf,lsup)
#plt.ylim(0,410)
#plt.xticks([])
#
#ax111 = ax13.twinx()
#plt.hist(ISIS(spkf = SPKF(Var_t)), bins=binss, range = (linf, lsup))
#plt.ylim(0,410)
#ax111.yaxis.tick_right()
#y = [200]
#labels = [r'$T=14^{\circ} C$']
#plt.yticks(y, labels, fontsize= 14, rotation='horizontal')



np.savetxt('system_data/time_voltage_trace.txt', time[bb:cc])
np.savetxt('system_data/var_t1.txt', Var_t[bb:cc,0])
np.savetxt('system_data/var_t2.txt', Var_t_2[bb:cc,0])
np.savetxt('system_data/var_t3.txt', Var_t_3[bb:cc,0])
np.savetxt('system_data/var_t4.txt', Var_t_4[bb:cc,0])
np.savetxt('system_data/var_t5.txt', Var_t_5[bb:cc,0])
np.savetxt('system_data/var_t6.txt', Var_t_6[bb:cc,0])



#plt.savefig('46v17.png', dpi = 300)
#plt.savefig('001.png', dpi = 300)
#plt.savefig('002.png', dpi = 300)
