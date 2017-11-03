"""  
keshengxuu@gmail.com 
@author: keshengxu
"""

from PyDSTool import *

#introduction
#DSargs = args()                   # create an empty object instance of the args class, call it DSargs
#DSargs.name = 'SHM'               # name our model
#DSargs.ics = icdict               # assign the icdict to the ics attribute
#DSargs.pars = pardict             # assign the pardict to the pars attribute
#DSargs.tdata = [0, 20]            # declare how long we expect to integrate for
#DSargs.varspecs = vardict         # assign the vardict dictionary to the 'varspecs' attribute of DSargs
# "Computational Cell Biology", Fall (Type II)
pars = {'gsd': 0.0,
        'gsr': 0.28,
        'gl': 0.06,
        'gh': 0.25,
        'V0d': -25,
        'V0r': -25,
        'zd':0.25,
        'zr': 0.25,
        'tr': 2.,
        'V0sd': -40,
        'zsd': 0.11,
        'tsd': 10,
        'eta': 0.014,
        'kappa': 0.18,
        'tsr': 35,
        'V0h': -85,
        'zh':-0.14,
        'th': 125,
        'Ed': 50,
        'Er': -90,
        'El': -80,
        'Eh': -30,
        'temp': 36}

icdict = {'v': 0.,'asd': 0.06,'ca': 0.002,'ah': 0.027}

# Set up model
#defines two auxiliary functions for use in those specifications
auxfndict = {'asdinf': (['v'], '1/(1+exp(-zsd*(v-V0sd)))'), \
			 'ahinf': (['v'], '1/(1+exp(-zh*(v-V0h)))')\
			}

vstr = '-1.3**((temp-25.)/10)*(gsd*asd*(v - Ed)+ \
gsr*(ca**2)/(ca**2+0.4**2)*(v-Er) + gl*(v - El)+gh*ah*(v - Eh))'
asdstr = '3**((temp-25.)/10)*(asdinf(v) - asd)/tsd'
castr = '3**((temp-25.)/10)*(-eta*1.3**((temp-25.)/10)*gsd*asd*(v - Ed) - kappa*ca)/tsr'
ahstr = '3**((temp-25.)/10)*(ahinf(v)-ah)/th'

DSargs = args(name='HB_Ih')
DSargs.pars = pars
DSargs.varspecs = {'v': vstr, 'asd': asdstr, 'ca': castr, 'ah': ahstr}
DSargs.fnspecs = auxfndict
DSargs.ics = icdict

ode = Generator.Vode_ODEsystem(DSargs)

#The following curve types are handled by PyCont:
#EP-C Equilibrium point curve
#LP-C Limit point curve/Fold curve/Saddle-node curve
#H-C1 Hopf point curve (method 1)
#H-C2 Hopf point curve (method 2)
#FP-C Fixed point curve (discrete systems)
#LC-C Limit cycle curve (interface to AUTO)
# Set up continuation class
PyCont = ContClass(ode)

PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['gsd']
PCargs.StepSize = 1e-2
PCargs.MaxNumPoints = 200
PCargs.MaxStepSize = 1.
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PCargs.SaveEigen = True
PyCont.newCurve(PCargs)

print('Computing curve...')
start = clock()
PyCont['EQ1'].forward()
print('done in %.3f seconds!' % (clock()-start))
PCargs.name = 'LC1'
PCargs.type = 'LC-C'
PCargs.initpoint = 'EQ1:H1'
PCargs.MinStepSize = 0.005
PCargs.MaxStepSize = 1.0
PCargs.StepSize = 0.01
PCargs.MaxNumPoints = 220
PCargs.NumSPOut = 40;
PCargs.LocBifPoints = 'LPC'
PCargs.SolutionMeasures = 'avg'
PCargs.SaveEigen = True
PyCont.newCurve(PCargs)

print('Computing curve...')
start = clock()
PyCont['LC1'].forward()
print('done in %.3f seconds!' % (clock()-start))

# Limit cycle curve
# Plot
PyCont.display(('gsd','v'),stability=True)
#PyCont['LC1'].display(('gsd','v_min'),stability=True)

plt.xlim([0.1, 0.5])
plt.ylim([-75, 10])
PyCont.plot.fig1.axes1.axes.set_title('Bifurcation Diagram')
plt.savefig('Bifurcation_Diagram.pdf')


#PyCont['LC1'].plot_cycles(figure='fig2', method='stack', exclude='P2', tlim='5T')
#PyCont.plot.fig2.axes1.axes.set_title('Cycles')
#plt.savefig('Cycles.pdf')
#PyCont['EQ1'].display(('v','asd'),stability=True, figure='new')
#
#PyCont['LC1'].plot_cycles(coords=('v','asd'), figure='fig3', exclude='P2')
#PyCont.plot.fig3.axes1.axes.set_title('Phase Space')
#plt.savefig('Phase_Space.pdf')
#
#PyCont.plot.toggleAll('off', bytype='P')
#PyCont.plot.fig3.refresh()
#plt.legend(loc=2)
show()
