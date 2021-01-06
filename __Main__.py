import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
from random import gauss
from itertools import *
import time

from Networks import Network
from Components.Source import Source
from Components.Heat_exchanger import HEX
from Components.Pipes import Pipe
from Components.Ressources import *


dx = 90  # spatial discretization step (m)
dt = 60  # discretization time step (s)
##Simulation
def simulation( RES, Ts2_h):
    '''Ts2_h:  list of list of temperature per hour (ie scheduled demand for each substation)'''
    nb_SS = len(Ts2_h)
    t = []
    T_in = []
    t_tot = 0
    m_dot = []
    mdot_SS = []
    T_return = []
    T_return_SS = []
    T_supplySS = []
    T_supply_secondary = []
    T_secondary_demand = []
    
    time1 = time.time()
    for j in range(24):
        for p, T_h in enumerate(Ts2_h):
            RES.substations[p][0].Ts2 = T_h[j]
        
        for m in range(60):
            RES.iteration()
            t.append(t_tot)
            T_in.append(RES.supplyT)
            T_return.append(RES.returnT)
            T_return_SS.append([X[0].Tr1 for X in RES.substations])
            T_supplySS.append([X[0].Ts1 for X in RES.substations])
            m_dot.append(RES.m_dot)
            mdot_SS.append([X[0].m_dot1 for X in RES.substations])
            T_secondary_demand.append([X[0].Ts2 for X in RES.substations])
            T_supply_secondary.append([X[0].Ts2_vrai for X in RES.substations])
            t_tot += 1
            
    time2 = time.time()
    plt.figure()
    plt.plot(t, T_return, label = 'Network return temperature at the source')
    for i in range(len(T_return_SS[0])):
        plt.plot(t, [a[i] for a in T_return_SS], label= f'return T substation{i}')
    plt.legend()
    
    plt.figure()

    plt.plot(t, m_dot, label = 'total Network mass flow')
    for i in range(len(mdot_SS[0])):
        plt.plot(t, [a[i] for a in mdot_SS], label= f'mass flow {i}')
    plt.legend()
    
    plt.figure()
    plt.plot(t, T_in, label = 'supply temperature at source')
    for i in range(len(T_supplySS[0])):
        plt.plot(t, [a[i] for a in T_supplySS], label = f'T_supply SS {i}')
    plt.legend()
    
    for i in range(len(T_supply_secondary[0])):
        plt.figure()
        plt.plot(t, [a[i] for a in T_supply_secondary], label = f'T_supply_secondary {i}')
        plt.plot(t, [a[i] for a in T_secondary_demand], label = 'T_secondary_demand')
        plt.plot(t, [a[i] for a in T_supplySS], label = 'T_supply SS network side')
        plt.plot(t, [a[i] for a in T_return_SS], label = 'T_return SS network side')
        plt.legend()     
        plt.figure()
        plt.plot(t, [np.abs(a[i] - b[i]) for a, b in zip(T_supply_secondary,T_secondary_demand) ], label = f'Secondary supply default {i}')
        plt.legend() 
    
    plt.show()
    return time2 - time1
    
    
##Model parameters
#Inlet temperature of the network
Tint = [gauss(65,4) for i in range(24)]
#Temperature demand (profiles)
Ts2_1 = [41, 38, 37, 37, 40, 47, 52, 49, 45, 43, 42, 43, 46, 45, 41, 42, 42, 43, 50, 49, 48, 46, 45, 42]

Ts2_2 = [42, 42.5, 43, 43.5, 44, 48.5, 53, 53, 52, 52, 50, 49, 48,47, 46.5, 46, 49, 49.5, 48, 49, 49.5, 50, 50, 42]
Ts2_2bis = [x - 3 for x in Ts2_2]
Ts2_3 = [41, 38, 37, 37, 40, 44, 48, 47, 40, 39, 39, 38.5, 39, 39.5, 39, 38.5, 39.5, 39, 50, 49.5, 49, 46, 46, 42]

noise = [gauss(0,2) for i in range(24)]
t = [i for i in range(24)]

Tr2 = 30
m_dot2 = 1.6
    
#Model components
SS1 = [HEX(70, 45, 44, 28, 1.6, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
SS2 = [HEX(70, 45, 44, 32, 1.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
SS3 = [HEX(70, 45, 44, 25, 1.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]

SRC1 = Source(70, 20)

NET = Network(75, SRC1, [SS1])
NET2 = Network(75, SRC1, [SS1, SS2])
NET3 = Network(75, SRC1, [SS1, SS2, SS3])

#simulation(NET3, [Ts2_1, Ts2_2, Ts2_3])
#simulation(NET3, [Ts2_1, Ts2_2bis, Ts2_3])