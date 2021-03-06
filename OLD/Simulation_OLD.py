import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
from random import gauss
from itertools import *
import time

from OLD.Network_twopipesSS_OLD import *
from Components.Source import Source
from Components.Heat_exchanger import HEX
from Components.Pipes import Pipe
from Components.Ressources import *


dx = 90  # spatial discretization step (m)
dt = 60  # discretization time step (s)
##Simulation
def initialise(NET):
    '''init = list->tuple(dict{'Tr1' : hex.Tr1, 'Ts2_vrai' : hex.Ts2_vrai, 'Ts1' : hex.Ts1}, pipe1.pipeS_T, pipe2.pipeS_T)''' #En réalité le dict ne sert pas à grand chose puisque ce sont des variables que l'on ajuste lors de l'initialisation
    
    for i, (hex, pipe1, pipe2) in enumerate(NET.substations):
        hex.Tr1, hex.Ts2_vrai, hex.Ts1 = init[i][0]['Tr1'], init[i][0]['Ts2_vrai'], init[i][0]['Ts1']
        pipe1.pipeS_T = list(init[i][1]) #Temperature en Kelvin!
        pipe2.pipeS_T = list(init[i][2])

def simulation(RES, Ts2_h):
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
    P_source = []
    P_boiler = []
    
    K = []
    UA = []
    
    #SIMULATION
    # # First loop for initialisation
    # if not RES.alreadyrun:
    #     RES.alreadyrun = True
    #     for j in range(24):
    #         for p, T_h in enumerate(Ts2_h):
    #             RES.substations[p][0].Ts2 = T_h[j]
    #         
    #         for m in range(60):
    #             RES.iteration()
            
    time1 = time.time()
    for j in range(24):
        for p, T_h in enumerate(Ts2_h):
            RES.substations[p][0].Ts2 = T_h[j]
        
        for m in range(60):
            #print(t_tot)
            try :
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
                P_source.append(Cp*RES.src.m_dot*(RES.src.Ts_Geo - RES.src.Tr_Geo))
                if RES.NETtype == 'boiler':
                    P_boiler.append(RES.boiler * RES.P_boiler)
                t_tot += 1
                
                # UA.append([X[0].UA() for X in RES.substations])
                # K.append([X[0].K for X in RES.substations])
                 
            except ValueError:
                plt.figure()
                plt.title('UA')
                for i in range(len(UA[0])):
                    plt.plot(t, [a[i] for a in UA])
                
                plt.show()
                raise ValueError
            
    time2 = time.time()
    
    #PLOT
    t = [x/60 for x in t]
    
    plt.figure()
    plt.title('Return Temperature (°C)')
    plt.plot(t, T_return, label = 'Network return temperature')
    for i in range(len(T_return_SS[0])):
        plt.plot(t, [a[i] for a in T_return_SS], label= f'Return_T_net SS_{i}')
    plt.legend()
    
    plt.figure()
    plt.title('Mass flow (kg/s)')
    plt.plot(t, m_dot, label = 'Total network mass flow')
    for i in range(len(mdot_SS[0])):
        plt.plot(t, [a[i] for a in mdot_SS], label= f'Massflow SS_{i}')
    plt.legend()
    
    plt.figure()
    plt.title('Supply Temperature (°C)')
    plt.plot(t, T_in, label = 'Network supply temperature')
    for i in range(len(T_supplySS[0])):
        plt.plot(t, [a[i] for a in T_supplySS], label = f'Supply_T_net SS_{i}')
    plt.legend()
    
    for i in range(len(T_supply_secondary[0])):
        plt.figure()
        plt.title(f'SS_{i}')
        plt.plot(t, [a[i] for a in T_secondary_demand], label = 'Demand Ts2')
        plt.plot(t, [a[i] for a in T_supply_secondary], label = 'Supplied Ts2')
        
        plt.plot(t, [a[i] for a in T_supplySS], label = 'Supply_T_net SS')
        plt.plot(t, [a[i] for a in T_return_SS], label = 'Return_T_net SS')
        plt.legend()   
    
    plt.figure()
    plt.title('Secondary supply default (°C)')
    for i in range(len(T_supply_secondary[0])):
        plt.plot(t, [np.abs(a[i] - b[i]) for a, b in zip(T_supply_secondary,T_secondary_demand) ], label = f'Secondary supply default SS_{i}')
    plt.legend() 
    
    plt.figure()
    plt.title('Source Power (kW)')
    plt.plot(t, [x/1000 for x in P_source], label = 'Geothermal Source' )
    if RES.NETtype == 'boiler':
        plt.plot(t, [x/1000 for x in P_boiler], label = 'Gas boiler')
    plt.legend() 
    
    # plt.figure()
    # plt.title('UA')
    # for i in range(len(UA[0])):
    #     plt.plot(t, [a[i] for a in UA], label = f'UA_{i}')
    # plt.legend()
    # 
    # plt.figure()
    # plt.title('K')
    # for i in range(len(K[0])):
    #     plt.plot(t, [a[i] for a in K], label = f'K_{i}')
        
    plt.legend()
    plt.show()
    return time2 - time1
    
    
##Model parameters
#Inlet temperature of the network
Tint = [gauss(65,4) for i in range(24)]
#Temperature demand (profiles)
Ts2_1 = [41, 38, 37, 37, 40, 47, 52, 49, 45, 43, 42, 43, 46, 45, 41, 42, 42, 43, 50, 49, 48, 46, 45, 42]

Ts2_2 = [42, 42.5, 43, 43.5, 44, 48.5, 53, 53, 52, 52, 50, 49, 48,47, 46.5, 46, 49, 49.5, 48, 49, 49.5, 50, 50, 42]
Ts2_2bis = [x - 5 for x in Ts2_2]
Ts2_3 = [41, 38, 37, 37, 40, 44, 48, 47, 40, 39, 39, 38.5, 39, 39.5, 39, 38.5, 39.5, 39, 50, 49.5, 49, 46, 46, 42]

noise = [gauss(0,2) for i in range(24)]
t = [i for i in range(24)]

Tr2 = 30
m_dot2 = 1.6
    
#Model components
SS1 = [HEX(70, 45, 44, 28, 1.6, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
SS2 = [HEX(70, 45, 44, 32, 1.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
SS3 = [HEX(70, 45, 44, 25, 1.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]

SRC1 = Source(70, 20,150000, 3)

NET = Network(SRC1, [SS1])
NET2 = Network(SRC1, [SS1, SS2])
NET3 = Network(SRC1, [SS1, SS2, SS3])

NETb = Network_boiler(SRC1, 265000/3, [SS1])
NETb3 = Network_boiler(SRC1, 265000, [SS1, SS2, SS3])
NETb3bis = Network_boiler(SRC1, 300000, [SS1, SS2, SS3])

#simulation(NET3, [Ts2_1, Ts2_2, Ts2_3])
#simulation(NET3, [Ts2_1, Ts2_2bis, Ts2_3])


        

def test(hex):
    a = 0.3275
    Tr1, Ts2_vrai, Ts2, mdot = hex.Tr1, hex.Ts2_vrai,hex.Ts2, hex.m_dot1
    Q = Cp * hex.m_dot2 * (Ts2 - hex.Tr2)
    UA = hex.UA()
    print(UA)
    TS1 = hex.Ts1
    A = []
    b = []
    print(UA)
    
    l = np.linspace(Ts2 + 10, Ts2 + 100)
    #print(l)
    for Ts1 in l:
        A.append(2 *(Q/(UA))**a - (Ts1 - Ts2)**a)
        b.append(Ts1 - Ts2)
    plt.figure()
    plt.plot(b, A)    

    c = []
    d = []
    for Ts2 in np.linspace(Ts2_vrai -5, Ts2_vrai + 15 ):
        c.append(2 * (Q/(UA))**a - (TS1 - Ts2)**a)
        d.append(TS1 - Ts2)
    plt.figure()
    plt.plot(d, c)    

    plt.show()