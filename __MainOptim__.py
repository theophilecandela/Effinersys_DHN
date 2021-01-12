import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
from random import gauss
from itertools import *
import time

from geneticalgorithm import geneticalgorithm as ga

from Networks import Network_bOptim
from Components.Source import Source
from Components.Heat_exchanger import HEX
from Components.Pipes import Pipe
from Components.Ressources import *


dx = 90  # spatial discretization step (m)
dt = 60  # discretization time step (s)
##Simulation
class Simulation():

    def __init__(self, RES, Ts2_h):
        '''Ts2_h:  list of list of temperature per hour (ie scheduled demand for each substation)'''
        self.Ts2_h = Ts2_h
        self.RES = RES
        self.nb_SS = len(Ts2_h)
        
        self.t = []
        self.E_Geo = 0
        self.E_boiler = 0
        self.E_default = 0
        
        self.cost_Tdefault_SS = [0] * self.nb_SS
        
        self.initialised = False
        self.save_init_ = []

    def initialisation(self):
        # First loop for initialisation
        if not self.initialised:
            t_tot = 0
            self.initialised = True
            self.RES.alreadyrun = True
            for j in range(24):
                for p, T_h in enumerate(self.Ts2_h):
                    self.RES.substations[p][0].Ts2 = T_h[j]
                for m in range(60):
                    self.RES.iteration()
                    self.t.append(t_tot)
                    t_tot += 1
                    
            self.save_init_ = [self.RES.supplyT, self.RES.m_dot, self.RES.returnT] #Initial conditions
            save =  []
            for i, (hex, pipe1, pipe2) in enumerate(self.RES.substations):
                save.append(({'Tr1' : hex.Tr1, 'Ts2_vrai' : hex.Ts2_vrai, 'Ts1' : hex.Ts1, 'm1':hex.m_dot1}, np.copy(pipe1.pipeS_T), np.copy(pipe2.pipeS_T), np.copy(pipe1.pipeR_T), np.copy(pipe2.pipeR_T) ))
            self.save_init_.append(save)
            
        else:
            self.RES.supplyT = self.save_init_[0]
            self.RES.m_dot = self.save_init_[1]
            self.RES.returnT = self.save_init_[2]
            save = self.save_init_[3]
            
            for i, (hex, pipe1, pipe2) in enumerate(self.RES.substations):
                hex.Tr1, hex.Ts2_vrai, hex.Ts1, hex.m_dot1 = save[i][0]['Tr1'], save[i][0]['Ts2_vrai'], save[i][0]['Ts1'], save[i][0]['m1']
                pipe1.pipeS_T = list(save[i][1])
                pipe2.pipeS_T = list(save[i][2])
                pipe1.pipeR_T = list(save[i][3])
                pipe2.pipeR_T = list(save[i][4])
                
                
    def simulation(self, P_boiler_instruction):
        self.initialisation()
        
        nbiter_instruc = len(self.t)//len(P_boiler_instruction)
        self.E_Geo = 0
        self.E_boiler = 0
        self.E_default = 0
        self.cost_Tdefault_SS = [0] * self.nb_SS
        
        time1 = time.time()
        t_tot = 0
        i_instruct = 0
        
        t_unsupplied = [0] * self.nb_SS
        
        for j in range(24):
            for p, T_h in enumerate(self.Ts2_h):
                self.RES.substations[p][0].Ts2 = T_h[j]
            
            for m in range(60):
                if t_tot % nbiter_instruc == 0:
                    self.RES.P_boiler = P_boiler_instruction[i_instruct]
                    i_instruct += 1
                    
                self.RES.iteration()

                self.E_Geo += self.RES.P_Geo
                self.E_boiler += self.RES.P_boiler
                self.E_default += self.RES.P_supplied - self.RES.P_demand
                
                #Calculation of the cost of supply default at time t
                for i, Tdefault in enumerate(self.RES.Tsupply_default_SS):
                    if Tdefault <= 10**(-5):
                        t_unsupplied[i] = t_tot
                    else:
                        ct = 1/(1 + np.exp((-2/10)*((t_tot - t_unsupplied[i]) -15)))
                        cT = 1/(1 + np.exp((-3/2)*(Tdefault -1.5)))
                    
                        self.cost_Tdefault_SS[i] += ct*cT
                
                t_tot+= 1
                
        time2 = time.time()
    
    def objective_function(self, P_boiler_instruction):
        time1 = time.time()
        self.simulation(P_boiler_instruction)
        
        #kilowatt-hour
        E_boiler = self.E_boiler/60
        E_Geo = self.E_Geo/60
        E_tot = E_boiler + E_Geo
        E_default = self.E_default/60
        
        C1 = E_boiler/(E_tot) #* C_kWh * E_boiler
        
        C2 = sum(self.cost_Tdefault_SS)/(len(self.t)*self.nb_SS)

        #print(C1, C2)
        
        time2 = time.time()
        #print(time2-time1)
        return C1 + C2
        
    
    def plot(self, P_boiler_instruction):
        self.initialisation()
        nbiter_instruc = len(self.t)//len(P_boiler_instruction)
        
        T_return_SS = []
        T_supplySS = []
        T_supply_secondary = []
        T_secondary_demand = []
        
        time1 = time.time()
        i_instruct = 0
        t_tot = 0
        for j in range(24):
            for p, T_h in enumerate(self.Ts2_h):
                self.RES.substations[p][0].Ts2 = T_h[j]
            
            for m in range(60):
                if t_tot % nbiter_instruc == 0:
                    self.RES.P_boiler = P_boiler_instruction[i_instruct]
                    i_instruct += 1
                
                self.RES.iteration()
                
                T_return_SS.append([X[0].Tr1 for X in self.RES.substations])
                T_supplySS.append([X[0].Ts1 for X in self.RES.substations])
                T_secondary_demand.append([X[0].Ts2 for X in self.RES.substations])
                T_supply_secondary.append([X[0].Ts2_vrai for X in self.RES.substations])
                
                t_tot += 1
                    
        time2 = time.time()
        
        t = [x/60 for x in self.t]
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
        
        plt.show()
        
        
    
    
    
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

SRC1 = Source(70, 20)

# NET = Network(SRC1, [SS1])
# NET2 = Network(SRC1, [SS1, SS2])
# NET3 = Network(SRC1, [SS1, SS2, SS3])
# 
# NETb = Network_boiler(SRC1, 265000/3, [SS1])
# NETb3 = Network_boiler(SRC1, 265000, [SS1, SS2, SS3])
# NETb3bis = Network_boiler(SRC1, 300000, [SS1, SS2, SS3])

Boilerexemple = [88333 * (x - 42)/(53-42) for x in Ts2_2]
Bex = [x*1.4 for x in Boilerexemple]

NET = Network_bOptim(SRC1, 0, [SS1])
NET3 = Network_bOptim(SRC1, 0, [SS1, SS2, SS3])

#simulation(NET3, [Ts2_1, Ts2_2, Ts2_3])
#simulation(NET3, [Ts2_1, Ts2_2bis, Ts2_3])

##
A1 = Simulation(NET, [Ts2_2])
A3 = Simulation(NET3,[Ts2_1, Ts2_2, Ts2_3])

algorithm_param={'max_num_iteration': 50, 'population_size': 50, 'mutation_probability': 0.1, 'elit_ratio': 0.01, 'crossover_probability': 0.5, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': 15}

varbound=np.array([[0,90000]]*24)
model=ga(function=A1.objective_function,dimension=24,variable_type='int',variable_boundaries=varbound, algorithm_parameters=algorithm_param)

# varbound=np.array([[90000,90001]]*3)
# model=ga(function=A1.objective_function,dimension=3,variable_type='int',variable_boundaries=varbound, algorithm_parameters=algorithm_param)




