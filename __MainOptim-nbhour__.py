import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
from random import gauss
from itertools import *
import time

from geneticalgorithm import geneticalgorithm as ga
from Components.Data_meteo_reading import Ts2_A, Ts2_B, Ts2_C
from Components.Data_meteo_reading import Ta as ext_T
from Networks import Network_bOptim
from Components.Source import Source
from Components.Heat_exchanger import HEX, HEX_nom
from Components.Pipes import Pipe
from Components.Ressources import *


##Simulation
class Simulation():

    def __init__(self, RES, Ts2_h):
        '''Ts2_h:  list of list of temperature per hour (ie scheduled demand for each substation)'''
        self._Ts2_h = Ts2_h  #Ts2_h is defined as a property so we can control access to it and therfore simulation lenght
        nbh = len(Ts2_h[0])
        for SS_instruct in Ts2_h: #verifying that each substation received the same number of instructions
            if len(SS_instruct) != nbh:
                raise ValueError('Every substation should receive the same number of instruction')
        self.nb_hour = nbh
        self.RES = RES
        self.nb_SS = len(RES.substations)
        
        self.t = list(range(0, self.nb_hour*(3600//dt)))
        self.E_Geo = 0
        self.E_boiler = 0
        self.E_default = 0
        
        self.cost_Tdefault_SS = [0] * self.nb_SS
        self.cost_constraintT = 0
        
        self.initialised = False
        self.save_init_ = []

    def _get_Ts2_h(self):
        return self._Ts2_h 
        
    def _set_Ts2_h(self, Ts2_h):
        '''when user wants to modify Ts2_h, this method is called and consequently modifies the duration of simulation'''
        if len(Ts2_h) != self.nb_SS: 
            raise ValueError('Instruction should be given for the corresponding number of substations in the network')
        else:
            nbh = len(Ts2_h[0])
            for SS_instruct in Ts2_h: #verifying that each substation received the same number of instructions
                if len(SS_instruct) != nbh:
                    raise ValueError('Every substation should receive the same number of instruction')
                    
            self._Ts2_h = Ts2_h
            self.nb_hour = nbh
            self.t = list(range(0, self.nb_hour*(3600//dt)))
            
    Ts2_h = property(_get_Ts2_h, _set_Ts2_h)  #Definition of the Ts2_h attribut of Simulation class

    def initialisation(self):
        # First loop for initialisation
        if not self.initialised:
            t_tot = 0
            self.initialised = True
            self.RES.alreadyrun = True
            for j in range(self.nb_hour):
                for p, T_h in enumerate(self.Ts2_h):
                    self.RES.substations[p][0].Ts2 = T_h[j]
                for m in range(3600//dt):
                    self.RES.iteration()
                    t_tot += 1
                    
            self.save_init_ = [self.RES.supplyT, self.RES.m_dot, self.RES.returnT, self.RES.Ts_nodes, self.RES.Tr_nodes] #Initial conditions
            save =  []
            for i, (hex, pipe1, pipe2) in enumerate(self.RES.substations):
                save.append(({'Tr1' : hex.Tr1, 'Ts2_vrai' : hex.Ts2_vrai, 'Ts1' : hex.Ts1, 'm1':hex.m_dot1}, np.copy(pipe1.pipeS_T), np.copy(pipe2.pipeS_T), np.copy(pipe1.pipeR_T), np.copy(pipe2.pipeR_T) ))
            self.save_init_.append(save)
            
        else:
            self.RES.supplyT = self.save_init_[0]
            self.RES.m_dot = self.save_init_[1]
            self.RES.returnT = self.save_init_[2]
            self.RES.Ts_nodes = self.save_init_[3]
            self.RES.Tr_nodes = self.save_init_[4]
            save = self.save_init_[5]
            
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
        self.cost_constraintT = 0
        
        time1 = time.time()
        t_tot = 0
        i_instruct = 0
        
        t_unsupplied = [0] * self.nb_SS
        
        for j in range(self.nb_hour):
            for p, T_h in enumerate(self.Ts2_h):
                self.RES.substations[p][0].Ts2 = T_h[j]
            
            for m in range((3600//dt)):
                if t_tot % nbiter_instruc == 0:
                    self.RES.P_boiler = P_boiler_instruction[i_instruct]
                    i_instruct += 1
                    
                self.RES.iteration()

                self.E_Geo += self.RES.P_Geo
                self.E_boiler += self.RES.P_boiler
                self.E_default += self.RES.P_supplied - self.RES.P_demand
                maxT = self.RES.maxT
                
                if maxT > 95:
                    self.cost_constraintT += 3*(np.exp(0.2*(maxT-105)) - np.exp(-2)) #+ np.exp(100 * (maxT - 95))
                
                #Calculation of the cost of supply default at time t
                for i, Tdefault in enumerate(self.RES.Tsupply_default_SS):
                    if Tdefault <= 10**(-5):
                        t_unsupplied[i] = t_tot
                    else:
                        #print(i, len(t_unsupplied))
                        ct = 2*(np.exp(((t_tot - t_unsupplied[i])*dt/60 -25)/10) - np.exp(-2.5))
                        cT = 2*(np.exp((Tdefault -1)/2) - np.exp(-1/2))
                        
                        #ct = 1/(1 + np.exp((-2/10)*((t_tot - t_unsupplied[i]) -15)))
                        #cT = 1/(1 + np.exp((-3/2)*(Tdefault -1.5)))
                    
                        self.cost_Tdefault_SS[i] += ct*cT
                
                t_tot+= 1
                
        time2 = time.time()
    
    def objective_function(self, P_boiler_instruction, exe_time = False):
        time1 = time.time()

        #self.simulation(P_boiler_instruction)
        try:
            self.simulation(P_boiler_instruction)
        except ValueError as e:
            return float(e.__str__())*1e200
            
        #kilowatt-hour
        E_boiler = self.E_boiler*dt/3600
        E_Geo = self.E_Geo*dt/3600
        E_tot = E_boiler + E_Geo
        E_default = self.E_default*dt/3600
        
        C1 = E_boiler/(E_tot) #* C_kWh * E_boiler
        
        C2 = sum(self.cost_Tdefault_SS)/(len(self.t)*self.nb_SS)
        
        C_constraint = self.cost_constraintT/(len(self.t)*self.nb_SS)
        
        #print(C1, C2, C_constraint)
        
        time2 = time.time()
        if exe_time:
            print(time2-time1)
            
        return C1 + C2 + C_constraint
        
    
    def plot(self, P_boiler_instruction):
        self.initialisation()
        nbiter_instruc = len(self.t)//len(P_boiler_instruction)
        
        T_return_SS = []
        T_supplySS = []
        T_supply_secondary = []
        T_secondary_demand = []
        P_boiler = []
        water_flow = []
        
        time1 = time.time()
        i_instruct = 0
        t_tot = 0
        for j in range(self.nb_hour):
            for p, T_h in enumerate(self.Ts2_h):
                self.RES.substations[p][0].Ts2 = T_h[j]
            
            for m in range(3600//dt):
                if t_tot % nbiter_instruc == 0:
                    self.RES.P_boiler = P_boiler_instruction[i_instruct]
                    i_instruct += 1
                
                self.RES.iteration()
                
                P_boiler.append(self.RES.P_boiler)
                T_return_SS.append([X[0].Tr1 for X in self.RES.substations])
                T_supplySS.append([X[0].Ts1 for X in self.RES.substations])
                T_secondary_demand.append([X[0].Ts2 for X in self.RES.substations])
                T_supply_secondary.append([X[0].Ts2_vrai for X in self.RES.substations])
                water_flow.append(self.RES.m_dot)
                
                t_tot += 1
                    
        time2 = time.time()
        if len(self.t)//(3600//dt) > 72:
            t = [x/(3600*24)*dt for x in self.t]
        else:
            t = [x*dt/3600 for x in self.t]
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
        plt.title(f'Boiler Power')
        plt.plot(t, [a for a in P_boiler])
        plt.show()
        
        plt.figure()
        plt.title(f'Network water flow')
        plt.plot(t, [a for a in water_flow])
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
#
SS1 = [HEX(70, 45, 44, 28, 1.6, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
SS1_bis = [HEX(70, 65, 44, 40, 1.6, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
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
Bex = [x*0.4 for x in Boilerexemple]

NET = Network_bOptim(SRC1, 0, [SS1])
NET_bis = Network_bOptim(SRC1, 0, [SS1_bis])
NET3 = Network_bOptim(SRC1, 0, [SS1, SS2, SS3])

#simulation(NET3, [Ts2_1, Ts2_2, Ts2_3])
#simulation(NET3, [Ts2_1, Ts2_2bis, Ts2_3])


##OPTIM
A1 = Simulation(NET, [Ts2_2])
A1_bis = Simulation(NET, [[x+ 20 for x in Ts2_1]])
A3 = Simulation(NET3,[Ts2_1, Ts2_2, Ts2_3])

algorithm_param={'max_num_iteration': 100, 'population_size': 75, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}

def optim(dim, MOD, max_P = 90000):
    time1 = time.time()
    #max_Powers = [0]*dim   #Si dim == len(Ta) ok! sinon non mais ça n'aurait pas d'interet
    
    # Ta_ind = list(enumerate(MOD.Ta))
    # sortedTa = sorted(Ta_ind, key = lambda x: x[1])
    # i_Ta = [x[0] for x in sortedTa[::-1]]
    # a_max = 0
    # for i in i_Ta:
    #     P = [0]*dim
    #     for a in range(a_max, max_P, 10000):
    #         P[i] = a
    #         if MOD.objective_function(P) < 1e199:
    #             a_max = a
    #         else:
    #             break
    #     
    #     for a in range(a_max, max_P, 1000):
    #         P[i] = a
    #         if MOD.objective_function(P) < 1e199:
    #             a_max = a
    #         else:
    #             break
    #             
    #     for a in range(a_max, max_P, 100):
    #         P[i] = a
    #         if MOD.objective_function(P) < 1e199:
    #             a_max = a
    #         else:
    #             break   
    #             
    #     max_Powers[i]= a_max
    #     print(max_Powers)
    # 
    # print(max_Powers)
    # varbound=np.array([[0,max_P] for max_P in max_Powers])
    
    print(max_P*MOD.nb_SS)
    P = max_P*MOD.nb_SS
    
    if len(MOD.PBoiler_Ta) ==0:
        a = MOD.objective_function_Ta()
    n = len(MOD.PBoiler_Ta)
    step = n//dim
    varbound=np.array([[x/2, 2*x] for x in MOD.PBoiler_Ta[step-1::step]])
    #varbound=np.array([[0,100000]]*dim)
    
    
    time2 = time.time()
    #print(f'temps etape 1: {time2- time1}')
    
    model=ga(function=MOD.objective_function,dimension= dim,variable_type='int',variable_boundaries=varbound, algorithm_parameters=algorithm_param)

    t1 = time.time()
    model.run()
    t2 = time.time()
    print(f'time = {t2-t1}')
    
    Boiler_instructP = list(model.output_dict['variable'])
    
    MOD.plot(Boiler_instructP)
    
    return Boiler_instructP
    

## Ts1 linear function of Ta

class Simulation_Ta(Simulation):

    def __init__(self, RES, Ta):
        '''Ta:  list of external temperature per hour'''
        Ts2_h = []
        for p, SS in enumerate(RES.substations):
            Ts2_h.append(SS[0].f_Ts2(Ta))
        Simulation.__init__(self, RES, Ts2_h)
        self._Ta = Ta
        self.f_Ts1 = lambda x: 70
        self.PBoiler_Ta = [] #This attribute allows us to directly compare the instructions given by the linear law and by the optimisation process, as well as reducing the dimension of the problem for optimisation.
        
    def _get_Ta(self):
        return self._Ta
        
    def _set_Ta(self, Ta):
        '''when user wants to modify Ta, this method is called and consequently modifies the duration of simulation'''
        Ts2_h = []
        for p, SS in enumerate(self.RES.substations):
            Ts2_h.append(SS[0].f_Ts2(Ta))
        self.Ts2_h = Ts2_h

    Ta = property(_get_Ta, _set_Ta)
            
    def refined_Ts1(self):
        upper_limit = 90
        lower_limit = 60
        optimf_Ts1 = lambda x: 70
        opt_lT, opt_uT = 70, 70
        min_score = self.objective_function_Ta(optimf_Ts1)
        for uT in np.linspace(lower_limit, upper_limit, 30):
            for lT in np.linspace(lower_limit, uT, (uT -lower_limit)+1):
                f = lambda Ta: uT - ((uT - lT)/27)* (Ta + 7)
                try:
                    f_score= self.objective_function_Ta(f)

                except ValueError:
                    #print(f'Error, upperT = {uT}, lowerT = {lT} ')
                    f_score = min_score+1
                if f_score < min_score:
                    #print(f_score, min_score)
                    opt_lT, opt_uT = lT, uT
                    min_score = f_score
        
        optimf_Ts1 = lambda Ta: opt_uT - ((opt_uT - opt_lT)/27)* (Ta + 7)
        self.f_Ts1 = optimf_Ts1
        print(f'optim score = {self.objective_function_Ta()}')
        
        
    def simulation_Ta(self, f_Ts1 = None, nbiter_per_instruc = 1):
        '''f_Ts1: linear rule between ambiant Temperature and instruction Ts1, if None it is the function self.f_Ts1
            nbiter_per_instruc : defines the precision/time step of instruction
        '''
            
        self.initialisation()
        if f_Ts1 == None:
            f_Ts1 = self.f_Ts1
            
        self.E_Geo = 0
        self.E_boiler = 0
        self.E_default = 0
        self.cost_Tdefault_SS = [0] * self.nb_SS
        self.cost_constraintT = 0
        
        time1 = time.time()
        t_tot = 0
        t_unsupplied = [0] * self.nb_SS
        
        for j in range(self.nb_hour):
            for p, T_h in enumerate(self.Ts2_h):
                self.RES.substations[p][0].Ts2 = T_h[j]
            
            #Calculation of wanted Ts1 with linear function
            instructTs1 = f_Ts1(self.Ta[j])
            
            for m in range((3600//dt)):
                Ts1 = self.RES.supplyT
                mdot = self.RES.m_dot
                
                if t_tot % nbiter_per_instruc == 0:
                    #Calculation of Pboiler required:
                    if instructTs1 > Ts1:
                        P = (instructTs1 - Ts1)*Cp*mdot 
                    else: 
                        P = 0
                    self.RES.P_boiler = P
                    self.PBoiler_Ta.append(P)
                    
                self.RES.iteration()

                self.E_Geo += self.RES.P_Geo
                self.E_boiler += self.RES.P_boiler
                self.E_default += self.RES.P_supplied - self.RES.P_demand
                maxT = self.RES.maxT
                
                if maxT > 85:
                    self.cost_constraintT += 3*(np.exp(0.2*(maxT-95)) - np.exp(-2)) #+ np.exp(100 * (maxT - 95))
                    
                #Calculation of the cost of supply default at time t
                for i, Tdefault in enumerate(self.RES.Tsupply_default_SS):
                    if Tdefault <= 10**(-5):
                        t_unsupplied[i] = t_tot
                    else:
                        ct = 2*(np.exp(((t_tot - t_unsupplied[i]) -25)/10) - np.exp(-2.5))
                        cT = 2*(np.exp((Tdefault -1)/2) - np.exp(-1/2))

                        self.cost_Tdefault_SS[i] += ct*cT
                
                t_tot+= 1
                
        time2 = time.time()
        
        
    def objective_function_Ta(self, f_Ts1 = None, dim_perday = (1440*60)//dt):
        #Reinitialisation of self.PBoiler_Ta
        self.PBoiler_Ta = []
        time1 = time.time()
        
        nb_instruct = dim_perday * self.nb_hour//24
        if nb_instruct  > len(self.t):
            raise ValueError("Instructions can't be more precise than time step of pipe model")
        else:
            nbiter_per_instruc = len(self.t)//nb_instruct
        
        self.simulation_Ta(f_Ts1, nbiter_per_instruc)
        
        
        #kilowatt-hour
        E_boiler = self.E_boiler*dt/3600
        E_Geo = self.E_Geo*dt/3600
        E_tot = E_boiler + E_Geo
        E_default = self.E_default*dt/3600
        
        C1 = E_boiler/(E_tot) #* C_kWh * E_boiler
        
        C2 = sum(self.cost_Tdefault_SS)/(len(self.t)*self.nb_SS)
        
        C_constraint = self.cost_constraintT/(len(self.t)*self.nb_SS)
        
        #print(C1, C2, C_constraint)
        time2 = time.time()
        #print(time2-time1)
        return C1 + C2 + C_constraint
    
    
    def plot_Ta(self):
        self.initialisation()
        
        T_return_SS = []
        T_supplySS = []
        T_supply_secondary = []
        T_secondary_demand = []
        P_boiler = []
        
        time1 = time.time()
        t_tot = 0
        for j in range(self.nb_hour):
            for p, T_h in enumerate(self.Ts2_h):
                self.RES.substations[p][0].Ts2 = T_h[j]
            
            instructTs1 = self.f_Ts1(self.Ta[j])
            for m in range((3600//dt)):
                Ts1 = self.RES.supplyT
                mdot = self.RES.m_dot
                #Calculation of Pboiler required:
                if instructTs1 > Ts1:
                    P = (instructTs1 - Ts1)*Cp*mdot 
                else: 
                    P = 0
                    
                self.RES.P_boiler = P
                self.RES.iteration()
                
                P_boiler.append(self.RES.P_boiler)
                T_return_SS.append([X[0].Tr1 for X in self.RES.substations])
                T_supplySS.append([X[0].Ts1 for X in self.RES.substations])
                T_secondary_demand.append([X[0].Ts2 for X in self.RES.substations])
                T_supply_secondary.append([X[0].Ts2_vrai for X in self.RES.substations])
                
                t_tot += 1
                    
        time2 = time.time()
        
        if len(self.t)//(3600//dt) > 72:
            t = [x*dt/(3600*24) for x in self.t]
        else:
            t = [x*dt/3600 for x in self.t]
            
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
        plt.title(f'Boiler Power')
        plt.plot(t, [a for a in P_boiler])
        plt.legend()
        
        plt.show()
    
        
##Models
SS1_nom = [HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
SS2_nom = [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 75)]
SS3_nom = [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
SS4_nom = [HEX_nom(100e3, 43, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 700), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 20)]
NET_nom = Network_bOptim(SRC1, 0, [SS1_nom])
NET3_nom = Network_bOptim(SRC1, 0, [SS1_nom, SS2_nom, SS3_nom])


Ta = [14.212499999999999, 14.212499999999999, 13.837499999999999, 13.387500000000001, 13.075500000000002, 12.637500000000001, 12.1875, 11.7375, 11.4255, 10.9875, 12.874500000000001, 14.562000000000001, 16.212, 18.138, 19.5375, 19.512, 18.674999999999997, 17.5005, 16.987499999999997, 16.5375, 16.2255, 15.787500000000001, 15.337499999999999, 15.0255]
A_nom = Simulation_Ta(NET_nom, Ta)
A3_nom = Simulation_Ta(NET3_nom, Ta) #[Ts2_Ta])

SRC2 = Source(75, 20)
NET4_nom = Network_bOptim(SRC2, 0, [SS1_nom, SS2_nom, SS3_nom, SS4_nom])
A4_nom = Simulation_Ta(NET4_nom, ext_T)

 # [160810.  65242. 174214. 130763.  56980. 116297. 198465. 193662. 199989.
 # 198039. 199616. 199754. 191696.  81751.  33179.  17442.  48478.  81413.
 # 140092. 184138. 135619.  79352. 197408. 194840.]

 #   Objective function:
 # 1.660859759466892e+30
 
BP = [174214,  174214, 174214, 174214,  174214, 174214, 198465, 193662, 199989, 198039, 199616, 199754,191696,  81751,  81751,  81751,  81751,  81413, 140092, 184138, 135619 , 81751, 197408, 194840]


## GA-Parallel

# A4_nom_1 = Simulation_Ta(Network_bOptim(Source(75, 20), 0, [[HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 75)], [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 43, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 700), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 20)]]), ext_T)
# 
# A4_nom_2 = Simulation_Ta(Network_bOptim(Source(75, 20), 0, [[HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 75)], [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 43, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 700), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 20)]]), ext_T)
# 
# A4_nom_3 = Simulation_Ta(Network_bOptim(Source(75, 20), 0, [[HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 75)], [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 43, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 700), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 20)]]), ext_T)
# 
# A4_nom_4 = Simulation_Ta(Network_bOptim(Source(75, 20), 0, [[HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 75)], [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 43, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 700), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 20)]]), ext_T)
# 
# A4_nom_5 = Simulation_Ta(Network_bOptim(Source(75, 20), 0, [[HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 75)], [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)], [HEX_nom(100e3, 43, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 700), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 20)]]), ext_T)
# 
# #algorithm_param={'max_num_iteration': 100, 'population_size': 75, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}
# 
# MOD = [A4_nom_1, A4_nom_2, A4_nom_3, A4_nom_4, A4_nom_5]
# def optim_A4(dim, MOD, max_P = 90000):
#     varbound=np.array([[0,200000]]*dim)
#     
#     model=ga(function=[A.objective_function for A in MOD],dimension= dim,variable_type='int',variable_boundaries=varbound, algorithm_parameters=algorithm_param)
# 
#     t1 = time.time()
#     model.run()
#     t2 = time.time()
#     print(f'time = {t2-t1}')
#     
#     Boiler_instructP = list(model.output_dict['variable'])
#     
#     MOD[0].plot(Boiler_instructP)
#     
#     return Boiler_instructP

        













