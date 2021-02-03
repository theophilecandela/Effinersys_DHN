import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
from random import gauss
from itertools import *
import time
import sys

from Components.Ressources import *


##Simulation
class Simulation():
    def __init__(self, RES, Ta):
        '''Ta:  list of external temperature per hour'''

        self._Ta = Ta
        self.f_Ts1 = lambda x: 70
        
        self.nb_hour =len(Ta)
        self.RES = RES
        self.nb_SS = len(RES.substations)
        
        self.t = list(range(0, self.nb_hour*(3600//dt)))
        
        self.E_Geo = 0
        self.E_boiler = 0
        self.E_default = 0
        self.P_boiler = []
         
        self.cost_Tdefault_SS = [0] * self.nb_SS
        self.cost_constraintT = 0
        
        self.storage = self.RES.Storage is not None
        
        self.initialised = False
        self.save_init_ = []
        
    def _get_Ta(self):
        return self._Ta
        
    def _set_Ta(self, Ta):
        '''when user wants to modify Ta, this method is called and consequently modifies the duration of simulation'''
        self.nb_hour =len(Ta)
        self.t = list(range(0, self.nb_hour*(3600//dt)))

    Ta = property(_get_Ta, _set_Ta)
    

    def initialisation(self, T_instruct = None, Storage_instruct = None):
        '''T_instruct: instruct T_supply, step 1hour'''
        # First loop for initialisation
        if not self.initialised:
            t_tot = 0
            
            if T_instruct is not None:
                if len(T_instruct) != self.nb_hour:
                    raise IndexError('There shall be one instruction for each hour')
            
            if Storage_instruct is not None:
                if not self.storage:
                    raise ValueError('There is no storage buffer in this DH Network')
                elif len(Storage_instruct) < self.nb_hour:
                    raise ValueError('There shall be at least on instruction for each hour')
                nbiter_instruc = len(self.t)//len(Storage_instruct)
                i_instruct = 0
            
            for j in range(self.nb_hour):
                for p in range(self.nb_SS):
                    self.RES.substations[p][0].Ts2 = self.RES.substations[p][0].f_Ts2(self.Ta[j])
                if T_instruct is not None:
                    self.RES.Boiler_Tinstruct = T_instruct[j]
                    
                for m in range(3600//dt):
                    if Storage_instruct is not None:
                        if t_tot % nbiter_instruc == 0:
                            self.RES.storage_flow = Storage_instruct[i_instruct]
                            i_instruct += 1
                            
                    self.RES.iteration()
                    t_tot += 1
                    
            self.save_init_ = [self.RES.supplyT, self.RES.m_dot, self.RES.returnT, self.RES.Ts_nodes.copy(), self.RES.Tr_nodes.copy()] #Initial conditions
            save =  []
            for i, (hex, pipe1) in enumerate(self.RES.substations):
                save.append(({'Tr1' : hex.Tr1, 'Ts2_vrai' : hex.Ts2_vrai, 'Ts1' : hex.Ts1, 'm1':hex.m_dot1}, np.copy(pipe1.pipeS_T), np.copy(pipe1.pipeR_T) ))
            self.save_init_.append(save)
            
            if self.storage:
                storage = [self.RES.Storage.hot_T, self.RES.Storage.low_T, self.RES.Storage.hot_V, self.RES.Storage.low_V]
                self.save_init_.append(storage)
                
            self.initialised = True
            self.RES.alreadyrun = True
            
        else:
            self.RES.Boiler_Tinstruct = None
            self.RES.supplyT = self.save_init_[0]
            self.RES.m_dot = self.save_init_[1]
            self.RES.returnT = self.save_init_[2]
            self.RES.Ts_nodes = self.save_init_[3].copy()
            self.RES.Tr_nodes = self.save_init_[4].copy()
            save = self.save_init_[5]
            if self.storage:
                self.RES.storage_flow = 0
                self.RES.Storage.m_dot = 0
                self.RES.Storage.hot_T, self.RES.Storage.low_T, self.RES.Storage.hot_V, self.RES.Storage.low_V = self.save_init_[6]
                
            for i, (hex, pipe1) in enumerate(self.RES.substations):
                hex.Tr1, hex.Ts2_vrai, hex.Ts1, hex.m_dot1 = save[i][0]['Tr1'], save[i][0]['Ts2_vrai'], save[i][0]['Ts1'], save[i][0]['m1']
                pipe1.pipeS_T = list(save[i][1]).copy()
                pipe1.pipeR_T = list(save[i][2]).copy()
                
    def refined_Ts1(self):
        '''Calculates the optimal linear law (function of external temperature) for Boiler T_supply instruction'''
        upper_limit = 90
        lower_limit = 60
        optimf_Ts1 = lambda x: 70
        opt_lT, opt_uT = 70, 70
        min_score = self.objective_function_Ta(optimf_Ts1)
        total = (upper_limit - lower_limit + 1)
        k_iter = 0
        for uT in np.linspace(lower_limit, upper_limit, 30):
            for lT in np.linspace(lower_limit, uT, (uT -lower_limit)+1):
                f = lambda Ta: uT - ((uT - lT)/27)* (Ta + 7)
                try:
                    f1 = self.objective_function_Ta(f)
                    f2 = self.objective_function_Ta(f)
                    f3 = self.objective_function_Ta(f)
                    f_score = max(f1, f2, f3)
                    #f_score= self.objective_function_Ta(f)

                except ValueError:
                    #print(f'Error, upperT = {uT}, lowerT = {lT} ')
                    f_score = min_score+1
                if f_score < min_score:
                    #print(f_score, min_score)
                    opt_lT, opt_uT = lT, uT
                    min_score = f_score
                    
            k_iter+=1
            percent = (k_iter*100)/total
            i = int((percent * 20) // 100)
            percent = ((percent * 100) //10) /10
            string = '|'+(u"\u2588"*2)*(i+1) + '--'*(19-i) +'| '
            sys.stdout.write(string + f'Refining f_Ts1 -- {percent}%\r')
            sys.stdout.flush()
            
        
        self.f_Ts1 = lambda Ta: opt_uT - ((opt_uT - opt_lT)/27)* (Ta + 7)
        score = self.objective_function_Ta()
        
        print(f'optim score = {score} \n')
        
    def simulation(self, T_instruct, Storage_instruct = None):
        self.initialisation()
        
        self.E_Geo = 0
        self.E_boiler = 0
        self.E_default = 0
        self.cost_Tdefault_SS = [0] * self.nb_SS
        self.cost_constraintT = 0
        self.P_boiler = [] 
        time1 = time.time()
        t_tot = 0
        t_unsupplied = [0] * self.nb_SS
        
        if T_instruct is not None:
            if len(T_instruct) != self.nb_hour:
                raise IndexError('There shall be one instruction for each hour')
        if Storage_instruct is not None:
            nbiter_instruc = len(self.t)//len(Storage_instruct)
            i_instruct = 0
                 
        for j in range(self.nb_hour):
            for p in range(self.nb_SS):
                self.RES.substations[p][0].Ts2 = self.RES.substations[p][0].f_Ts2(self.Ta[j])
            if T_instruct is not None:
                self.RES.Boiler_Tinstruct = T_instruct[j]

            for m in range((3600//dt)):
                if Storage_instruct is not None:
                    if t_tot % nbiter_instruc == 0:
                        self.RES.storage_flow = Storage_instruct[i_instruct]
                        i_instruct += 1
                
                self.RES.iteration()
                
                self.P_boiler.append(self.RES.P_Boiler)
                self.E_Geo += self.RES.P_Geo
                self.E_boiler += self.RES.P_Boiler
                self.E_default += self.RES.P_supplied - self.RES.P_demand
                maxT = self.RES.maxT
                
                if maxT > 95:
                    self.cost_constraintT += (maxT-95)/10
                
                #Calculation of the cost of supply default at time t
                for i, Tdefault in enumerate(self.RES.Tsupply_default_SS):
                    if Tdefault <= 10**(-5):
                        t_unsupplied[i] = t_tot
                    else:
                        ct = ((t_tot - t_unsupplied[i])*dt/60)/15
                        cT = Tdefault / 5 
                        
                        self.cost_Tdefault_SS[i] += ct*cT
                
                t_tot+= 1
                
        time2 = time.time()
    
    def objective_function_optim(self, Instructions):
        T_instruct = Instructions[0:self.nb_hour]
        Stor_instruct = None
        if len(Instructions) > self.nb_hour:
            Stor_instruct = Instructions[24::]
        return self.objective_function(T_instruct, Storage_instructions = Stor_instruct)
            
    def objective_function(self, T_boiler_instruction, Storage_instructions = None, exe_time = False):
        time1 = time.time()
        
        if Storage_instructions is not None:
            if not self.storage:
                raise ValueError('There is no storage buffer in this DH Network')
            elif len(Storage_instructions) < self.nb_hour:
                raise ValueError('There shall be at least on instruction for each hour')

        #self.simulation(T_boiler_instruction)
        try:
            self.simulation(T_boiler_instruction, Storage_instruct = Storage_instructions)
        except ValueError as e:
            return max(10, float(e.__str__())*20)
            
        #kilowatt-hour
        E_boiler = self.E_boiler*dt/3600
        E_Geo = self.E_Geo*dt/3600
        E_tot = E_boiler + E_Geo
        E_default = self.E_default*dt/3600
        
        C1 = E_boiler/(E_tot) #* C_kWh * E_boiler
        
        C2 = sum(self.cost_Tdefault_SS)/(len(self.t)*self.nb_SS)
        
        C_constraint = self.cost_constraintT/(len(self.t)*self.nb_SS)
        
        print(f'coût conso = {C1}, coût ecart consigne = {C2}, coût T_maximale = {C_constraint}')
        print(f'E_boiler = {E_boiler}, E_Geothermie = {E_Geo}')
        print(f'Coût total = {C1 + C2 + C_constraint}')
        
        time2 = time.time()
        
        if exe_time:
            print(time2-time1)
        #return C1 + C2 + C_constraint
        
    
    def plot(self, T_instruct, Storage_instruct = None):
        self.initialisation()
        
        T_return_SS = []
        T_supplySS = []
        T_supply_secondary = []
        T_secondary_demand = []
        P_boiler = []
        water_flow = []
        storage_flow = []
        E_boiler = 0
        
        storage_hV = []
        storage_lV = []
        storage_hT = []
        
        time1 = time.time()
        t_tot = 0
        
        if T_instruct is not None:
            if len(T_instruct) != self.nb_hour:
                raise IndexError('There shall be one instruction for each hour')
        if Storage_instruct is not None:
            if not self.storage:
                raise ValueError('There is no storage buffer in this DH Network')
            elif len(Storage_instruct) < self.nb_hour:
                raise ValueError('There shall be at least on instruction for each hour')
            nbiter_instruc = len(self.t)//len(Storage_instruct)
            i_instruct = 0
                
        for j in range(self.nb_hour):
            for p in range(self.nb_SS):
                self.RES.substations[p][0].Ts2 = self.RES.substations[p][0].f_Ts2(self.Ta[j])
            if T_instruct is not None:
                self.RES.Boiler_Tinstruct = T_instruct[j]
            
            for m in range(3600//dt):
                if Storage_instruct is not None:
                    if t_tot % nbiter_instruc == 0:
                        self.RES.storage_flow = Storage_instruct[i_instruct]
                        i_instruct += 1
                
                self.RES.iteration()
                
                P_boiler.append(self.RES.P_Boiler)
                E_boiler += self.RES.P_Boiler
                T_return_SS.append([X[0].Tr1 for X in self.RES.substations])
                T_supplySS.append([X[0].Ts1 for X in self.RES.substations])
                T_secondary_demand.append([X[0].Ts2 for X in self.RES.substations])
                T_supply_secondary.append([X[0].Ts2_vrai for X in self.RES.substations])
                water_flow.append(self.RES.m_dot)
                storage_flow.append(self.RES.storage_flow)
                
                storage_hV.append(self.RES.Storage.hot_V)
                storage_lV.append(self.RES.Storage.low_V)
                
                storage_hT.append(self.RES.Storage.hot_T)
                
                t_tot += 1
              
        time2 = time.time()
        
        print(E_boiler * dt)
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
        plt.plot(t, [a for a in storage_flow])
        plt.show()
        
        plt.figure()
        plt.title(f'Storage Volume')
        plt.plot(t, [a for a in storage_hV], label = 'hV')
        plt.plot(t, [a for a in storage_lV], label = 'lV')
        plt.legend()
        plt.show()
        
        plt.figure()
        plt.title(f'Storage hT')
        plt.plot(t, [a for a in storage_hT], label = 'hT')
        plt.legend()
        plt.show()
        
        
    def objective_function_Ta(self, f_Ts1 = None, print_time = False):
        '''evaluates a given linear function (or if not furnished self.f_Ts1) as instruction for T_supply instead of a list of instructions calling self.objective_function()'''
        
        if f_Ts1 == None:
            T_instruct = [self.f_Ts1(T) for T in self.Ta]
        else:
            T_instruct = [f_Ts1(T) for T in self.Ta]
            
        return self.objective_function(np.array(T_instruct), exe_time = print_time)
        
        
    def plot_Ta(self, f_Ts1 = None):
        if f_Ts1 == None:
            T_instruct = [self.f_Ts1(T) for T in self.Ta]
        else:
            T_instruct = [f_Ts1(T) for T in self.Ta]
        
        self.plot(np.array(T_instruct))
        



