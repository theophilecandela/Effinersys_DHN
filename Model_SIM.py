import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
from random import gauss
from itertools import *
import time

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

    def initialisation(self, P_instruct = [0]):
        # First loop for initialisation
        if not self.initialised:
            t_tot = 0
            nbiter_instruc = len(self.t)//len(P_instruct)
            i_instruct = 0
            
            for j in range(self.nb_hour):
                for p, T_h in enumerate(self.Ts2_h):
                    self.RES.substations[p][0].Ts2 = T_h[j]
                for m in range(3600//dt):
                    if t_tot % nbiter_instruc == 0:
                        self.RES.P_boiler = P_instruct[i_instruct]
                        i_instruct += 1
                        
                    self.RES.iteration()
                    t_tot += 1
                    
            self.save_init_ = [self.RES.supplyT, self.RES.m_dot, self.RES.returnT, self.RES.Ts_nodes, self.RES.Tr_nodes] #Initial conditions
            save =  []
            for i, (hex, pipe1) in enumerate(self.RES.substations):
                save.append(({'Tr1' : hex.Tr1, 'Ts2_vrai' : hex.Ts2_vrai, 'Ts1' : hex.Ts1, 'm1':hex.m_dot1}, np.copy(pipe1.pipeS_T), np.copy(pipe1.pipeR_T) ))
            self.save_init_.append(save)
            self.initialised = True
            self.RES.alreadyrun = True
            
        else:
            self.RES.supplyT = self.save_init_[0]
            self.RES.m_dot = self.save_init_[1]
            self.RES.returnT = self.save_init_[2]
            self.RES.Ts_nodes = self.save_init_[3]
            self.RES.Tr_nodes = self.save_init_[4]
            save = self.save_init_[5]
            
            for i, (hex, pipe1) in enumerate(self.RES.substations):
                hex.Tr1, hex.Ts2_vrai, hex.Ts1, hex.m_dot1 = save[i][0]['Tr1'], save[i][0]['Ts2_vrai'], save[i][0]['Ts1'], save[i][0]['m1']
                pipe1.pipeS_T = list(save[i][1])
                pipe1.pipeR_T = list(save[i][2])
                
                
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
                    #self.cost_constraintT += 3*(np.exp(0.2*(maxT-105)) - np.exp(-2)) 
                    self.cost_constraintT += (maxT-95)/10
                
                #Calculation of the cost of supply default at time t
                for i, Tdefault in enumerate(self.RES.Tsupply_default_SS):
                    if Tdefault <= 10**(-5):
                        t_unsupplied[i] = t_tot
                    else:
                        #print(i, len(t_unsupplied))
                        ct = ((t_tot - t_unsupplied[i])*dt/60)/15
                        cT = Tdefault / 5 
                        
                        self.cost_Tdefault_SS[i] += ct*cT
                
                t_tot+= 1
                
        time2 = time.time()
    
    def objective_function(self, P_boiler_instruction, exe_time = False):
        time1 = time.time()

        #self.simulation(P_boiler_instruction)
        try:
            self.simulation(P_boiler_instruction)
        except ValueError as e:
            return float(e.__str__())*20
            
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
        

## Ts1 linear function of Ta
# The next class, herited from original class Simulation, have new methods and attributes, unlocking the ability to calculate and a linear law between outside Temperature and supplied water temperature instruction. 

class Simulation_Ta(Simulation):

    def __init__(self, RES, Ta):
        '''Ta:  list of external temperature per hour'''
        Ts2_h = []
        for p, SS in enumerate(RES.substations):
            Ts2_h.append(SS[0].f_Ts2(Ta))
        Simulation.__init__(self, RES, Ts2_h)
        self._Ta = Ta
        self.f_Ts1 = lambda x: 70
        self.PBoiler_Ta = [] #This attribute, containing the power instruction resulting from the water law allows us to directly compare the instructions given by the linear law and by the optimisation process, as well as reducing the dimension of the problem for optimisation.
        
    def _get_Ta(self):
        return self._Ta
        
    def _set_Ta(self, Ta):
        '''when user wants to modify Ta, this method is called and consequently modifies the duration of simulation'''
        Ts2_h = []
        for p, SS in enumerate(self.RES.substations):
            Ts2_h.append(SS[0].f_Ts2(Ta))
        self.Ts2_h = Ts2_h

    Ta = property(_get_Ta, _set_Ta)
            
    def refined_Ts1(self, dim = (1440*60)//dt):
        upper_limit = 90
        lower_limit = 60
        optimf_Ts1 = lambda x: 70
        opt_lT, opt_uT = 70, 70
        min_score = self.objective_function_Ta(optimf_Ts1)
        for uT in np.linspace(lower_limit, upper_limit, 30):
            for lT in np.linspace(lower_limit, uT, (uT -lower_limit)+1):
                f = lambda Ta: uT - ((uT - lT)/27)* (Ta + 7)
                try:
                    f_score= self.objective_function_Ta(f, dim_perday = dim)

                except ValueError:
                    #print(f'Error, upperT = {uT}, lowerT = {lT} ')
                    f_score = min_score+1
                if f_score < min_score:
                    #print(f_score, min_score)
                    opt_lT, opt_uT = lT, uT
                    min_score = f_score
        
        optimf_Ts1 = lambda Ta: opt_uT - ((opt_uT - opt_lT)/27)* (Ta + 7)
        self.f_Ts1 = optimf_Ts1
        score = self.objective_function_Ta(dim_perday = dim)
        
        print(f'optim score = {score} \n')
        
        
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
                
                if maxT > 95:
                    #self.cost_constraintT += 3*(np.exp(0.2*(maxT-105)) - np.exp(-2))
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
        
        
    def objective_function_Ta(self, f_Ts1 = None, dim_perday = (1440*60)//dt):
        #Reinitialisation of self.PBoiler_Ta
        self.PBoiler_Ta = []
        time1 = time.time()
        
        nb_instruct = dim_perday * self.nb_hour//24
        if nb_instruct  > len(self.t):
            raise ValueError("Instructions can't be more precise than time step of pipe model")
        else:
            nbiter_per_instruc = len(self.t)//nb_instruct
        
        try:
            self.simulation_Ta(f_Ts1, nbiter_per_instruc)
        except ValueError as e:
            return float(e.__str__())*10
        
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
    
    
    def plot_Ta(self, dim_perday = (1440*60)//dt):
        self.initialisation()
        
        nb_instruct = dim_perday * self.nb_hour//24
        if nb_instruct  > len(self.t):
            raise ValueError("Instructions can't be more precise than time step of pipe model")
        else:
            nbiter_per_instruc = len(self.t)//nb_instruct
            
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
                
                if t_tot % nbiter_per_instruc == 0:
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




