import numpy as np
from numpy import log as ln
from itertools import *
from Components.Ressources import *

class Network:
    def __init__(self, source, liste_substations):
        ''' inlet_T initialization supply Temperature of the network
        liste_substations : [[HEX_i, Pipe_nodetopreviousnode_i, Pipe_hextonode_i]]'''
        self.supplyT = source.Ts_Net
        self.returnT = 40
        self.subm_dot = [X[0].m_dot1 for X in liste_substations]
        self.m_dot = sum([X[0].m_dot1 for X in liste_substations])
        self.Ts_nodes = [X[1].TS_ext() for X in liste_substations]
        self.Tr_nodes = [X[1].TR_ext() for X in liste_substations]
        self.nb_substations =  len(liste_substations)
        self.substations = liste_substations
        self.src = source
        
        self.alreadyrun = False
        self.NETtype = 'basic'
        
    def iter_returnside(self):
        Tr_nodes = self.Tr_nodes
        Tr_nodes_new = [0] * len(self.Tr_nodes)
        Tr_node_network_upstream = 0
        m_dotR = 0 
        
        for i in range(len(self.substations)):
            #Calculation of the return nodes temperature in the network at time (t+1)
            m = -(i+1) #we must iter from further node to closer
            h, p1, p2 = self.substations[m]
            m_dotSSr = h.m_dot1
            Tr = h.Tr1
            p2.evolR_T(m_dotSSr,Tr) 
            Tr_SSconfluence = p2.TR_ext() #return temperature from SS at confluence with network for time t based on Tr at time (t-1) and m_dot at time (t)
            
            Tr_nodes_new[m] = (m_dotSSr*Tr_SSconfluence + m_dotR * Tr_node_network_upstream)/(m_dotSSr + m_dotR)
            m_dotR += m_dotSSr
            p1.evolR_T(m_dotR, Tr_nodes[m])
            Tr_node_network_upstream = p1.TR_ext()
            
        self.returnT = Tr_node_network_upstream
        self.Tr_nodes = np.copy(Tr_nodes_new)
        
    def iter_supplyside(self):
        m_dot = self.m_dot
        T_node = [self.supplyT] + list(self.Ts_nodes)
        
        for i, (hex, pipe1, pipe2) in enumerate(self.substations):
            #Calculation of Temperatures in the network at time (t+1) (for next iteration)
            pipe1.evolS_T(m_dot, T_node[i])
            
            #Calculation of the new mass flow and return temperature at substation for time t+1
            m_dotSS = hex.m_dot1
            pipe2.evolS_T(m_dotSS, T_node[i+1])
            Ts1 = pipe2.TS_ext()
            hex.solve(Ts1)
            m_dot -= m_dotSS
            
            
    def iteration(self):
        '''we consider that mass flow is established much faster than temperature flow''' #?
        
        #RETURN SIDE
        self.iter_returnside()
        #SUPPLY SIDE
        self.iter_supplyside()
            
        self.subm_dot = [X[0].m_dot1 for X in self.substations] #substations mass flows at time t+1
        self.m_dot = sum(self.subm_dot) #Network mass flow at time t+1
        self.Ts_nodes = [X[1].TS_ext() for X in self.substations]
        
        #Calculation of the supply temperature re-heated by the source
        self.src.solve(self.m_dot, self.returnT)
        self.supplyT = self.src.Ts_Net
        
        
        
class Network_iter(Network):
    def __init__(self, source, liste_substations):
        Network.__init__(self, source, liste_substations)
        self.NETtype = 'iter'
    
    def iter_supplyside(self):
        m_dot = self.m_dot
        T_node = [self.supplyT] + list(self.Ts_nodes)
        
        save = [] #Initial conditions
        for i, (hex, pipe1, pipe2) in enumerate(self.substations):
            save.append(({'Tr1' : hex.Tr1, 'Ts2_vrai' : hex.Ts2_vrai, 'Ts1' : hex.Ts1}, np.copy(pipe1.pipeS_T), np.copy(pipe2.pipeS_T)))
        mdot0 = self.m_dot
        m_dot = mdot0
        
        #first iteration outside of the 'while' loop
        for i, (hex, pipe1, pipe2) in enumerate(self.substations):
            #Calculation of Temperatures in the network at time (t+1) (for next iteration)
            pipe1.evolS_T(m_dot, T_node[i])
            
            #Calculation of the new mass flow and return temperature at substation for time t+1
            m_dotSS = hex.m_dot1
            pipe2.evolS_T(m_dotSS, T_node[i+1])
            Ts1 = pipe2.TS_ext()
            hex.solve(Ts1)
            m_dot -= m_dotSS
            
        m_dot1 = sum(X[0].m_dot1 for X in self.substations)
        n = 0
        while np.abs(m_dot1 - mdot0) > 0.01 and n<1000:
            mdot0 = m_dot1
            #print(mdot0)
            m_dot = mdot0
            for i, (hex, pipe1, pipe2) in enumerate(self.substations):
                hex.Tr1, hex.Ts2_vrai, hex.Ts1 = save[i][0]['Tr1'], save[i][0]['Ts2_vrai'], save[i][0]['Ts1']
                pipe1.pipeS_T = list(save[i][1])
                pipe2.pipeS_T = list(save[i][2])
            
            for i, (hex, pipe1, pipe2) in enumerate(self.substations):
                #Calculation of Temperatures in the network at time (t+1) (for next iteration)
                pipe1.evolS_T(m_dot, T_node[i])
                
                #Calculation of the new mass flow and return temperature at substation for time t+1
                m_dotSS = hex.m_dot1
                pipe2.evolS_T(m_dotSS, T_node[i+1])
                Ts1 = pipe2.TS_ext()
                hex.solve(Ts1)
                m_dot -= m_dotSS
            m_dot1 = sum(X[0].m_dot1 for X in self.substations)
            n+=1
        print(n)
        
    
    
class Network_bOptim(Network):
    def __init__(self, source, P_boiler, liste_substations):
        Network.__init__(self, source, liste_substations)
        self.NETtype = 'Optim_boiler'
        self.P_boiler = P_boiler
        self.P_Geo = 0
        self.P_demand = 0
        self.P_supplied = 0
        self.Tsupply_default_SS = []
        self.maxT = 0
        
    def iter_supplyside(self):
        self.maxT = 0
        m_dot = self.m_dot
        
        supplyT_reheated = self.supplyT + self.P_boiler/(Cp * m_dot)
        T_node = [supplyT_reheated] + list(self.Ts_nodes)
        self.maxT = supplyT_reheated
            
        self.P_Geo = self.src.P
        self.P_demand = 0
        self.P_supplied = 0
        self.Tsupply_default_SS = []
        
        for i, (hex, pipe1, pipe2) in enumerate(self.substations):
            #Calculation of Temperatures in the network at time (t+1) (for next iteration)
            pipe1.evolS_T(m_dot, T_node[i])
            
            #Calculation of the new mass flow and return temperature at substation for time t+1
            m_dotSS = hex.m_dot1
            pipe2.evolS_T(m_dotSS, T_node[i+1])
            Ts1 = pipe2.TS_ext()
            hex.solve(Ts1)
            m_dot -= m_dotSS
            
            Pd = hex.m_dot2 * Cp * (hex.Ts2 - hex.Tr2)
            Ps = hex.m_dot2 * Cp * (hex.Ts2_vrai - hex.Tr2)
            self.Tsupply_default_SS.append(hex.Ts2 - hex.Ts2_vrai)
            self.P_demand += Pd
            self.P_supplied += Ps
            
        
class Network_boiler(Network):
    def __init__(self, source, P_boiler, liste_substations):
        Network.__init__(self, source, liste_substations)
        self.NETtype = 'boiler'
        self.P_boiler = P_boiler
        self.boiler = False
        
    def iter_supplyside(self):
        m_dot = self.m_dot
        T_node = [self.supplyT] + list(self.Ts_nodes)
        
        save = [] #Initial conditions
        for i, (hex, pipe1, pipe2) in enumerate(self.substations):
            save.append(({'Tr1' : hex.Tr1, 'Ts2_vrai' : hex.Ts2_vrai, 'Ts1' : hex.Ts1, 'm1':hex.m_dot1}, np.copy(pipe1.pipeS_T), np.copy(pipe2.pipeS_T)))
        mdot0 = self.m_dot
        m_dot = mdot0
        
        BoilerOn = False
        for i, (hex, pipe1, pipe2) in enumerate(self.substations):
            #Calculation of Temperatures in the network at time (t+1) (for next iteration)
            pipe1.evolS_T(m_dot, T_node[i])
            
            #Calculation of the new mass flow and return temperature at substation for time t+1
            m_dotSS = hex.m_dot1
            pipe2.evolS_T(m_dotSS, T_node[i+1])
            Ts1 = pipe2.TS_ext()
            hex.solve(Ts1)
            m_dot -= m_dotSS
            
            #Boiler Ignition
            if (hex.Ts2 - hex.Ts2_vrai) > 1:
                BoilerOn = True
        
        if BoilerOn:
            m_dot = mdot0
            supplyT_reheated = self.supplyT + self.P_boiler/(Cp * m_dot)
            T_node[0] = supplyT_reheated
            
            for i, (hex, pipe1, pipe2) in enumerate(self.substations):
                hex.Tr1, hex.Ts2_vrai, hex.Ts1, hex.m_dot1 = save[i][0]['Tr1'], save[i][0]['Ts2_vrai'], save[i][0]['Ts1'], save[i][0]['m1']
                pipe1.pipeS_T = list(save[i][1])
                pipe2.pipeS_T = list(save[i][2])
                
            for i, (hex, pipe1, pipe2) in enumerate(self.substations):
                #Calculation of Temperatures in the network at time (t+1) (for next iteration)
                pipe1.evolS_T(m_dot, T_node[i])
                
                #Calculation of the new mass flow and return temperature at substation for time t+1
                m_dotSS = hex.m_dot1
                pipe2.evolS_T(m_dotSS, T_node[i+1])
                Ts1 = pipe2.TS_ext()
                hex.solve(Ts1)
                m_dot -= m_dotSS
            
        self.boiler = BoilerOn
        