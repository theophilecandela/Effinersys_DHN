import numpy as np
from numpy import log as ln
from itertools import *
from Components.Ressources import *

class Network:
    def __init__(self, source, list_substations):
        ''' inlet_T initialization supply Temperature of the network
        list_substations : [[HEX_i, Pipe_nodetopreviousnode_i]]'''
        self.supplyT = source.Ts_Net
        self.returnT = list_substations[0][1].TR_ext()
        self.subm_dot = [X[0].m_dot1 for X in list_substations]
        self.m_dot = sum([X[0].m_dot1 for X in list_substations])
        self.Ts_nodes = [X[1].TS_ext() for X in list_substations]
        self.Tr_nodes = [X[1].TR_ext() for X in list_substations[::-1]] #from further node to closer
        self.nb_substations =  len(list_substations)
        self.substations = list_substations
        self.src = source
        
        self.alreadyrun = False
        self.NETtype = 'basic'
        
    def iter_returnside(self):
        Tr_node_network_upstream = 0
        m_dotR = 0 
        
        for i, (h, p1) in enumerate(self.substations[::-1]):
            #Calculation of the return nodes temperature in the network at time (t+1)
            m_dotSSr = h.m_dot1
            Tr_SS = h.Tr1
            
            p1.evolR_T(m_dotR + m_dotSSr, self.Tr_nodes[i])
            self.Tr_nodes[i] = ((m_dotSSr*Tr_SS + m_dotR * Tr_node_network_upstream)/(m_dotSSr + m_dotR))
            m_dotR += m_dotSSr
            
            Tr_node_network_upstream = p1.TR_ext()
            
        self.returnT = Tr_node_network_upstream
       
        
    def iter_supplyside(self):
        m_dot = self.m_dot
        T_node = [self.supplyT] + list(self.Ts_nodes)
        
        for i, (hex, pipe1) in enumerate(self.substations):
            #Calculation of Temperatures in the network at time (t+1) (for next iteration)
            pipe1.evolS_T(m_dot, T_node[i])
            
            #Calculation of the new mass flow and return temperature at substation for time t+1
            m_dotSS = hex.m_dot1
            hex.solve(T_node[i+1])
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
        
        
        
class Network_bOptim(Network):
    def __init__(self, source, P_boiler, list_substations):
        Network.__init__(self, source, list_substations)
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
        
        for i, (hex, pipe1) in enumerate(self.substations):
            #Calculation of Temperatures in the network at time (t+1) (for next iteration)
            pipe1.evolS_T(m_dot, T_node[i])
            
            #Calculation of the new mass flow and return temperature at substation for time t+1
            m_dotSS = hex.m_dot1
            
            hex.solve(T_node[i+1])
            m_dot -= m_dotSS
            
            Pd = hex.m_dot2 * Cp * (hex.Ts2 - hex.Tr2)
            Ps = hex.m_dot2 * Cp * (hex.Ts2_vrai - hex.Tr2)
            self.Tsupply_default_SS.append(hex.Ts2 - hex.Ts2_vrai)
            self.P_demand += Pd
            self.P_supplied += Ps
            