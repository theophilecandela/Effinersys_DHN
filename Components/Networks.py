import numpy as np
from numpy import log as ln
from itertools import *
from Components.Ressources import *

class Network():
    def __init__(self, source, list_substations, storage_buffer = None):
        ''' inlet_T initialization supply Temperature of the network
        list_substations : [[HEX_i, Pipe_nodetopreviousnode_i]]'''
        self.supplyT = source.Ts_Net
        self.returnT = list_substations[0][1].TR_ext()
        self.subm_dot = [X[0].m_dot1 for X in list_substations]
        self.m_dot = sum([X[0].m_dot1 for X in list_substations])
        self.Ts_nodes = [X[1].TS_ext() for X in list_substations]
        self.Tr_nodes = [0] + [X[1].TR_ext() for X in list_substations[::-1]] #from further node to closer
        self.nb_substations =  len(list_substations)
        
        self.substations = list_substations
        self.src = source
        self.Storage = storage_buffer
        self.alreadyrun = False
        
        self.storage_flow = 0
        self.Boiler_Tinstruct = None
        
        self.P_Boiler = 0
        self.P_Geo = 0
        self.P_demand = 0
        self.P_supplied = 0
        self.Tsupply_default_SS = []
        self.maxT = 0
        
        self.NETtype = 'Optim_boiler_storage'
        
    def iter_returnside(self):
        Tr_node_network_upstream = 0
        m_dotR = 0 
        
        for i, (h, p1) in enumerate(self.substations[::-1]):
            #Calculation of the return nodes temperature in the network at time (t+1)
            m_dotSSr = h.m_dot1
            Tr_SS = h.Tr1
            
            T_mix_node = (m_dotR*self.Tr_nodes[i] + m_dotSSr*Tr_SS)/(m_dotR + m_dotSSr)
            p1.evolR_T(m_dotR + m_dotSSr, T_mix_node)
            self.Tr_nodes[i] = Tr_node_network_upstream
            m_dotR += m_dotSSr
            
            Tr_node_network_upstream = p1.TR_ext()
        
        self.Tr_nodes[-1] = Tr_node_network_upstream #= self.returnT
        self.returnT = Tr_node_network_upstream    
        
        
    def iter_supplyside(self):
        self.maxT = 0
        m_dot = self.m_dot
        self.P_Boiler = 0
        self.P_Geo = self.src.P
        self.P_demand = 0
        self.P_supplied = 0
        self.Tsupply_default_SS = []
        
        if self.Boiler_Tinstruct != None:
            if self.Boiler_Tinstruct > self.supplyT:
                self.P_Boiler = (self.Boiler_Tinstruct - self.supplyT)*Cp*m_dot 
            else: 
                self.P_Boiler = 0
                        
        supplyT_reheated = self.supplyT + self.P_Boiler/(Cp * m_dot) #=self.Boiler_Tinstruct
        T_node = [supplyT_reheated] + list(self.Ts_nodes)
        self.maxT = supplyT_reheated
            

        
        for i, (hex, pipe1) in enumerate(self.substations):
            #Calculation of Temperatures in the network at time (t+1) (for next iteration)
            pipe1.evolS_T(m_dot, T_node[i])
            
            #Calculation of the new mass flow and return temperature at substation for time t+1
            m_dotSS = hex.m_dot1
            
            hex.solve(pipe1.TS_ext())
            self.Ts_nodes[i] = pipe1.TS_ext()
            
            m_dot -= m_dotSS
            
            Pd = hex.m_dot2 * Cp * (hex.Ts2 - hex.Tr2)
            Ps = hex.m_dot2 * Cp * (hex.Ts2_vrai - hex.Tr2)
            self.Tsupply_default_SS.append(hex.Ts2 - hex.Ts2_vrai)
            self.P_demand += Pd
            self.P_supplied += Ps
    
    
    def storage_cold(self):
        mdot = 0
        if self.storage_flow > 0:
            if np.abs(self.storage_flow) >= self.m_dot:
                raise ValueError((np.abs(self.storage_flow)/self.m_dot)*20)
                
            mdot = self.Storage.intake_cold_water(np.abs(self.storage_flow), self.returnT)
            self.m_dot -= mdot
            
        elif self.storage_flow < 0:
            T, mdot = self.Storage.delivery_cold_water(np.abs(self.storage_flow))
            self.returnT = (self.returnT * self.m_dot + T * mdot)/(self.m_dot + mdot)
            self.m_dot += mdot
            mdot = -mdot
            
        self.Storage.m_dot = mdot
            
            
    def storage_hot(self):
        if self.Storage.m_dot > 0:
            T , mdot = self.Storage.delivery_hot_water()
            self.supplyT = (self.supplyT * self.m_dot + T * mdot)/(self.m_dot + mdot)
            self.m_dot += mdot
            
        elif self.Storage.m_dot < 0:
            self.Storage.intake_hot_water(self.supplyT)
            self.m_dot += self.Storage.m_dot
            
            
    def iteration(self):
        '''we consider that mass flow is established much faster than temperature flow
           heat exchanges in HEX are considered instantaneous''' 
        
        #STORAGE AND GEOTHERMAL HEX
        if self.Storage is not None:
            self.storage_cold()
        
        self.src.solve(self.m_dot, self.returnT)
        self.supplyT = self.src.Ts_Net
        
        # print(self.src.P)
        # print(Cp*self.m_dot*(self.supplyT - self.returnT))
        
        if self.Storage is not None:
            self.storage_hot()
        
        #RETURN SIDE
        self.iter_returnside()
        
        #SUPPLY SIDE
        self.iter_supplyside()
            
        self.subm_dot = [X[0].m_dot1 for X in self.substations] #substations mass flows at time t+1
        self.m_dot = sum(self.subm_dot) #Network mass flow at time t+1
         
        #Calculation of the supply temperature re-heated by the source
        
        
        
            