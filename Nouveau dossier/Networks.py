import numpy as np
from numpy import log as ln
from itertools import *

class Network:
    def __init__(self, inlet_T, source, liste_substations):
        ''' inlet_T initialization supply Temperature of the network
        liste_substations : [HEX, Pipe_nodetopreviousnode, Pipe_hextonode]'''
        self.supplyT = inlet_T
        self.returnT = 40
        self.subm_dot = [X[0].m_dot1 for X in liste_substations]
        self.m_dot = sum([X[0].m_dot1 for X in liste_substations])
        self.Ts_nodes = [X[1].TS_ext() for X in liste_substations]
        self.Tr_nodes = [X[1].TR_ext() for X in liste_substations]
        self.nb_substations =  len(liste_substations)
        self.substations = liste_substations
        self.src = source
        
        
    def iteration(self):
        #we consider that mass flow is established much faster than temperature flow
        m_dot = self.m_dot
        T_node = [self.supplyT] + list(self.Ts_nodes)
        Tr_nodes = self.Tr_nodes
        Ts_nodes_new = []
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
            
        
        for i, (hex, pipe1, pipe2) in enumerate(self.substations):
            #Calculation of Temperatures in the network at time (t+1) (for next iteration)
            pipe1.evolS_T(m_dot, T_node[i])
            Ts_nodes_new.append(pipe1.TS_ext()) 

            #Calculation of the new mass flow and return temperature at substation for time t+1
            m_dotSS = hex.m_dot1
            pipe2.evolS_T(m_dotSS, T_node[i+1])
            Ts1 = pipe2.TS_ext()
            hex.solve(Ts1)
            m_dot -= m_dotSS
            
        self.returnT = Tr_node_network_upstream
        self.subm_dot = [X[0].m_dot1 for X in self.substations] #substations mass flows at time t+1
        self.m_dot = sum(self.subm_dot) #Network mass flow at time t+1
        self.Ts_nodes = np.copy(Ts_nodes_new)
        self.Tr_nodes = np.copy(Tr_nodes_new)
        
        #Calculation of the supply temperature re-heated by the source
        self.src.solve(self.m_dot, self.returnT)
        self.supplyT = self.src.Ts_Net
        
        
class Network_iter(Network):
    def iteration(self):
        pass
    
class Network_boiler(Network):
    def iteration(self):
        pass