import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
from random import gauss
from itertools import *
import time

from geneticalgorithm import geneticalgorithm as ga
from Components.Data_meteo_reading import Ts2_A, Ts2_B, Ts2_C, Ts2_D
from Components.Data_meteo_reading import Ta as ext_T
from Components.Networks import Network_bOptim
from Components.Source import Source
from Components.Heat_exchanger import HEX_nom
from Components.Pipes import Pipe
from Components.Ressources import *
from Model_SIM import Simulation_Ta



##Test iter_returnside
def f1():
    v1_11 = []
    v1_12 = []
    v1_3 = []
    v1_31 = []
    for jas in range(100):
    
        SS_A = [HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
        SS_B = [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 75)]
        SS_C = [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
        
        liste_substations = [SS_A, SS_B, SS_C]
        
        
        Tr_nodes = [X[1].TR_ext() for X in liste_substations]
        t1 = time.time()
        for k in range(1000):
            Tr_nodes_new = [0] * len(Tr_nodes)
            Tr_node_network_upstream = 0
            m_dotR = 0 
            TR_SS= []
            for i in range(len(liste_substations)):
                #Calculation of the return nodes temperature in the network at time (t+1)
                m = -(i+1) #we must iter from further node to closer
                h, p1, p2 = liste_substations[m]
                m_dotSSr = h.m_dot1
                Tr = h.Tr1
                p2.evolR_T(m_dotSSr,Tr) 
                Tr_SSconfluence = p2.TR_ext() #return temperature from SS at confluence with network for time t based on Tr at time (t-1) and m_dot at time (t)
                TR_SS.append((Tr,Tr_SSconfluence))
                Tr_nodes_new[m] = (m_dotSSr*Tr_SSconfluence + m_dotR * Tr_node_network_upstream)/(m_dotSSr + m_dotR)
                m_dotR += m_dotSSr
                p1.evolR_T(m_dotR, Tr_nodes[m])
                Tr_node_network_upstream = p1.TR_ext()
                
            Tr_nodes = np.copy(Tr_nodes_new)
        t2 = time.time()
        a = (t2-t1)
        # print((t2-t1)/1000)
        # print(Tr_node_network_upstream,'\n', Tr_nodes_new[::-1], '\n',  TR_SS)
        
        
        ## Comparasion entre version initiale et version avec enumerate elegante
        SS_A = [HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
        SS_B = [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 75)]
        SS_C = [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
        
        liste_substations = [SS_A, SS_B, SS_C]
        
        Tr_nodes = [X[1].TR_ext() for X in liste_substations]
        t1 = time.time()
        for k in range(1000):
            Tr_nodes_new = []
            Tr_node_network_upstream = 0
            m_dotR = 0 
            TR_SS= []
        
            for i, (h, p1, p2) in enumerate(liste_substations[::-1]):
                #Calculation of the return nodes temperature in the network at time (t+1)
                
                m_dotSSr = h.m_dot1
                Tr = h.Tr1
                p2.evolR_T(m_dotSSr,Tr) 
                Tr_SSconfluence = p2.TR_ext()
                
                TR_SS.append((Tr,Tr_SSconfluence))
                Tr_nodes_new.append((m_dotSSr*Tr_SSconfluence + m_dotR * Tr_node_network_upstream)/(m_dotSSr + m_dotR))
                m_dotR += m_dotSSr
                p1.evolR_T(m_dotR, Tr_nodes[i])
                Tr_node_network_upstream = p1.TR_ext()
            
            Tr_nodes = np.copy(Tr_nodes_new)
            
        t2 = time.time()
        b = (t2-t1)
        # print('________')
        # print((t2-t1)/1000)
        # print(Tr_node_network_upstream,'\n', Tr_nodes_new, '\n',  TR_SS)
        # 
        # print(b/a)
        v1_11.append(b/a)
        
    ## Comparaison version init avec un ou deux tuyau
        SS_A = [HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000)]
        SS_B = [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350)]
        SS_C = [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500)]
        
        liste_substations = [SS_A, SS_B, SS_C]
        
        Tr_nodes = [X[1].TR_ext() for X in liste_substations]
        t1 = time.time()
        for k in range(1000):
            Tr_nodes_new = [0] * len(Tr_nodes)
            Tr_node_network_upstream = 0
            m_dotR = 0 
            TR_SS= []
            for i in range(len(liste_substations)):
                #Calculation of the return nodes temperature in the network at time (t+1)
                m = -(i+1) #we must iter from further node to closer
                h, p1 = liste_substations[m]
                m_dotSSr = h.m_dot1
                Tr_SS = h.Tr1
                
                TR_SS.append(Tr_SS)
                Tr_nodes_new[m] = (m_dotSSr*Tr_SS + m_dotR * Tr_node_network_upstream)/(m_dotSSr + m_dotR)
                m_dotR += m_dotSSr
                p1.evolR_T(m_dotR, Tr_nodes[m])
                Tr_node_network_upstream = p1.TR_ext()
                
            Tr_nodes = np.copy(Tr_nodes_new)
        t2 = time.time()
        b = (t2-t1)
        # print('________')
        # print((t2-t1)/1000)
        # print(Tr_node_network_upstream,'\n', Tr_nodes_new, '\n',  TR_SS)
        # 
        # print(b/a)
        v1_12.append(b/a)
        
        ## Comparasion entre version initiale et nouvelle version avec enumerate et un seul tuyau
        
        SS_A = [HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000)]
        SS_B = [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350)]
        SS_C = [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500)]
        
        liste_substations = [SS_A, SS_B, SS_C]
        
        Tr_nodes = [X[1].TR_ext() for X in liste_substations]
        t1 = time.time()
        for k in range(1000):
            Tr_nodes_new = []
            Tr_node_network_upstream = 0
            m_dotR = 0 
            TR_SS= []
        
            for i, (h, p1) in enumerate(liste_substations[::-1]):
                #Calculation of the return nodes temperature in the network at time (t+1)
                
                m_dotSSr = h.m_dot1
                Tr_SS = h.Tr1
                
                TR_SS.append(Tr_SS)
                Tr_nodes_new.append((m_dotSSr*Tr_SS + m_dotR * Tr_node_network_upstream)/(m_dotSSr + m_dotR))
                m_dotR += m_dotSSr
                p1.evolR_T(m_dotR, Tr_nodes[i])
                Tr_node_network_upstream = p1.TR_ext()
            
            Tr_nodes = np.copy(Tr_nodes_new)
            
        t2 = time.time()
        c = (t2-t1)
        
        # print('________')
        # print((t2-t1)/1000)
        # print(Tr_node_network_upstream,'\n', Tr_nodes_new, '\n',  TR_SS)
        # 
        # print(b/a)
        v1_3.append(c/a)
        
        ## Comparaison entre version initiale et nouvelle version avec enumerate et un seul tuyau avec complexite en espace moindre
        
        SS_A = [HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000)]
        SS_B = [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350)]
        SS_C = [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500)]
        
        liste_substations = [SS_A, SS_B, SS_C]
        
        Tr_nodes = [X[1].TR_ext() for X in liste_substations]
        t1 = time.time()
        for k in range(1000):
            Tr_node_network_upstream = 0
            m_dotR = 0 
            TR_SS= []
        
            for i, (h, p1) in enumerate(liste_substations[::-1]):
                #Calculation of the return nodes temperature in the network at time (t+1)
                
                m_dotSSr = h.m_dot1
                Tr_SS = h.Tr1
                
                TR_SS.append(Tr_SS)
                p1.evolR_T(m_dotR+m_dotSSr, Tr_nodes[i])
                Tr_nodes[i]= (m_dotSSr*Tr_SS + m_dotR * Tr_node_network_upstream)/(m_dotSSr + m_dotR)
                m_dotR += m_dotSSr
                Tr_node_network_upstream = p1.TR_ext()
            
        t2 = time.time()
        b = (t2-t1)
        
        # print('________')
        # print((t2-t1)/1000)
        # print(Tr_node_network_upstream,'\n', Tr_nodes_new, '\n',  TR_SS)
        # 
        # print(b/a)
        v1_31.append(b/c)
    
    print(np.sum(v1_11)/100)
    print(np.sum(v1_12)/100)
    print(np.sum(v1_3)/100)
    print(np.sum(v1_31)/100)
    
##TEST 2

def f2():
    v1_11 = []
    v1_12 = []

    for jas in range(100):
        
        SS_A = [HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000)]
        SS_B = [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350)]
        SS_C = [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500)]
        
        liste_substations = [SS_A, SS_B, SS_C]
        Ts_nodes = [X[1].TS_ext() for X in liste_substations]
        subm_dot = [X[0].m_dot1 for X in liste_substations] #substations mass flows at time t+1
        m_dot = sum(subm_dot)
        
        t1 = time.time()
        
        for p in range(1000):
            
            T_node = [70] + list(Ts_nodes)
            for i, (hex, pipe1) in enumerate(liste_substations):
                #Calculation of Temperatures in the network at time (t+1) (for next iteration)
                pipe1.evolS_T(m_dot, T_node[i])
                
                #Calculation of the new mass flow and return temperature at substation for time t+1
                m_dotSS = hex.m_dot1
                hex.solve(T_node[i+1])
                m_dot -= m_dotSS
                
            Ts_nodes = [X[1].TS_ext() for X in liste_substations]
            subm_dot = [X[0].m_dot1 for X in liste_substations] #substations mass flows at time t+1
            m_dot = sum(subm_dot)
            
        t2 = time.time()
        a = (t2-t1)
        
        
        
        ##
        
        SS_A = [HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000)]
        SS_B = [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350)]
        SS_C = [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500)]
        
        liste_substations = [SS_A, SS_B, SS_C]
        Ts_nodes = [X[1].TS_ext() for X in liste_substations]
        subm_dot = [X[0].m_dot1 for X in liste_substations] #substations mass flows at time t+1
        mdot_new = sum(subm_dot)
        
        t1 = time.time()
        for p in range(1000):
            m_dot = mdot_new
            mdot_new = 0
            T_node = [70] + list(Ts_nodes)
            
            for i, (hex, pipe1) in enumerate(liste_substations):
                
                #Calculation of Temperatures in the network at time (t+1) (for next iteration)
                pipe1.evolS_T(m_dot, T_node[i])
                Ts_nodes[i] = pipe1.TS_ext()
                
                #Calculation of the new mass flow and return temperature at substation for time t+1
                m_dotSS = hex.m_dot1
                hex.solve(T_node[i+1])
                m_dot -= m_dotSS
                
                subm_dot[i] = hex.m_dot1
                mdot_new += hex.m_dot1
            
        t2 = time.time()
        b = (t2-t1)
        
        v1_11.append(b/a)
        
    print(np.sum(v1_11)/100)
        
    