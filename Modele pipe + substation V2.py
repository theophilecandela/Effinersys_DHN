import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
from random import gauss
from itertools import *

## Datas
#Pipe parameters
#ambient temperature
Ta = 283.15
#Thermal conductivity of pipe isolation (W/m/K)
lam_i = 0.025
#Thermal conductivity of steel (W/m/K)
lam_p = 50.2
#Thermal conductivity of the ground (concrete) (W/m/K) 
lam_s = 1.6
#Thermal conductivity of water (W/m/K)
lam_e = 0.6071
#Thermal capacity of water (J/K/kg)
Cp = 4185
#Dynamic viscosity of water (Pa.s)
mu_e = 0.535e-3
#volumetric mass density of water (kg/m^3)
rho = 1000 
#Depth of buried pipe (m)
z = 2
#Pipe radius (m)
R_int = 200e-3
R_p = 250e-3
R_i = 400e-3

#Numerical parameters
dx = 90  #spatial discretization step (m)
dt = 60  #discretization time step (s)

##Pipeline
def step(L, dx0=dx):
    '''determine the most convenient spatial step for discretization, such as dx <90m and L%dx = 0'''
    nb_controlvolume = np.ceil(L/dx0)
    dx_new = L/ nb_controlvolume
    return int(nb_controlvolume), dx_new
    
class Pipe:
    #Every pipe is double pipe, with both a supply and a return pipes
    def __init__(self, lam_i, lam_p, lam_s, R_int, R_p, R_i, z, L):
        self.param = {'lam_i': lam_i,
        'lam_p': lam_p,
        'lam_s': lam_s,
        'R_int' : R_int,
        'R_p': R_p,
        'R_i': R_i,
        'prof': z,
        'Rth': (ln(R_p/R_int)/lam_p + ln(R_i/R_p)/lam_i + ln(2*z/R_i)/lam_s)/(2*np.pi),
        'S': np.pi * (R_int)**2}
        
        self.k = 0.027 * (lam_e **0.67) * (Cp ** 0.33) / ((2*(R_int) **0.2) * (mu_e**0.47) * ((np.pi * (R_int)**2) ** 0.8))
        
        self.length = L
        self.nb_controlvolume, self.dx = step(self.length)
        #Arbitratry temperatures initialization
        self.pipeS_T = [343.15]*self.nb_controlvolume
        self.pipeR_T = [323.15]*self.nb_controlvolume
        
    def R(self, m_dot):
        '''calculates the thermal resistance of the pipe, taking convection into consideration, as a function of mass flow only'''
        H_n = self.k * m_dot**0.8
        R = self.param['Rth'] + 1/(2*np.pi*H_n)
        return R
    
    def evolS_T(self, m_dot, T_in): 
        ''' Calculates the evolution of temperatures in the supply pipe during on minute with entrance temperature = T_in.
        Temperatures must be given in °C '''
        T_in = T_in + 273.15
        A = self.param['S']
        T_n = np.copy(self.pipeS_T)
        k = 0
        R = self.R(m_dot)
        for i, Ti in enumerate(T_n):
            if i == 0:
                T_n[0] = (T_n[0] + m_dot*dt/dx *T_in + dt/(A * rho * Cp * R)*Ta)/(1 + m_dot*dt/dx + dt/(A * rho * Cp * R))
            else:
                T_n[i] = (Ti + m_dot*dt/dx *T_n[i-1] + dt/(A * rho * Cp * R)*Ta)/(1 + m_dot*dt/dx + dt/(A * rho * Cp * R))
            
        self.pipeS_T = T_n
        
    def evolR_T(self, m_dot, Tr_in):
        ''' Calculates the evolution of temperatures in the return pipe during on minute with entrance temperature = T_in.
        Temperatures must be given in °C'''
        T_in = Tr_in + 273.15
        A = self.param['S']
        T_n = np.copy(self.pipeR_T)
        k = 0
        R = self.R(m_dot)
        for i, Ti in enumerate(T_n):
            if i == 0:
                T_n[0] = (T_n[0] + m_dot*dt/dx *T_in + dt/(A * rho * Cp * R)*Ta)/(1 + m_dot*dt/dx + dt/(A * rho * Cp * R))
            else:
                T_n[i] = (Ti + m_dot*dt/dx *T_n[i-1] + dt/(A * rho * Cp * R)*Ta)/(1 + m_dot*dt/dx + dt/(A * rho * Cp * R))
            
        self.pipeR_T = T_n
    
    def TS_ext(self):
        '''return downstream outlet temperature of the supply pipe in °C'''
        return (self.pipeS_T[-1]-273.15)
    
    def TR_ext(self):
        '''return downstream outlet temperature of the return pipe in °C'''
        return (self.pipeR_T[-1]-273.15)

def test_evol_pipe(pipe, T_in, mdot, n):
    k = 0
    plt.plot(list(range(pipe.nb_controlvolume)), pipe.pipeS_T)
    plt.show()
    while k < n :
        pipe.evolS_T(mdot, T_in)
        if k%5 == 0:
            plt.plot(list(range(pipe.nb_controlvolume)), pipe.pipeS_T)
            plt.show()
        k += 1

        
##Heat exchanger
# def UA(Ts1, Tr1, Ts2, Tr2, m2):
#     a = 0.3275
#     Q = Cp * m2 * (Ts2 - Tr2)
#     deltaTlm = ((Ts1-Tr1) - (Ts2-Tr2))/ln((Ts1 - Tr1)/(Ts2-Tr2))
#     UA = Q/deltaTlm
#     return UA

def eff(N, R):
    return (1 - np.exp((R-1)*N))/(1-R*np.exp((R-1)*N))
    
def newton(f, f_prime, x0, eps = 0.0001):
    x = x0
    n = 15
    k = 0
    while np.abs(f(x)) > eps and k < n:
        #print(f(x), f_prime(x), x)
        x = max(0.0001, x - f(x)/f_prime(x))
        k += 1
    if k>= n:
        raise ConvergenceError('Ne converge pas')
    else:
        return x

class HEX: 
    #U = constant method
    def __init__(self, Ts1_init, Ts2, Tr1_init, Tr2, m_dot2, qmax):
        ''' Ts1 : inlet temperature primary side
            Tr1: outlet temperature primary side
            Tr2: inlet temperature secondary side
            Ts2: output command temperature secondary side
            m_dot2: mass flow rate secondary side'''
        self.Ts1 = Ts1_init
        self.Tr1 = Tr1_init
        self.Ts2 = Ts2
        self.Ts2_vrai = Ts2 #Ts2_vrai == Ts2 unless uncovered heat 
        self.Tr2 = Tr2
        self.m_dot2 = m_dot2
        self.m_dot1 = 1.1
        self.qmax = qmax   #maximal mass flow going through the heat exchanger
    
    def UA(self):
        '''calculates the UA in Q = UA.deltaTlog  '''
        a = 0.3275
        Q = Cp * self.m_dot2 * (self.Ts2_vrai - self.Tr2)
        #deltaTlm = ((self.Ts1-self.Tr1) - (self.Ts2_vrai-self.Tr2))/ln((self.Ts1 - self.Tr1)/(self.Ts2_vrai-self.Tr2))
        deltaTlm = ((self.Ts1-self.Ts2_vrai) - (self.Tr1-self.Tr2))/ln((self.Ts1 - self.Ts2_vrai)/(self.Tr1-self.Tr2))
        UA = Q/deltaTlm
        return UA
        
    def solve(self, Ts1):
        '''For a primary side supply temperature, calculates the return temperature and primary side mass flow, as well as determining whether the network can cover the heat demand'''
        a = 0.3275
        Ts2c = True
        UA = self.UA() #Attention si on modifie Ts2 entre l'itération prédente et la nouvelle?
        Ts2 = self.Ts2
        Tr2 = self.Tr2
        m2 = self.m_dot2
        Q = Cp * self.m_dot2 * (Ts2 - Tr2)
        m1 = self.qmax
        
        if Ts2 <= Ts1:
            Tr1 = Tr2 + (2 * (Q/(UA))**a - (Ts1 - Ts2)**a)**(1/a)
            m1 = m2*(Ts2 - Tr2)/(Ts1 - Tr1)
            Ts2c = (Tr1 > Tr2)*(Tr1 < Ts1)
            Ts2_vrai = Ts2
        
        if Ts2>Ts1 or not Ts2c or m1 > self.qmax:
            #print('Uncovered heat')
                
            m1 = self.qmax
            NUT = UA/(Cp*m2)
            R = m2/m1
            E = eff(NUT, R)
            #print(E)
            Ts2_vrai = Tr2 + E*(Ts1-Tr2)
            Tr1 = Ts1 - (m2/m1)*(Ts2_vrai-Tr2)
                
        # #Case in which the approximation for DeltaTlm doesn't work
        # if Ts2_vrai > Ts2: 
        #     #We use another method
        #     q = 0.8
        #     k = Q/(UA* (Ts1 - Ts2))
        #     
        #     def f(x):
        #         return np.exp(k*(x-1)) - x

        #     def f_prime(x):
        #         return k*np.exp(k*(x-1)) - 1
        #             
        #     x0 = (self.Tr1 - Tr2)/(Ts1 -Ts2)
        #     x1 = newton(f, f_prime, x0)
        #     #print(x1, f(x1))
        #     Tr1 = Tr2 + x1 * (Ts1-Ts2)
        #     m1 = Q/(Cp * (Ts1 - Tr1))
        #     Ts2_vrai = Ts2
            
        self.Ts1 = Ts1
        self.Tr1 = Tr1
        self.m_dot1 = m1
        self.Ts2_vrai = Ts2_vrai
        
## Réseau

class Network:
    def __init__(self, inlet_T, liste_substations):
        ''' inlet_T ) initialization supply Temperature of the network
        liste_substations : [HEX, Pipe_nodetopreviousnode, Pipe_hextonode]'''
        self.supplyT = inlet_T
        self.returnT = 0
        self.subm_dot = [X[0].m_dot1 for X in liste_substations]
        self.m_dot = sum([X[0].m_dot1 for X in liste_substations])
        self.Ts_nodes = [X[1].TS_ext() for X in liste_substations]
        self.Tr_nodes = [X[1].TR_ext() for X in liste_substations]
        self.nb_substations =  len(liste_substations)
        self.substations = liste_substations
        
        
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
        

## Simulation
dx = 90  # spatial discretization step (m)
dt = 60  # discretization time step (s)
#Inlet temperature of the network
Tint = [gauss(60,1) for i in range(24)]

def simulation( RES, Ts2_h, Tsupply_h = Tint):
    '''Ts2_h:  list of list of temperature per hour (ie scheduled demand for each substation)
    Tsupply_h: list of source water temperature thoughout the day (per hour)'''
    nb_SS = len(Ts2_h)
    t = []
    T_in = []
    t_tot = 0
    m_dot = []
    T_return = []
    T_return_SS = []
    T_supplySS = []
    mdot_SS = []
    T_supply_secondary = []
    T_secondary_demand = []
    for j in range(24):
        RES.supplyT = Tsupply_h[j]
        print(j)
        for p, T_h in enumerate(Ts2_h):
            RES.substations[p][0].Ts2 = T_h[j]
        
        for m in range(60):
            RES.iteration()
            #print(j, Tsupply_h[j], RES.m_dot, RES.substations[0][0].Ts1, RES.substations[0][0].Tr1, RES.substations[0][0].Ts2_vrai)
            t.append(t_tot)
            T_in.append(RES.supplyT)
            T_return.append(RES.returnT)
            T_return_SS.append([X[0].Tr1 for X in RES.substations])
            T_supplySS.append([X[0].Ts1 for X in RES.substations])
            mdot_SS.append([X[0].m_dot1 for X in RES.substations])
            m_dot.append(RES.m_dot)
            T_secondary_demand.append([X[0].Ts2 for X in RES.substations])
            T_supply_secondary.append([X[0].Ts2_vrai for X in RES.substations])
            t_tot += 1
            
    
    plt.figure()
    plt.plot(t, T_return, label = 'Network return temperature at the source')
    for i in range(len(T_return_SS[0])):
        plt.plot(t, [a[i] for a in T_return_SS], label= f'return T substation{i}')
    plt.legend()
    
    plt.figure()
    plt.plot(t, m_dot, label = 'total Network mass flow')
    for i in range(len(mdot_SS[0])):
        plt.plot(t, [a[i] for a in mdot_SS], label= f'mass flow {i}')
    plt.legend()
    
    plt.figure()
    plt.plot(t, T_in, label = 'supply temperature at source')
    for i in range(len(T_supplySS[0])):
        plt.plot(t, [a[i] for a in T_supplySS], label = f'T_supply SS {i}')
    plt.legend()
    
    for i in range(len(T_supply_secondary[0])):
        plt.figure()
        plt.plot(t, [a[i] for a in T_supply_secondary], label = f'T_supply_secondary {i}')
        plt.plot(t, [a[i] for a in T_secondary_demand], label = 'T_secondary_demand')
        plt.plot(t, [a[i] for a in T_supplySS], label = 'T_supply SS network side')
        plt.plot(t, [a[i] for a in T_return_SS], label = 'T_return SS network side')
        plt.legend()   
        # plt.figure()
        # plt.plot(t, [np.abs(a[i] - b[i]) for a, b in zip(T_supply_secondary,T_secondary_demand) ], label = 'Secondary supply default')
        plt.legend() 
    
    plt.show()
    return t, mdot_SS
    
#Inlet temperature of the network
Tint = [gauss(75,4) for i in range(24)]
#Temperature demand (profiles)
Ts2_1 = [41, 38, 37, 37, 40, 47, 52, 49, 45, 43, 42, 43, 46, 45, 41, 42, 42, 43, 50, 49, 48, 46, 45, 42]

Ts2_2 = [42, 42.5, 43, 43.5, 44, 48.5, 53, 53, 52, 52, 50, 49, 48,47, 46.5, 46, 49, 49.5, 48, 49, 49.5, 50, 50, 42]
Ts2_2bis = [x - 3 for x in Ts2_2]

Ts2_3 = [41, 38, 37, 37, 40, 44, 48, 47, 40, 39, 39, 38.5, 39, 39.5, 39, 38.5, 39.5, 39, 50, 49.5, 49, 46, 46, 42]

noise = [gauss(0,2) for i in range(24)]
t = [i for i in range(24)]

Tr2 = 30
m_dot2 = 1.6
    
    
SS1 = [HEX(70, 45, 44, 30, 1.6, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
SS2 = [HEX(70, 45, 44, 25, 1.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]
SS3 = [HEX(70, 45, 44, 28, 1.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 2500), Pipe(lam_i, lam_p, lam_s, 120e-3, 145e-3, 195e-3, z, 50)]

NET = Network(75, [SS1])
NET2 = Network(75, [SS1, SS2])
NET3 = Network(75, [SS1, SS2, SS3])


#simulation(NET3, [Ts2_1, Ts2_2, Ts2_3])
#simulation(NET3, [Ts2_1, Ts2_2bis, Ts2_3])








        
