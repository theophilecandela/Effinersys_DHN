import numpy as np
from numpy import log as ln
from itertools import *
from Components.Ressources import *

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
        
##Test

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
