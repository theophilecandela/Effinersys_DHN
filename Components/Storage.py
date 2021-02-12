import numpy as np
from numpy import log as ln
import math
from itertools import *
from Components.Ressources import *

class Buffer:
    def __init__(self, hT, lT, hV, lV):
        self.hot_T = hT
        self.low_T = lT
        self.hot_V = hV
        self.low_V = lV
        self.m_dot = 0
        self.P_in_hot = 0
        self.P_in_cold = 0
        
    def intake_hot_water(self, T):
        mdot = np.abs(self.m_dot)
        V_new = self.hot_V + mdot * dt /rho
        T_new = (self.hot_V*self.hot_T + mdot*T * dt/rho)/V_new
        
        self.hot_V = V_new
        self.hot_T = T_new
        self.P_in_hot = mdot*T*Cp
      
        
    def intake_cold_water(self, mdot, T):
        mdot = np.abs(mdot)
        if self.hot_V - mdot * dt/rho >= 0:
            V_new = self.low_V + mdot * dt /rho
            T_new = (self.low_V*self.low_T + mdot*T * dt/rho)/V_new
            
            self.low_V = V_new
            self.low_T = T_new
            self.P_in_cold = mdot*T*Cp
            return mdot
            
        else:
            self.P_in_cold = 0
            return 0
        
    
    
    def delivery_hot_water(self):
        mdot = np.abs(self.m_dot)
        self.hot_V = self.hot_V - mdot * dt/rho
        self.P_in_hot = - mdot*self.hot_T*Cp
        return self.hot_T, mdot
        
       
        
    def delivery_cold_water(self, mdot):
        mdot = np.abs(mdot)
        if self.low_V - mdot * dt/rho > 0:
            self.low_V = self.low_V - mdot * dt/rho
            self.P_in_cold = -mdot*self.low_T*Cp
            return self.low_T, mdot
        else:
            self.P_in_cold = 0
            return self.low_T, 0

# There is no possibility to pump and inject in the same compartment at the same time