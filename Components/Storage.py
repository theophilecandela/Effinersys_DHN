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
        
    def intake_hot_water(self, T):
        #appel à delivery_cold_water --> no because we have to deal with it in the previous iteration
        mdot = np.abs(self.m_dot)
        V_new = self.hot_V + mdot * dt /rho
        T_new = (self.hot_V*self.hot_T + mdot*T * dt/rho)/V_new
        
        self.hot_V = V_new
        self.hot_T = T_new
        
    def intake_cold_water(self, mdot, T):
        # appel à delivery_hot_water --> no because we have to deal with it in the next iteration
        V_new = self.low_V + mdot * dt /rho
        T_new = (self.low_V*self.low_T + mdot*T * dt/rho)/V_new
        
        self.low_V = V_new
        self.low_T = T_new
    
    def delivery_hot_water(self):
        mdot = np.abs(self.m_dot)
        if self.hot_V > 0:
            self.hot_V = self.hot_V - mdot * dt/rho
            return self.hot_T, mdot
        else:
            self.hot_V = 0
            self.low_V = self.low_V - mdot * dt/rho
            return self.low_T, mdot
        
    def delivery_cold_water(self, mdot):
        if self.low_V > 0:
            self.low_V = self.low_V - mdot * dt/rho
            return self.low_T
        else:
            self.low_V  = 0
            self.hot_V = self.hot_V - mdot * dt/rho
            return self.hot_T

# There is no possibility to pump and inject in the same compartment at the same time