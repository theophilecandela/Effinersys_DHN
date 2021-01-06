import numpy as np
from numpy import log as ln
from itertools import *
from Components.Ressources import *

class Source: 
    def __init__(self, geoT, geoMdot):
        self.m_dot = geoMdot
        self.Ts_Geo= geoT
        self.Tr_Geo = 55
        self.Ts_Net = 70
        self.Tr_Net = 40
    
    def UA(self):
        '''calculates the UA in Q = UA.deltaTlog, with the approximation given in J.J.J. Chen, Comments on improvements on a replacement for the logarithmic mean  '''
        a = 0.3275
        Q = Cp * self.m_dot * (self.Ts_Geo - self.Tr_Geo)
        deltaTlm = ((self.Ts_Geo-self.Tr_Geo) - (self.Ts_Net-self.Tr_Net))/ln((self.Ts_Geo - self.Tr_Geo)/(self.Ts_Net-self.Tr_Net))
        UA = Q/deltaTlm
        return UA
        
    def solve(self, m_dotNET, TrNET):
        '''For a secondary side (network side) return temperature and mass flow, calculates the supply temperature of the network'''
        a = 0.3275
        UA = self.UA() 
        if m_dotNET < self.m_dot:
            NUT = UA/(Cp*m_dotNET)
            R = m_dotNET/self.m_dot
            E = eff(NUT, R)
            self.Ts_Net = self.Tr_Net + E*(self.Ts_Geo - self.Tr_Net)
            self.Tr_Geo = self.Ts_Geo - R*(self.Ts_Net-self.Tr_Net)
        
        else:
            NUT = UA/(Cp*self.m_dot)
            R = self.m_dot/m_dotNET
            E = eff(NUT, R)
            self.Tr_Geo = self.Ts_Geo - E*(self.Ts_Geo - self.Tr_Net)
            self.Ts_Net = self.Tr_Net + R*(self.Ts_Geo - self.Tr_Geo)
            