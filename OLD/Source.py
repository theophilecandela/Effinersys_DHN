import numpy as np
from numpy import log as ln
from itertools import *
from Components.Ressources import *

class Source: 
    def __init__(self, geoT, geoMdot, Q_nom, deltaTlm_nom):
        self.m_dot = geoMdot
        self.Ts_Geo= geoT
        self.P = Q_nom
        self.Tr_Geo = geoT - Q_nom/(geoMdot * Cp)
        self.Ts_Net = geoT - 2
        #self.Tr_Net = 40
        if self.Ts_Geo - self.Ts_Net > deltaTlm_nom:
            x0 = deltaTlm_nom/2
        elif self.Ts_Geo - self.Ts_Net < deltaTlm_nom:
            x0 = 2*deltaTlm_nom
        f = lambda x: (x - deltaTlm_nom*ln(x) - (self.Ts_Geo - self.Ts_Net - deltaTlm_nom*ln(self.Ts_Geo - self.Ts_Net)))
        f_prime = lambda x: 1 - deltaTlm_nom/x
        self.Tr_Net = self.Tr_Geo - newton(f, f_prime, x0)
        
    
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
        #print(f'UA = {UA}')
        self.Tr_Net  = TrNET
        if m_dotNET < self.m_dot:
            NUT = UA/(Cp*m_dotNET)
            R = m_dotNET/self.m_dot
            E = eff(NUT, R)
            print(f'E = {E}')
            self.Ts_Net = self.Tr_Net + E*(self.Ts_Geo - self.Tr_Net)
            self.Tr_Geo = self.Ts_Geo - R*(self.Ts_Net-self.Tr_Net)
        
        else:
            NUT = UA/(Cp*self.m_dot)
            R = self.m_dot/m_dotNET
            E = eff(NUT, R)
            self.Tr_Geo = self.Ts_Geo - E*(self.Ts_Geo - self.Tr_Net)
            self.Ts_Net = self.Tr_Net + R*(self.Ts_Geo - self.Tr_Geo)
        
        self.P = (self.Ts_Geo - self.Tr_Geo)*self.m_dot*Cp
        