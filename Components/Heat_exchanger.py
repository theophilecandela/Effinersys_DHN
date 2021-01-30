import numpy as np
from numpy import log as ln
import math
from itertools import *
from Components.Ressources import *

class HEX: 
    #(U = constant) method
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
        '''calculates the UA in Q = UA.deltaTlog, with the approximation given in J.J.J. Chen, Comments on improvements on a replacement for the logarithmic mean  '''
        a = 0.3275
        TS2 = self.Ts2_vrai #min(self.Ts2_vrai,self.Ts2)
        Q = Cp * self.m_dot2 * (TS2 - self.Tr2)
        deltaTlm = ((1/2)*((self.Ts1-TS2)**a + (self.Tr1-self.Tr2)**a))**(1/a)
        #deltaTlm_vrai = ((self.Ts1-TS2) - (self.Tr1-self.Tr2))/ln((self.Ts1 - TS2)/(self.Tr1-self.Tr2))
        UA = Q/deltaTlm
        
        #print(Q,deltaTlm, self.m_dot1, UA)
        return UA
        
    def solve(self, Ts1):
        '''For a primary side supply temperature, calculates the return temperature and primary side mass flow, as well as determining whether the network can cover the heat demand'''
        a = 0.3275
        Ts2c = True
        UA = self.UA() 
        Ts2 = self.Ts2
        Tr2 = self.Tr2
        m2 = self.m_dot2
        Q = Cp * self.m_dot2 * (Ts2 - Tr2)
        m1 = self.qmax
        if Ts2 <= Ts1:
            Tr1 = Tr2 + (2 * (Q/(UA))**a - (Ts1 - Ts2)**a)**(1/a)
            m1 = m2*(Ts2 - Tr2)/(Ts1 - Tr1)
            
            #if approximation does not give an intended result
            if math.isnan(Tr1):
                raise ValueError(((Ts1 - Ts2)/(Q/(UA)))/(15/3.5))

            Ts2c = (Tr1 > Tr2)*(Tr1 < Ts1)
        Ts2_vrai = Ts2
                
        if Ts2>Ts1 or not Ts2c or m1 > self.qmax:
            m1 = self.qmax
            NUT = UA/(Cp*m2)
            R = m2/m1
            E = eff(NUT, R)
            Ts2_vrai = Tr2 + E*(Ts1-Tr2)
            Tr1 = Ts1 - (m2/m1)*(Ts2_vrai-Tr2)

        self.Tr1 = Tr1
        self.m_dot1 = m1
        self.Ts2_vrai = Ts2_vrai
        self.Ts1 = Ts1
        
class HEX_nom(HEX):
    def __init__(self, Qnom, Tr2, f_Ts2, hnom, deltaTlm_nom , qmax):
        mdot_2 = Qnom/(Cp * (f_Ts2(-7) - Tr2))
        UAnom = Qnom/deltaTlm_nom
        Unom = hnom/2
        self.Aex = UAnom/Unom
        self.hnom = hnom
        self.Ts1 = f_Ts2(-7) + 3
        
        self.f_Ts2 = f_Ts2
        self.Ts2 = f_Ts2(-7)
        self.Ts2_vrai = f_Ts2(-7) #Ts2_vrai == Ts2 unless uncovered heat 
        self.Tr2 = Tr2
        self.m_dot2 = mdot_2
        self.m_dot1 = mdot_2
        self.qmax = qmax 
        
        #Why? We do know Q, m_dot and Ts1
        if self.Ts1 - self.Ts2 > deltaTlm_nom:
            x0 = deltaTlm_nom/2
        elif self.Ts1 - self.Ts2 < deltaTlm_nom:
            x0 = 2*deltaTlm_nom
        f = lambda x: (x - deltaTlm_nom*ln(x) - (self.Ts1 - self.Ts2 - deltaTlm_nom*ln(self.Ts1 - self.Ts2)))
        f_prime = lambda x: 1 - deltaTlm_nom/x
        self.Tr1 = Tr2 + newton(f, f_prime, x0)
        
    def UA(self):
        '''calculates the UA in Q = UA.deltaTlog, with the approximation given in J.J.J. Chen, Comments on improvements on a replacement for the logarithmic mean  '''
        U = self.hnom*self.m_dot2**0.8/((self.m_dot2/self.m_dot1)**0.8 + 1)
        UA = U*self.Aex
        return UA
        
        