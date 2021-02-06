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
    
    def solve_bis(self, Ts1, Ts2, Q, UA):
        raise ValueError(((Ts1 - Ts2)/(Q/(UA)))/(15/3.5))
        return(0,0)
        
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
        
        # P = 50000
        # Tr1 = Ts1 - P/(Cp * self.m_dot1)
        # Ts2_vrai = self.Ts2
        # self.Tr2 = self.Ts2 - P/(Cp * m2)
        if Ts2 <= Ts1:
            Tr1 = Tr2 + (2 * (Q/(UA))**a - (Ts1 - Ts2)**a)**(1/a)
            m1 = m2*(Ts2 - Tr2)/(Ts1 - Tr1)
            
            #if approximation does not give an intended result
            if math.isnan(Tr1):
                Tr1, m1 = self.solve_bis(Ts1, Ts2, Q, UA)

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
        #self.m_dot1 = m1
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
        
    # def solve_bis(self, Ts1, Ts2, Q, UA):
    #     q = 0.8
    #     kA = self.hnom*self.Aex
    #     Tr2 = self.Tr2
    #     Q1 = Q/(Ts1 - Ts2)
    #     Q2 = Q/(Cp * (Ts2 - Tr2))
    #     c1 = kA*(Q2**q)
    #     c2 = (Ts1-Ts2)/(Ts2-Tr2)
    #     
    #     def g1(x):
    #         return x - 1
    #     def g2(x):
    #         return 1 + (1 - c2 * (x-1))**q
    #     def g2_prime(x):
    #         return -c2 * q * (1 - c2 * (x-1))**(q-1)
    #     def g3(x):
    #         return np.log(x)
    #         
    #     def f(x):
    #         if x == 1:
    #             return 1/2 - Q1/c1
    #         elif x == 0:
    #             return - Q1/c1
    #         else:
    #             return g1(x)/(g2(x) * g3(x))  - Q1/c1
    #     
    #     def f_prime(x):
    #         if x == 1:
    #             return 1/2 * (1/2 + c2 * q /2)
    #         elif x == 0:
    #             return 1
    #         else:
    #             return ((f(x)+Q1/c1) * ((g3(x) - g1(x)/x)/(g1(x)*g3(x)) - g2_prime(x)/g2(x)))
    #     
    #     x0 = 1/c2 + 0.8
    #     try:
    #         x1 = newton(f, f_prime, x0)
    #         
    #     except ValueError:
    #         #x1 = 0
    #         raise ValueError(((Ts1 - Ts2)/(Q/(UA)))/(15/3.5))
    #         
    #     Tr1 = Tr2 + x1 * (Ts1-Ts2)
    #     m1 = Q /(Cp * (Ts1 - Tr1))
    #     
    #     return Tr1, m1
     
     
        
class Source(): 

    def __init__(self, geoT, geoMdot, Q_nom, deltaTlm_nom, hnom):
        self.m_dot = geoMdot
        self.Ts_Geo= geoT
        self.P = Q_nom
        self.Tr_Geo = geoT - Q_nom/(geoMdot * Cp)
        self.Ts_Net = geoT - 2
        
        UAnom = Q_nom/deltaTlm_nom
        Unom = hnom/2
        self.Aex = UAnom/Unom
        self.hnom = hnom
        
        if self.Ts_Geo - self.Ts_Net > deltaTlm_nom:
            x0 = deltaTlm_nom/2
        elif self.Ts_Geo - self.Ts_Net < deltaTlm_nom:
            x0 = 2*deltaTlm_nom
        f = lambda x: (x - deltaTlm_nom*ln(x) - (self.Ts_Geo - self.Ts_Net - deltaTlm_nom*ln(self.Ts_Geo - self.Ts_Net)))
        f_prime = lambda x: 1 - deltaTlm_nom/x
        self.Tr_Net = self.Tr_Geo - newton(f, f_prime, x0)
        
    
    def UA(self, m_dotNET ):
        '''calculates the UA in Q = UA.deltaTlog, with the approximation given in J.J.J. Chen, Comments on improvements on a replacement for the logarithmic mean  '''
        U = self.hnom/(self.m_dot**(-0.8) + m_dotNET**(-0.8))
        UA = U*self.Aex
        return UA
        
    def solve(self, m_dotNET, TrNET):
        '''For a secondary side (network side) return temperature and mass flow, calculates the supply temperature of the network'''
        a = 0.3275
        UA = self.UA(m_dotNET) 
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
        
        self.P = (self.Ts_Geo - self.Tr_Geo)*self.m_dot*Cp