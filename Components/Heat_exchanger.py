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
        
        self.UA_backup = 4000
    
    def UA(self):
        '''calculates the UA in Q = UA.deltaTlog, with the approximation given in J.J.J. Chen, Comments on improvements on a replacement for the logarithmic mean  '''
        a = 0.3275
        TS2 = self.Ts2_vrai #min(self.Ts2_vrai,self.Ts2)
        Q = Cp * self.m_dot2 * (TS2 - self.Tr2)
        deltaTlm = ((self.Ts1-TS2) - (self.Tr1-self.Tr2))/ln((self.Ts1 - TS2)/(self.Tr1-self.Tr2))
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
            # print(UA)
            # print(Q/(Ts1 - Ts2))
            Tr1 = Tr2 + (2 * (Q/(UA))**a - (Ts1 - Ts2)**a)**(1/a)
            m1 = m2*(Ts2 - Tr2)/(Ts1 - Tr1)
            
            #if approximation does not give an intended result
            if math.isnan(Tr1):
                UA = self.UA_backup
                Tr1 = Tr2 + (2 * (Q/(UA))**a - (Ts1 - Ts2)**a)**(1/a)
                m1 = m2*(Ts2 - Tr2)/(Ts1 - Tr1)
            # if math.isnan(Tr1):

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
            #     print(x1, f(x1))
            #     Tr1 = Tr2 + x1 * (Ts1-Ts2)
            #     m1 = Q/(Cp * (Ts1 - Tr1))
            #     Ts2_vrai = Ts2
            
            # #if approximation does not give an intended result
            # if math.isnan(Tr1):
            #     q = 0.8
            #     kA = UA*(self.m_dot1**(-q) + m2**(-q))
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
            #         x1 = 0
            #     Tr1 = Tr2 + x1 * (Ts1-Ts2)
            #     m1 = Q /(Cp * (Ts1 - Tr1))
                
            Ts2c = (Tr1 > Tr2)*(Tr1 < Ts1)
            Ts2_vrai = Ts2
                
        if Ts2>Ts1 or not Ts2c or m1 > self.qmax:
            #print('Uncovered heat')
            m1 = self.qmax
            NUT = UA/(Cp*m2)
            R = m2/m1
            E = eff(NUT, R)
            
            #Artificial safety:
            if E > 0.50:
                print(E)
                E = 0.42
                
            Ts2_vrai = Tr2 + E*(Ts1-Tr2)
            Tr1 = Ts1 - (m2/m1)*(Ts2_vrai-Tr2)

        self.Ts1 = Ts1
        self.Tr1 = Tr1
        self.m_dot1 = m1
        self.Ts2_vrai = Ts2_vrai
        
        if Ts2_vrai <= Ts2:
            self.UA_backup = self.UA()