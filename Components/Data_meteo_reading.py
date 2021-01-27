import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
import pickle

#File containing outside Temperature measure, every ten minutes
with open(f'Components/Data_meteo.txt', 'rb') as data:
    data.readline()
    data.readline()
    f = data.readlines()

with open(f'Components/Data_meteo_oneday.txt', 'rb') as data:
    data.readline()
    data.readline()
    f = data.readlines()
    
    
    
lignes = [str(l.rstrip()) for l in f if len(l.rstrip()) != 0]
lignes = [l[2:-1:] for l in lignes]
lignes = [l.replace('+', '').replace(' ', '').split('\\t') for l in lignes]
#lignes = [(l[0] + ')').replace('+', '(') for l in lignes]
lignes = [(float(l[0]), float(l[1])) for l in lignes]


Ta = [l[1] for l in lignes[0::6]]
#plt.plot([i for i in range(len(Ta))], Ta)
#Ta_oneday = Ta[240:265]
# plt.plot([i for i in range(len(Ta_oneday))], Ta_oneday)
# plt.show()

def Ts2_A(Ta):
    return 90 - (30/27)* (Ta + 7)

def Ts2_B(Ta):
    return 85 - (23/27)* (Ta + 7)
    
def Ts2_C(Ta):
    return 87 - (15/27)* (Ta + 7)
    
def Ts2_D(Ta):
    return 78 - (8/27)* (Ta + 7)
    
Ts2_A = np.vectorize(Ts2_A)
Ts2_B = np.vectorize(Ts2_B)
Ts2_C = np.vectorize(Ts2_C)
Ts2_D = np.vectorize(Ts2_D)

Tr2_A = 42
Tr2_2 = 39
Tr2_3 = 36
Tr2_4 = 43

def test():
    T = list(range(-7, 21))
    plt.figure()
    plt.plot(T, Ts2_A(T), label = 'Ts2_A')
    plt.plot(T, Ts2_B(T), label = 'Ts2_B')
    plt.plot(T, Ts2_C(T), label = 'Ts2_C')
    plt.plot(T, Ts2_D(T), label = 'Ts2_D')
    plt.legend()
    plt.show()