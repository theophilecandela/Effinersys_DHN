import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
import pickle
import random as rd

#File containing outside Temperature measure, every ten minutes
with open(f'Components/Data_meteo.txt', 'rb') as data:
    data.readline()
    data.readline()
    f = data.readlines()

# with open(f'Components/Data_meteo_oneday.txt', 'rb') as data:
#     data.readline()
#     data.readline()
#     f = data.readlines()
    
    
    
lignes = [str(l.rstrip()) for l in f if len(l.rstrip()) != 0]
lignes = [l[2:-1:] for l in lignes]
lignes = [l.replace('+', '').replace(' ', '').split('\\t') for l in lignes]
#lignes = [(l[0] + ')').replace('+', '(') for l in lignes]
lignes = [(float(l[0]), float(l[1])) for l in lignes]


Ta = [l[1] for l in lignes[0::6]]

#plt.plot([i/24 for i in range(len(Ta))], Ta)
Ta_oneday = Ta[240:264]
# plt.figure()
# plt.plot([i for i in range(len(Ta_oneday))], Ta_oneday)
# plt.show()

Ta_week = Ta[216:396] 
#180 hour = 7days + twelve hour in order to optimise the last twelve hour of the seventh day

# plt.figure()
# plt.plot([i/24 for i in range(len(Ta_week))], Ta_week)
# plt.show()


#parameters = [(90,30), (85, 23), (87,15), (77, 7), (83, 12), (88, 26), (82, 16), (80, 10), (83, 21), (79, 14), (84, 9)]

#parameters = [(rd.randrange(70, 85, 1), rd.randrange(60, 67, 1)) for i in range(10)]
parameters = [(75, 64), (79, 60), (80, 62), (82, 62), (70, 66), (70, 66), (80, 65), (73, 66), (81, 60), (73, 64)]
def functionize(a,b):
    c = a - b
    f= lambda x: a - (c/27) * (x + 7)
    return np.vectorize(f)
    
functions = [functionize(a, b) for a, b in parameters]
Tr2_SS = [rd.randrange(36, 43, 1) for i in range(10)]


def plot():
    T = list(range(-7, 21))
    plt.figure()
    for i, f in enumerate(functions):
        plt.plot(T, f(T), label = f'{i}')
    plt.legend()
    plt.show()




## Useless - Previous garbage
def Ts2_A(Ta):
    return 90 - (30/27)* (Ta + 7)

def Ts2_B(Ta):
    return 85 - (23/27)* (Ta + 7)
    
def Ts2_C(Ta):
    return 87 - (15/27)* (Ta + 7)
    
def Ts2_D(Ta):
    return 78 - (8/27)* (Ta + 7)
    
def Ts2_E(Ta):
    return 83 - (12/27)* (Ta + 7)
    
def Ts2_F(Ta):
    return 88 - (26/27)* (Ta + 7)
    
def Ts2_G(Ta):
    return 82 - (16/27)* (Ta + 7)
    
Ts2_A = np.vectorize(Ts2_A)
Ts2_B = np.vectorize(Ts2_B)
Ts2_C = np.vectorize(Ts2_C)
Ts2_D = np.vectorize(Ts2_D)
Ts2_E = np.vectorize(Ts2_E)
Ts2_F = np.vectorize(Ts2_F)
Ts2_G = np.vectorize(Ts2_G)

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
    plt.plot(T, Ts2_E(T), label = 'Ts2_E')
    plt.plot(T, Ts2_F(T), label = 'Ts2_F')
    plt.plot(T, [A3.f_Ts1(t) for t in T])
    plt.legend()
    plt.show()