import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
from random import gauss, randrange
from itertools import *
import time

#from GA_algorithm import geneticalgorithm as ga
from geneticalgorithm import geneticalgorithm as ga

from Components.Data_meteo_reading import Ts2_A, Ts2_B, Ts2_C, Ts2_D, Ta_twodays
from Components.Data_meteo_reading import Ta_oneday as ext_T
from Components.Data_meteo_reading import Ta_week, functions,fun1_highT, Tr2_SS, functions1, Tr2_SS1
from Components.Networks import Network
from Components.Source import Source
from Components.Heat_exchanger import HEX_nom #, Source
from Components.Pipes import Pipe
from Components.Ressources import *
from Components.Storage import Buffer
from Model_SIM import Simulation

ext_T = [ 3.65, 3.3, 2.85, 2.2, 1.6, 0.517, 0.317,  -0.183, -0.675, -0.883, -0.35, 1.808, 3.65, 5.542, 6.433, 7.325, 7.217, 6.283, 5.65, 4.958, 4.258, 3.65, 2.958, 2.35]

ext_T_twodays = [ 3.65, 2.958, 2.35, 1.8,  0.517, 0.317,  -0.183, -0.383, -0.675, -0.883, -0.35, 1.808, 3.65, 5.542, 6.433, 7.325, 8.3, 7.217, 6.8, 6.283, 5.65, 4.958, 4.258, 3.65, 3.65, 2.958, 2.35, 1.8,  0.517, 0.317,  -0.183, -0.383, -0.675, -0.883, -0.35, 1.808, 3.65, 5.542, 6.433, 7.325, 7.217, 6.283, 5.65, 4.958, 4.258, 3.65, 2.958, 2.35]
##OPTIM Function
algorithm_param={'max_num_iteration': 60, 'population_size': 100, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}

def optim(dim, MOD, param = algorithm_param, step = 1, Storage_dim = None, plot = False):
    '''Uses Genetic algorithm to calculate an optimal/sub-optimal list of instructions;
    dim: number of instruction throughout the considered period of time. From now on, cannot be different than MOD.nb_hour 
    MOD: an object of class Simulation_Ta
    param: a dict containing the parameters of the optimization process
    '''
    time1 = time.time()
    MOD.refined_Ts1()
    try:
        a = MOD.objective_function_Ta()
    except ValueError:
        print('ValueError')
        a = 20
    
    if a > 1.2:
        print(a)
        param = {'max_num_iteration': 60, 'population_size': 110, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}
        if a > 2:
            param = {'max_num_iteration': 60, 'population_size': 130, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}
    
    if Storage_dim is None:
        varbound=np.array([[60, 105]]*dim)
        var_type = None
        if dim == len(MOD.Ta):
            ex = [math.floor(MOD.f_Ts1(T)) for T in MOD.Ta]
        else:
            ex = None
            
    else:
        varbound=np.array([[60, 105]]*dim + [[-5, 5]]*Storage_dim)
        #varbound=np.array([[60, 105]]*dim + ([[-5, 0]]*4 + [[0, 5]]*4)*3)
        var_type = np.array(['int'] * dim + ['real']*Storage_dim)
        if dim == len(MOD.Ta):
            ex = [math.ceil(MOD.f_Ts1(T)) for T in MOD.Ta] + [0]*Storage_dim
        else:
            ex = None
            
        dim += Storage_dim
        
    print(f'ex = {ex} \n')
    
    model=ga(function=MOD.objective_function_optim,dimension= dim,variable_type='int',variable_boundaries=varbound, variable_type_mixed= var_type, algorithm_parameters=param, value_step = step, exemple = ex)
    
    model.run()
    time2 = time.time()
    print(f'[Optimization process] - Time: {time2-time1}')
    
    Instruct = list(model.output_dict['variable'])
    if plot:
        MOD.plot(Instruct[0:MOD.nb_hour])
    return Instruct
       
    
class Result:
    def __init__(self):
        self.state_init = None
        self.T_Boiler_optim = None
        self.storage_flow = None
        self.Tboiler_Ta = None
        self.Ta_week = None 
        self.current_state = None
        
        self.current_Ta = None
    
def optim_week(dim, MOD, Result_class, param = algorithm_param, Ta_w = Ta_week, step_optim_h = 24):
    
    Result_class.Ta_week = np.copy(Ta_week)
    days = [Ta_w[i:i + 24] for i in range(0, 168, step_optim_h)]   #24*7 = 168
    Result_class.T_Boiler_optim = []
    Result_class.storage_flow = []
    Result_class.Tboiler_Ta = []
    MOD.initialisation()
    Result_class.state_init = MOD.save_init_.copy()
    time1 = time.time()
    
    for day in days:
        MOD.Ta = day
        
        Result_class.current_state = MOD.save_init_.copy()
        Result_class.current_Ta = day
        
        Instruct_T_Boiler = optim(dim, MOD, param = {'max_num_iteration': 50, 'population_size': 90, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}, step = 1000)

        Result_class.Tboiler_Ta.extend(MOD.f_Ts1(T)  for T in day[0:step_optim_h])
        
        MOD.initialisation()  #Setting the parameters to init values
        MOD.initialised = False
        MOD.Ta = day[0:step_optim_h]
        
        MOD.initialisation(Instruct_T_Boiler[0:step_optim_h])
        Result_class.T_Boiler_optim.extend(Instruct_T_Boiler[0:step_optim_h])
        
    time2 = time.time()
    print(f'time = {time2 - time1}')

Result_3 = Result()
    
def improve_soluce_storage(soluce, MOD):
    param={'max_num_iteration': 60, 'population_size': 120, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}
    
    varbound=np.array([[60, 105]]*24+ [[-5, 5]]*24)
    
    var_type = np.array(['int'] * 24 + ['real']*24)
    
    step = 1
    dim = 48
    ex =  soluce
    model=ga(function=A3_storage.objective_function_optim,dimension= dim,variable_type='int',variable_boundaries=varbound, variable_type_mixed= var_type, algorithm_parameters=param, value_step = step, exemple = ex)
    
    model.run()
    
    Instruct = list(model.output_dict['variable'])
    return Instruct
    
def improve_soluce_without_storage(soluce):
    param={'max_num_iteration': 60, 'population_size': 120, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}
    
    varbound=np.array([[60, 105]]*24)
    
    step = 1
    dim = 24
    ex =  soluce
    model=ga(function=A3_storage.objective_function_optim,dimension= dim,variable_type='int',variable_boundaries=varbound,  algorithm_parameters=param, value_step = step, exemple = ex)
    
    model.run()
    
    Instruct = list(model.output_dict['variable'])
    return Instruct
   
##Model parameters
# 3 Substations
pipes_length_3 = [randrange(90, 360, 90) for i in range(3)]
SubStations_3 = [[HEX_nom(100e3, Tr2_SS[i], functions[i], 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, pipes_length_3[i])] for i in range(3)]

SRC_3 = Source(72, 15, 450000, 3)
NET_3 = Network(SRC_3, SubStations_3)
A3 = Simulation(NET_3, ext_T)

# 10 substations
pipes_length_10 = [randrange(90, 270, 90) for i in range(10)]
SubStations10 = [[HEX_nom(100e3, Tr2_SS[i], functions[i], 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, pipes_length_10[i])] for i in range(10)]

Stor10 = Buffer(68, 50, 25, 25)
SRC_10 = Source(72, 18, 1500000, 3)
NET_10 = Network(SRC_10, SubStations10, storage_buffer = Stor10)
A10_storage = Simulation(NET_10, ext_T)

#Stor1 = Buffer(95, 50, 100, 6)
Stor3 = Buffer(68, 50, 15.6, 6)
SRC_3_stor = Source(72, 8, 450000, 3)
NET_3_storage = Network(SRC_3_stor, SubStations_3, storage_buffer = Stor3)
A3_storage=  Simulation(NET_3_storage, ext_T)

# Test 8h
pipes_length_10bis = [randrange(90, 180, 90) for i in range(10)]
SubStations10bis = [[HEX_nom(135e3, Tr2_SS1[i], functions1[i], 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, pipes_length_10[i])] for i in range(10)]
# SubStations10bis = [[HEX_nom(135e3, Tr2_SS1[i], fun1_highT[i], 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, pipes_length_10[i])] for i in range(10)]

m_NET_nom = np.sum([SS[0].m_dot2 for SS in SubStations10bis])
Stor10bis = Buffer(68, 50, 0, 100)
SRC_10bis = Source(77, m_NET_nom, 850000, 3)
NET_10bis = Network(SRC_10bis, SubStations10bis, storage_buffer = Stor10bis)
Ta_ext = [7]*4+[-3]*4 #+ [3] + [7]*3+[-2]*4 + [3] +[7]*3+[-2]*4
A10_storagebis = Simulation(NET_10bis, Ta_ext)

pipes_length_20 = [randrange(90, 180, 90) for i in range(20)]
SubStations20 = [[HEX_nom(135e3, Tr2_SS1[int(i%10)], functions1[int(i%10)], 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, pipes_length_20[i])] for i in range(20)]
# SubStations10bis = [[HEX_nom(135e3, Tr2_SS1[i], fun1_highT[i], 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, pipes_length_10[i])] for i in range(10)]

m_NET_nom = np.sum([SS[0].m_dot2 for SS in SubStations20])
Stor20bis = Buffer(68, 50, 0, 100)
SRC_20bis = Source(77, m_NET_nom, 1300000, 3)
NET_20bis = Network(SRC_20bis, SubStations20, storage_buffer = Stor20bis)
Ta_ext = [7]*4+[-3]*4 #+ [3] + [7]*3+[-2]*4 + [3] +[7]*3+[-2]*4
A20_storage = Simulation(NET_20bis, Ta_ext)

##
#Ta_ext = [7]*4+[-3]*4 + [3] + [7]*3+[-2]*4 + [3] +[7]*3+[-2]*4
#Solution with storage
I = [65.0, 70.0, 70.0, 61.0, 76.0, 77.0, 75.0, 68.0, 66.0, 70.0, 60.0, 63.0, 75.0, 75.0, 72.0, 64.0, 72.0, 70.0, 70.0, 70.0, 75.0, 75.0, 75.0, 71.0, 1.7116467980861823, -2.3951194986934126, -3.255059154370651, 0.0, 3.6222945748365554, 1.3434695204356109, -0.014941080580210223, 2.6857902793287245, -0.067666165496698, -2.5582280229456136, -3.5351552297194075, 0.0, 1.9853164179129372, 1.7977156723104089, 1.3413236802456385, 0.8690740568330402, -0.4009007904923898, -1.7977590092644522, -1.996538330823078, -1.9061913881382657, 1.3061625719203676, 1.8137306101730906, 2.8532868329899603, 0.0]

#Solution without storage
[65.0, 66.0, 69.0, 68.0, 76.0, 76.0, 75.0, 65.0, 65.0, 69.0, 66.0, 75.0, 74.0, 75.0, 74.0, 61.0, 62.0, 69.0, 69.0, 69.0, 75.0, 75.0, 74.0, 64.0]














#

I = [88.0, 76.0, 82.0, 82.0, 82.0, 81.0, 83.0, 84.0, 78.0, 83.0, 81.0, 73.0, 79.0, 75.0, 77.0, 76.0, 75.0, 77.0, 74.0, 78.0, 79.0, 78.0, 80.0, 72.0, 0.06749779648428259, 0.0, -0.1905950863498598, -1.4169685901719342, -3.1800741803646204, 0.0, -1.5888768140061742, -0.19925243628441902, 0.0, -2.531910434861875, -0.9503893374404973, -2.002626649166742, -4.891614960612433, -4.540493480710127, -1.8654397681728296, -4.90200348947989, -4.82858084283675, -3.212678787620822, 0.0, -1.1167814259368125, -2.2152456128532263, -3.891175555146509, 0.0, -2.6621606755395324]

T = I[0:24]

S  = I[24::]

T_optim=  [84.0, 75.0, 83.0, 79.0, 82.0, 80.0, 82.0, 80.0, 80.0, 79.0, 79.0, 76.0, 74.0, 74.0, 74.0, 76.0, 76.0, 76.0, 77.0, 77.0, 78.0, 79.0, 78.0, 72.0]

#with very hot storage
#Stor1 = Buffer(95, 50, 100, 6)

I = [78.0, 71.0, 78.0, 68.0, 77.0, 78.0, 60.0, 79.0, 85.0, 79.0, 79.0, 83.0, 87.0, 66.0, 74.0, 73.0, 66.0, 74.0, 75.0, 75.0, 75.0, 76.0, 76.0, 63.0, 1.2507715712251894, 0.4568522276953537, 0.5279772996585588, 1.0818689360971456, 0.7307118222293019, 0.747382551248248, 0.6930046918207831, 1.2851659380740745, 0.9, 1.1302091290576346, 0.8224781079006753, 1.0, 0.9515608068962571, 0.40673323804008044, -0.13561737907128601, 0.41354790353336845, 0.5891006773473331, 0.45822650434190315, 0.4800497233266068, -2.8857567411602227, 1.3325380028279348, 0.6210466351656141, 1.2549719638704633, 0.8457257041939981]

T = I[0:24]

S  = I[24::]

#With greater variation in the linear law
I = [74.0, 74.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 76.0, 76.0, 75.0, 73.0, 68.0, 63.0, 68.0, 68.0, 61.0, 62.0, 61.0, 71.0, 72.0, 72.0, 73.0, 73.0, 0.0, 0.0, -1.4089965925196593, 0.0, 0.0, 0.0, -4.019208268091916, 0.0, -4.642156986998234, 0.1415883746134563, 0.17376161853098182, -0.0006426813118835251, 0.3833720220876069, 0.0, 0.0, 0.0, 0.0, -2.5163805645125734, -2.6737346798390336, 0.0, -1.3698219536517406, -3.726391027341057, 0.0, 0.0]

T = I[0:24]

S  = I[24::]

#Twodays, greater linear law, weaker source


#ten days, greater linear law

def print_Ta(Ta):
    t= [i/60 for i in range(60*len(Ta))]
    T =[]
    for Ti in Ta:
        T += [Ti]*60
        
    plt.figure()
    plt.plot(t, T, label = 'outside Temperature')
    plt.legend()
    plt.show()