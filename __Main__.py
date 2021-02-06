import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
from random import gauss, randrange
from itertools import *
import time

#from GA_algorithm import geneticalgorithm as ga
from geneticalgorithm import geneticalgorithm as ga

from Components.Data_meteo_reading import Ts2_A, Ts2_B, Ts2_C, Ts2_D
from Components.Data_meteo_reading import Ta_oneday as ext_T
from Components.Data_meteo_reading import Ta_week, functions, Tr2_SS
from Components.Networks import Network
from Components.Source import Source
from Components.Heat_exchanger import HEX_nom #, Source
from Components.Pipes import Pipe
from Components.Ressources import *
from Components.Storage import Buffer
from Model_SIM import Simulation

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
        
    #print(len(MOD.PBoiler_Ta), MOD.PBoiler_Ta)
    
    if a > 1.2:
        print(a)
        param = {'max_num_iteration': 60, 'population_size': 110, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}
        if a > 2:
            param = {'max_num_iteration': 60, 'population_size': 130, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}
    
    if Storage_dim is None:
        varbound=np.array([[60, 105]]*dim)
        var_type = None
        if dim == len(MOD.Ta):
            ex = [math.ceil(MOD.f_Ts1(T)) for T in MOD.Ta]
        else:
            ex = None
            
    else:
        varbound=np.array([[60, 105]]*dim + [[-5, 5]]*Storage_dim)
        var_type = np.array(['int'] * dim + ['real']*Storage_dim)
        if dim == len(MOD.Ta):
            ex = [math.ceil(MOD.f_Ts1(T)) for T in MOD.Ta] + [0]*Storage_dim
        else:
            ex = None
            
        dim += Storage_dim
        
        
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
        self.Tboiler_Ta = None
        self.Ta_week = None 
        self.current_state = None
        
        self.current_Ta = None
    
def optim_week(dim, MOD, Result_class, param = algorithm_param, Ta_w = Ta_week, step_optim_h = 24):
    
    Result_class.Ta_week = np.copy(Ta_week)
    days = [Ta_w[i:i + 24] for i in range(0, 168, step_optim_h)]   #24*7 = 168
    Result_class.T_Boiler_optim = []
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
        
        # try:
        #     MOD.initialisation(Instruct_T_Boiler[0:step_optim_h])
        #     Result_class.T_Boiler_optim.extend(Instruct_T_Boiler[0:step_optim_h])
        #     
        # except ValueError: #we try to refined the solution found 
        #     print('Attempt to refined the previous solution')
        #     MOD.save_init_ = Result_class.current_state.copy()
        #     MOD.Ta = day
        #     MOD.initialised = False
        #     MOD.initialisation()
        #     
        #     varbound=np.array([boundaries(x, MOD.nb_SS) for x in Instruct_P_Boiler] )

        #      model=ga(function=MOD.objective_function,dimension= dim,variable_type='int',variable_boundaries=varbound, algorithm_parameters=param, value_step = 100, exemple = Instruct_P_Boiler)
        #     
        #     model.run()
        #     
        #     Instruct_P_Boiler = list(model.output_dict['variable'])
        #     print('refined solution')
        #     
        #     
        #     MOD.initialisation()  #Setting the parameters to init values
        #     MOD.initialised = False
        #     MOD.Ta = day[0:step_optim_h]
        #     
        #     MOD.initialisation(Instruct_P_Boiler[0:step_optim_h])
        #     Result_class.P_Boiler_optim.extend(Instruct_P_Boiler[0:step_optim_h])
            
            
    
    time2 = time.time()
    print(f'time = {time2 - time1}')
    # MOD.save_init_= state_init
    # MOD.initialisation()
    # MOD.Ta = Ta_w
    #MOD.plot(Pboiler)
    
    #return Pboiler, Pboiler_Ta

Result_3 = Result()
    
##Model parameters
# SS_A = [HEX_nom(100e3, 42, Ts2_A, 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000)]
# SS_B = [HEX_nom(100e3, 39, Ts2_B, 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350)]
# SS_C = [HEX_nom(100e3, 36, Ts2_C, 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500)]
# SS_D = [HEX_nom(100e3, 43, Ts2_D, 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 700)]
# SRC_3 = Source(70, 20, 450000, 3)
# NET = Network_bOptim(SRC1, 0, [SS_A])
# NET3 = Network_bOptim(SRC1, 0, [SS_A, SS_B, SS_C])
# A = Simulation_Ta(NET, ext_T)
# A3 = Simulation_Ta(NET3, ext_T)

#3 Substations
pipes_length_3 = [randrange(350, 800, 50) for i in range(3)]
SubStations_3 = [[HEX_nom(100e3, Tr2_SS[i], functions[i], 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, pipes_length_3[i])] for i in range(3)]

SRC_3 = Source(72, 15, 450000, 3)
NET_3 = Network(SRC_3, SubStations_3)
A3 = Simulation(NET_3, ext_T)

# 10 substations
pipes_length_10 = [randrange(350, 800, 50) for i in range(10)]
SubStations10 = [[HEX_nom(100e3, Tr2_SS[i], functions[i], 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, pipes_length_10[i])] for i in range(10)]

SRC_10 = Source(72, 40, 1500000, 3)
NET_10 = Network(SRC_10, SubStations10)
A10 = Simulation(NET_10, ext_T)

Stor1 = Buffer(70, 50, 15.6, 6)
NET_3_storage = Network(SRC_3, SubStations_3, storage_buffer = Stor1)
A3_storage=  Simulation(NET_3_storage, ext_T)

def improve_soluce_storage(soluce):
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
    
    
##
I = [88.0, 76.0, 82.0, 82.0, 82.0, 81.0, 83.0, 84.0, 78.0, 83.0, 81.0, 73.0, 79.0, 75.0, 77.0, 76.0, 75.0, 77.0, 74.0, 78.0, 79.0, 78.0, 80.0, 72.0, 0.06749779648428259, 0.0, -0.1905950863498598, -1.4169685901719342, -3.1800741803646204, 0.0, -1.5888768140061742, -0.19925243628441902, 0.0, -2.531910434861875, -0.9503893374404973, -2.002626649166742, -4.891614960612433, -4.540493480710127, -1.8654397681728296, -4.90200348947989, -4.82858084283675, -3.212678787620822, 0.0, -1.1167814259368125, -2.2152456128532263, -3.891175555146509, 0.0, -2.6621606755395324]

T = I[0:24]

S  = I[24::]

T_optim=  [84.0, 75.0, 83.0, 79.0, 82.0, 80.0, 82.0, 80.0, 80.0, 79.0, 79.0, 76.0, 74.0, 74.0, 74.0, 76.0, 76.0, 76.0, 77.0, 77.0, 78.0, 79.0, 78.0, 72.0]