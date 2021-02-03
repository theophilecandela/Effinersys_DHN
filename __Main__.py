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
from Components.Heat_exchanger import HEX_nom
from Components.Pipes import Pipe
from Components.Ressources import *
from Model_SIM import Simulation

##OPTIM Function
algorithm_param={'max_num_iteration': 60, 'population_size': 90, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}

def optim(dim, MOD, param = algorithm_param, step = 1, plot = False):
    '''Uses Genetic algorithm to calculate an optimal/sub-optimal list of instructions;
    dim: number of instruction throughout the considered period of time. Cannot be higher than the number of time iterations during a simulation (ie total_t/dt)
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
            
    varbound=np.array([[60, 105]]*dim)
    if dim == len(MOD.Ta):
        ex = [MOD.f_Ts1(T) for T in MOD.Ta]
    else:
        ex = None
        
    model=ga(function=MOD.objective_function,dimension= dim,variable_type='int',variable_boundaries=varbound, algorithm_parameters=param, value_step = step, exemple = ex)
    
    model.run()
    time2 = time.time()
    print(f'[Optimization process] - Time: {time2-time1}')
    
    Boiler_instructT = list(model.output_dict['variable'])
    if plot:
        MOD.plot(Boiler_instructT)
    return Boiler_instructT
       
    
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

SRC_3 = Source(72, 40, 450000, 3)
NET_3 = Network(SRC_3, SubStations_3)
A3 = Simulation(NET_3, ext_T)

# 10 substations
pipes_length_10 = [randrange(350, 800, 50) for i in range(10)]
SubStations10 = [[HEX_nom(100e3, Tr2_SS[i], functions[i], 5000, 2.8, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, pipes_length_10[i])] for i in range(10)]

SRC_10 = Source(72, 40, 1500000, 3)
NET_10 = Network(SRC_10, SubStations10)
A10 = Simulation(NET_10, ext_T)

