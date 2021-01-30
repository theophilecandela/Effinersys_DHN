import numpy as np
from numpy import log as ln
from matplotlib import pyplot as plt
import math
from random import gauss
from itertools import *
import time

from geneticalgorithm import geneticalgorithm as ga
from Components.Data_meteo_reading import Ts2_A, Ts2_B, Ts2_C, Ts2_D
from Components.Data_meteo_reading import Ta_oneday as ext_T
from Components.Data_meteo_reading import Ta_week
from Components.Networks import Network_bOptim
from Components.Source import Source
from Components.Heat_exchanger import HEX_nom
from Components.Pipes import Pipe
from Components.Ressources import *
from Model_SIM import Simulation_Ta

##OPTIM Function
algorithm_param={'max_num_iteration': 100, 'population_size': 100, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}

def boundaries(x, nb_SS):
    sup = math.ceil(x)
    inf = math.floor(x)
    if sup < 90000 * nb_SS:
        return [inf//2, 2*sup]
    else:
        return [inf//4, sup]

def optim(dim, MOD, param = algorithm_param, step = 100, plot = False):
    '''Uses Genetic algorithm to calculate an optimal/sub-optimal list of instructions;
    dim: number of instruction throughout the considered period of time. Cannot be higher than the number of time iterations during a simulation (ie total_t/dt)
    MOD: an object of class Simulation_Ta
    param: a dict containing the parameters of the optimization process
    '''
    time1 = time.time()
    MOD.refined_Ts1(dim)
    try:
        a = MOD.objective_function_Ta(dim_perday = dim)
    except ValueError:
        a = 2
    n = len(MOD.PBoiler_Ta)
    step = n//dim
            
    varbound=np.array([boundaries(x, MOD.nb_SS) for x in MOD.PBoiler_Ta[step-1::step]] )

    model=ga(function=MOD.objective_function,dimension= dim,variable_type='int',variable_boundaries=varbound, algorithm_parameters=param, value_step = step, exemple = MOD.PBoiler_Ta[step-1::step] )
    
    model.run()
    time2 = time.time()
    print(f'[Optimization process] - Time: {time2-time1}')
    
    Boiler_instructP = list(model.output_dict['variable'])
    if plot:
        MOD.plot(Boiler_instructP)
    return Boiler_instructP
    
    
def optim_stepped(dim, MOD, param = algorithm_param):
    
    time1 = time.time()
    
    MOD.refined_Ts1(dim)
    try:
        a = MOD.objective_function_Ta(dim_perday = dim)
    except ValueError:
        a = 2
    n = len(MOD.PBoiler_Ta)
    step = n//dim
    
    def f(x):
        sup = math.ceil(x)
        inf = math.floor(x)
        if sup < 90000 * MOD.nb_SS:
            return [inf//2, 2*sup]
        else:
            return [inf//4, sup]

    varbound=np.array([f(x) for x in MOD.PBoiler_Ta[step-1::step]] )
    
    time2 = time.time()
    print(f'First step [Refining linear T] - Time: {time2- time1} \n')
    
    steps = [10000, 1000, 100]
    parameters = [(70, 60), (50,40), (50,40)]
    
    time1 = time.time()
    for i, step in enumerate(steps):
        a, b = parameters[i]
        param={'max_num_iteration': a, 'population_size': b, 'mutation_probability': 0.1, 'elit_ratio': 0.2, 'crossover_probability': 0.65, 'parents_portion': 0.3, 'crossover_type': 'uniform', 'max_iteration_without_improv': None}
        
        model=ga(function=MOD.objective_function,dimension= dim,variable_type='int',variable_boundaries=varbound, algorithm_parameters=param, value_step = step)
    
        t1 = time.time()
        model.run()
        t2 = time.time()
        print(f'Step = {step} [Genetic algorithm optimization process] - Time: {t2-t1}')
        
        Boiler_instructP = list(model.output_dict['variable'])
        varbound = np.array([[max(0, x - step), x + step] for x in Boiler_instructP] )
        
    time2 = time.time()
    print(f'Second step [GA-Total] - Time: {time2- time1} \n')   
    
    MOD.plot(Boiler_instructP)
    
    return Boiler_instructP
    
    
    
def optim_week(dim, MOD, param = algorithm_param, Ta_w = Ta_week, step_optim_h = 24):
    # T_begin = Ta_w[0]
    # init_day = list(np.linspace(10, T_begin, 24))
    # MOD.Ta = init_day
    # MOD.refined_Ts1()
    # Instruct_P_Boiler = MOD.PBoiler_Ta
    # MOD.initialisation()
    # MOD.initialised = False
    # MOD.initialisation(Instruct_P_Boiler)
    
    state_init = MOD.save_init_
    days = [Ta_w[i:i + 24] for i in range(0, len(Ta_w), step_optim_h)]
    Pboiler = []
    Pboiler_Ta = []
    MOD.initialisation()
    time1 = time.time()
    for day in days:
        print('___________ \n')
        print(Pboiler)
        
        MOD.Ta = day
        Instruct_P_Boiler = optim(dim, MOD,step = 1000)
        Pboiler_Ta.extend(MOD.PBoiler_Ta)
        
        MOD.initialisation()
        MOD.initialised = False
        MOD.Ta = day[0:step_optim_h]
        MOD.initialisation(Instruct_P_Boiler[0:step_optim_h])
        Pboiler.extend(Instruct_P_Boiler[0:step_optim_h])
    
    time2 = time.time()
    print(f'time = {time2 - time1}')
    MOD.save_init_= state_init
    MOD.initialisation()
    MOD.Ta = Ta_w
    MOD.plot(Pboiler)
    
    return Pboiler
        
    
##Model parameters
SS_A = [HEX_nom(100e3, 42, Ts2_A, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 1000)]
SS_B = [HEX_nom(100e3, 39, Ts2_B, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 350)]
SS_C = [HEX_nom(100e3, 36, Ts2_C, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 500)]
SS_D = [HEX_nom(100e3, 43, Ts2_D, 5000, 2, 2.1), Pipe(lam_i, lam_p, lam_s, R_int, R_p, R_i, z, 700)]


SRC1 = Source(70, 20, 450000, 3)

NET = Network_bOptim(SRC1, 0, [SS_A])
NET3 = Network_bOptim(SRC1, 0, [SS_A, SS_B, SS_C])

A = Simulation_Ta(NET, ext_T)
A3 = Simulation_Ta(NET3, ext_T)


SRC2 = Source(75, 20, 600000, 3)
NET4 = Network_bOptim(SRC2, 0, [SS_A, SS_B, SS_C, SS_D])
A4 = Simulation_Ta(NET4, ext_T)


