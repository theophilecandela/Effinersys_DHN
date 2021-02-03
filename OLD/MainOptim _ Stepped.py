    
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