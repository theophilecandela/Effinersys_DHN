import numpy as np
from numpy import log as ln

## Datas
#Pipe parameters
#ambient temperature
Ta = 283.15
#Thermal conductivity of pipe isolation (W/m/K)
lam_i = 0.025
#Thermal conductivity of steel (W/m/K)
lam_p = 50.2
#Thermal conductivity of the ground (concrete) (W/m/K) 
lam_s = 1.6
#Thermal conductivity of water (W/m/K)
lam_e = 0.6071
#Thermal capacity of water (J/K/kg)
Cp = 4185
#Dynamic viscosity of water (Pa.s)
mu_e = 0.535e-3
#volumetric mass density of water (kg/m^3)
rho = 1000 
#Depth of buried pipe (m)
z = 2
#Pipe radius (m)
R_int = 200e-3
R_p = 250e-3
R_i = 400e-3


#cost-kilowatt-hour
efficiency_gaz = 0.6
CarbonperkWh = 490 #g
carbonTtaxe = 50 #€
C_kWh = 1/efficiency_gaz * CarbonperkWh*(10**(-6)) * 50



#Numerical parameters
dx = 90  #spatial discretization step (m)
dt = 60  #discretization time step (s)

##Functions
def eff(N, R):
    return (1 - np.exp((R-1)*N))/(1-R*np.exp((R-1)*N))
    
def newton(f, f_prime, x0, eps = 0.0001):
    x = x0
    n = 15
    k = 0
    while np.abs(f(x)) > eps and k < n:
        x = max(0.001, x - f(x)/f_prime(x))
        k += 1
    if k>= n:
        raise ValueError('Ne converge pas')
    else:
        return x
        
