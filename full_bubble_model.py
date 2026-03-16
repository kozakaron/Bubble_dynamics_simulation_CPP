"""
This program contains the JIT compiled (high performance) ODE function (_f()), the numerical solvers (solve()), and several supplementary functions. 
Before use, make sure to set the variables below "___Settings___" to your needs. 

Usage:
 * Use inp_data_extractor.py to turn a .inp file into parameters.py
 * Set the parameters in the "___Settings___" section
 * Import this file: from Bubble_dynamics_simulation import full_bubble_model as de
 * Assemble the control parameter dictionary: cpar = de.example_cpar() | 
   You may print cpar: de.print_cpar(cpar) | 
   The result can be copied, modified and used as code. 
 * You may solve the differential equation: num_sol, error_code, elapsed_time = de.solve(cpar)
 * Retrieve post processing data: data = de.get_data(cpar, num_sol, error_code, elapsed_time)
 * Use de.plot(cpar) to plot the results. This function has several costumization options and extra plots.

 Notes:
  * See the documentation or run help(de.solve) for more information. Substitute any function name instead of de.solve.
  * Use the Make_dir class to save simulation results.
  * See the example files. Example files start with a capital letter, and have the format .ipynb
  * You may define any excitation function in excitation.py with a fix number of excitation control parameters.
  * You can plot any extra variables with plot() by setting plot_extra=True. 
    To do this, modify the bottom lines in the _f() function, and the first lines of plot().
"""

"""________________________________Settings________________________________"""

enable_heat_transfer = True
enable_evaporation = False
enable_reactions = True
enable_dissipated_energy = True
enable_reaction_rate_threshold = True
enable_time_evaluation_limit = False
target_specie = 'NH3' # Specie to calculate energy demand for
excitation_type = 'no_excitation' # function to calculate pressure excitation (see excitation.py for options)

"""________________________________Libraries________________________________"""

from termcolor import colored   # colored error messages
import matplotlib.pyplot as plt   # for plotting
import numpy as np   # matrices, math
from scipy.integrate import solve_ivp   # differential equation solver
from scipy.signal import argrelmin   # loc min finding
import time   # runtime measurement
from datetime import datetime   # for accessing current datetime
import socket   # for accessing computer name
import psutil   # get system information
from numba import njit   # Just In Time compiler
from numba.types import Tuple, unicode_type, float64, float32, int64, int32   # JIT types
from func_timeout import func_timeout, FunctionTimedOut   # for timeout
import os    # file management
import importlib   # for reloading your own files
import traceback   # for error handling

# import parameters.py as par:
try:
    import parameters as par   # numeric constants and coefficents
    importlib.reload(par)   # reload changes you made
except Exception as _error:
    print(print(colored('Error, \'parameters.py\' not found','red')))
    raise _error

# import excitation.py as excitation:
try:
    import excitation
    importlib.reload(excitation)
except ImportError as _error:
    try:
        from  Bubble_dynamics_simulation import excitation
        importlib.reload(excitation)
    except ImportError as _error:
        print(colored(f'Error, \'excitation.py\' not found', 'red'))
        raise _error
    except Exception as _error:
        print(colored(f'Error, \'excitation.py\' failed to load', 'red'))
        raise _error
except Exception as _error:
    print(colored(f'Error, \'excitation.py\' failed to load', 'red'))
    raise _error


"""________________________________General________________________________"""
previous_time = None
        
class dotdict(dict):
    """Dot notation access to dictionary attributes. 
    Instead of dictionary['key'] you can use dictionary.key"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def copy(input):
    """Deep copies the input. Use: copied = copy(original). 
    Deep copy means that the input and the output are independent from each other, nothing is copied by reference."""

    if type(input) == list:
        return [copy(element) for element in input]
    elif type(input) == dict:
        return {key: copy(value) for key, value in input.items()}
    elif type(input) == np.ndarray:
        return input.copy()
    elif type(input) == dotdict:
        return dotdict({key: copy(value) for key, value in input.items()})
    else:
        return input


"""________________________________Load excitation________________________________"""

def _colorTF(boolean):
    """Input: boolean. Colors it red if False, green if True. Returns the colored string."""
    return colored(str(boolean), 'green') if boolean else colored(str(boolean), 'red')

Excitation, excitation_args, excitation_units, excitation_defaults = excitation.getExcitation(excitation_type=excitation_type)
if par.indexOfWater == -1:
    enable_evaporation = False
print(f'model: {par.model}')
print(f'target specie: {target_specie}')
print(f'excitation: {excitation_type} (control parameters: {excitation_args})')
print(f'enable heat transfer: {_colorTF(enable_heat_transfer)}\tenable evaporation: {_colorTF(enable_evaporation)}\tenable reactions: {_colorTF(enable_reactions)}\tenable dissipated energy: {_colorTF(enable_dissipated_energy)}\tenable reaction rate threshold: {_colorTF(enable_reaction_rate_threshold)}\tenable_time_evaluation_limit: {_colorTF(enable_time_evaluation_limit)}')
if target_specie not in par.species:
    print(colored(f'Error, target specie \'{target_specie}\' not found in parameters.py', 'red'))


"""________________________________Control parameters________________________________"""

cpar_rules = dict(
    ID =         dict(default=0,                      type=int,           range=[0, 10**15],              comment='ID of control parameter (not used during calculation)'),
# Initial conditions:
    R_E =        dict(default=10.0e-6,                type=float,         range=[1e-9, 1.0e3],             comment='bubble equilibrium radius [m]'),
    ratio =      dict(default=1.0,                    type=float,         range=[1.0, 1000.0],            comment='initial radius / equilibrium radius R_0/R_E [-]'),
    gases =      dict(default=[0],                    type=(list, int),   range=[0, par.K-1],             comment='indexes of species in initial bubble (list of species indexes)'),
    fractions =  dict(default=[1.0],                  type=(list, float), range=[0.0, 1.0],               comment='molar fractions of species in initial bubble (list of fractions for every gas)'),
# Ambient parameters:
    P_amb =      dict(default=1.0*par.atm2Pa,         type=float,         range=[0.0, 1.0e10*par.atm2Pa], comment='ambient pressure [Pa]'),
    T_inf =      dict(default=20.0+par.absolute_zero, type=float,         range=[0.0, 100000.0],          comment='ambient temperature [K]'),
# Liquid parameters:
    alfa_M =     dict(default=par.alfa_M,             type=float,         range=[0.0, 10.0],              comment='water accommodation coefficient [-]'),
    P_v =        dict(default=par.P_v,                type=float,         range=[0.0, 1e8],               comment='vapour pressure [Pa]'),
    mu_L =       dict(default=par.mu_L,               type=float,         range=[0.0, 10000.0],           comment='dynamic viscosity [Pa*s]'),
    rho_L_ref =      dict(default=998.2,              type=float,         range=[0.0, 100000.0],          comment='density [kg/m^3]'),
    c_L_ref =        dict(default=1483.0,                type=float,         range=[0.0, 100000.0],          comment='sound speed [m/s]'),
    surfactant = dict(default=1.0,                    type=float,         range=[0.0, 1000.0],            comment='surfactant (surface tension modfier) [-]'),
)
# Excitation parameters:
for _arg, _unit, _default in zip(excitation_args, excitation_units, excitation_defaults):
    cpar_rules[_arg] = dict(default=_default,           type=float,         range=[-1e30, 1e30],            comment=f'[{_unit}]')

def check_cpar(cpar):
    """Checks the existence, type, and value of each required cpar keys according to cpar_rules. Prints colored error messages.
    Returns True if cpar is correct, False otherwise. Input: cpar - dict or dotdict of control parameters."""

    for key, rule in cpar_rules.items():
        if key not in cpar: # check existence
            if key == 'ID':
                cpar[key] = rule['default']
                continue
            print(colored(f'Error in cpar, key \'{key}\' not found. ', 'red'))
            return False

        if type(rule['type']) == type:
            if type(cpar[key]) != rule['type']: # check type
                try:
                    cpar[key] = rule['type'](cpar[key])
                except:
                    print(colored(f'Error in cpar, key \'{key}\' has wrong type. Expected {rule["type"].__name__}, recieved {type(cpar[key]).__name__}: cpar.{key} = {cpar[key]}', 'red'))
                    return False
            if cpar[key] < rule['range'][0] or rule['range'][1] < cpar[key]:    # check range
                print(colored(f'Error in cpar, key \'{key}\' out of range: cpar.{key} = {cpar[key]} not in {rule["range"]}', 'red'))
                return False
        else:
            list_type, element_type = rule['type']
            if type(cpar[key]) != list_type:    # check list type
                try:
                    cpar[key] = list_type(cpar[key])
                except:
                    print(colored(f'Error in cpar, key \'{key}\' has wrong type. Expected {list_type.__name__}, recieved {type(cpar[key]).__name__}: cpar.{key} = {cpar[key]}', 'red'))
                    return False
            if len(cpar[key]) == 0:   # check list length
                print(colored(f'Error in cpar, key \'{key}\' is empty: cpar.{key} = {cpar[key]}', 'red'))
                return False
            for i, element in enumerate(cpar[key]):
                if type(element) != element_type:   # check element type
                    try:
                        cpar[key][i] = element_type(element)
                    except:
                        print(colored(f'Error in cpar, key \'{key}\' has wrong type. Expected list of {element_type.__name__}, recieved list of {type(element).__name__}: cpar.{key} = {cpar[key]}', 'red'))
                        return False
                if cpar[key][i] < rule['range'][0] or rule['range'][1] < cpar[key][i]:    # check element range
                    print(colored(f'Error in cpar, key \'{key}\' out of range: cpar.{key} = {cpar[key]} not in {rule["range"]}', 'red'))
                    return False
    
    if round(sum(cpar['fractions']), 5) != 1.0:
        print(print(colored(f'Error in cpar, sum of cpar.fractions is not 1: {cpar["fractions"]}', 'red')))
        return False
    if len(cpar['gases']) != len(cpar['fractions']):
        print(print(colored(f'Error in cpar, len(cpar.gases) != len(cpar.fractions): {cpar["gases"]} != {cpar["fractions"]}', 'red')))
        return False
    return True
    
def example_cpar(normal_dict=False):
    """Provides an example of the control parameter dictionary. Use print_cpar() to print it. Parameters:
    * normal_dict: if True, returns a normal dictionary, else returns a dotdict
    
    Returns:
    * cpar: control parameter dictionary"""
    
    cpar = {key: rule['default'] for key, rule in cpar_rules.items()}

    if target_specie == 'NH3':
        cpar['gases'] = [par.index['H2'], par.index['N2']]
        cpar['fractions'] = [0.75, 0.25]
    else:
        cpar['gases'] = [par.index['O2']]
        cpar['fractions'] = [1.0]
    if excitation_type == 'no_excitation':
        cpar['ratio'] = 3.0

    if normal_dict:
        return cpar
    else:
        return dotdict(cpar)

@njit(float64(float64))
def vapour_pressure(T): # [K]
    """Calculates the vapour pressure of water as a function of temperature (T) in Kelvin. """
    T -= 273.15 # [°C]
    return 611.21 * np.exp( (18.678 - T / 234.5) * (T / (257.14 + T)) ) # [Pa]

@njit(float64(float64))
def viscosity(T):   # [K]
    """Calculates the dynamic viscosity of water as a function of temperature (T) in Kelvin. Pressure dependence is neglected. """
    return 1.856e-14 * np.exp(4209.0/T + 0.04527*T - 3.376e-5*T**2) # [Pa*s]

# To retain compatibility with the old code
def VapourPressure(T):
    """Legacy function, use vapour_pressure(T) instead"""
    return vapour_pressure(T)
def Viscosity(T):
    """Legacy function, use viscosity(T) instead"""
    return viscosity(T)


"""________________________________Preprocessing________________________________"""
def _initial_condition(cpar, evaporation=False, extra_dims=0):
    """Calculates the initial condition of the bubble from the control parameters. Arguments:
     * cpar: control parameters, dotdict
     * evaporation: if True, water vapour will be present in the initial bubble (the water's partial pressure equals to the saturated water pressure)
     * extra_dims: add extra dimensions to the initial condition array (initial value: 0.0) | 
                   use it to plot extra variables (e.g. energy) during the simulation

    Returns:
     * IC: initial condition array
     * lowpressure_error: if True, the pressure of the gas is negative
     * lowpressure_warning: if True, the pressure during the expansion is lower, than the saturated water pressure
    """    
    IC = np.zeros((par.K + 4 + extra_dims), dtype=np.float64)
    R_0 = cpar.ratio * cpar.R_E

    # Equilibrium state
    p_E = cpar.P_amb + 2.0 * cpar.surfactant * par.sigma / cpar.R_E # [Pa]
    p_gas = p_E - cpar.P_v if evaporation else p_E
    
    c_V = cpar.R_gas /(cpar.kappa-1.0)
    b=4.467e-3
    rho_gas = p_gas / ((cpar.kappa-1.0)*c_V*cpar.T_inf + b*p_gas)
    V_norhc = 4.0 / 3.0 * cpar.R_E**3.0 * np.pi # [m^3]
    cpar.r_hc=(b*3.0*rho_gas*V_norhc/(4.0*np.pi))**(1.0/3.0)
    
    V_E = 4.0 / 3.0 * cpar.R_E**3.0 * np.pi#(cpar.R_E**3.0-cpar.r_hc**3.0) * np.pi # [m^3]
    lowpressure_error = lowpressure_warning = False
    if p_gas < 0.0:
        #print(colored('Error in _initial_condition(), the pressure of the gas is negative.', 'red'))
        lowpressure_error = True
    n_gas = p_gas * V_E / (par.R_g * cpar.T_inf) # [mol]
    
    
    # Isotermic expansion
    V_0 = 4.0 / 3.0 * R_0**3.0 * np.pi#(R_0**3-cpar.r_hc**3) * np.pi    # [m^3]
    n_H2O = cpar.P_v * V_0 / (par.R_g * cpar.T_inf) if evaporation else 0.0 # [mol]
    c_H2O = n_H2O / V_0    # [mol/m^3]
    c_gas = n_gas / V_0    # [mol/m^3]
    p_gas = c_gas * par.R_g * cpar.T_inf # [Pa]
    P_amb_min = cpar.P_v if evaporation else 0.0 # [Pa]
    P_amb_min += p_gas - 2.0 * cpar.surfactant * par.sigma / R_0 # [Pa]
    if P_amb_min < cpar.P_v:
        #print(colored('Warning in _initial_condition(), the pressure during the expansion is lower, than the saturated water pressure.', 'yellow'))
        lowpressure_warning = True

    # Initial conditions
    IC[0] = R_0   # R_0 [m]
    IC[1] = 0.0    # dRdt_0 [m/s]
    IC[2] = cpar.T_inf   # T_0 [K]
    if evaporation and cpar.indexOfWater != -1:
        IC[3 + par.index['H2O']] = c_H2O * 1.0e-6    # [mol/cm^3]
    for index, fraction in zip(cpar.gases, cpar.fractions):
        IC[3 + index] = fraction * c_gas * 1.0e-6    # [mol/cm^3]
    IC[3 + par.K] = 0.0 # dissipated energy [J]

    return IC, lowpressure_error, lowpressure_warning


def _work(cpar, evaporation=False):
    """Calculates expansion work of the bubble, if ratio != 1"""

    R_0 = cpar.ratio * cpar.R_E # [m]
    V_E = 4.0 / 3.0 * cpar.R_E**3.0 * np.pi#(cpar.R_E**3-cpar.r_hc**3) * np.pi    # [m^3]
    V_0 = 4.0 / 3.0 * R_0**3.0 * np.pi#(R_0**3-cpar.r_hc**3) * np.pi  # [m^3]
    
    p_E = cpar.P_amb + 2 * cpar.surfactant * par.sigma / cpar.R_E # [Pa]
    p_gas = p_E - cpar.P_v if evaporation else p_E # [Pa]
    n_gas = p_gas * V_E / (par.R_g * cpar.T_inf) # [mol]
    
    W_gas0 = -(cpar.P_v * V_0 + n_gas * par.R_g * cpar.T_inf * np.log(V_0)) if evaporation else -(n_gas * par.R_g * cpar.T_inf * np.log(V_0))
    W_gasE = -(cpar.P_v * V_E + n_gas * par.R_g * cpar.T_inf * np.log(V_E)) if evaporation else -(n_gas * par.R_g * cpar.T_inf * np.log(V_E))
    W_gas = W_gas0 - W_gasE # [J]
    
    W_surface_tension = par.sigma * cpar.surfactant *  4.0 * np.pi * (R_0**2 - cpar.R_E**2) # [J]
    W_flow = cpar.P_amb * 4.0 / 3.0 * np.pi * (R_0**3.0 - cpar.R_E**3.0) # [J]
    return W_gas + W_surface_tension + W_flow    # [J]


"""________________________________pressures________________________________"""

@njit(Tuple((float64, float64))(float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64[:]))
def _pressure(t, R, R_dot, mu_L, surfactant, rho_L_ref, p_L_ref, B_L, b_L, Gamma_L, p, p_dot, P_amb, args):
    (p_Inf, p_Inf_dot) = Excitation(t, P_amb, args)
    p_L = p - (2.0 * surfactant * par.sigma + 4.0 * mu_L * R_dot) / R
    p_L_dot_e = p_dot + (2.0 * surfactant * par.sigma * R_dot + 4.0 * mu_L * R_dot ** 2) / (R ** 2)
    K_L=rho_L_ref/((p_L_ref+B_L)**(1.0/Gamma_L)*(1.0-b_L*rho_L_ref))
    rho_L=K_L*(p_L+B_L)**(1.0/Gamma_L)/(1.0+b_L*K_L*(p_L+B_L)**(1.0/Gamma_L))
    rho_Inf=K_L*(p_Inf+B_L)**(1.0/Gamma_L)/(1.0+b_L*K_L*(p_Inf+B_L)**(1.0/Gamma_L))
    H=(Gamma_L/(Gamma_L-1.0)*(p_L+B_L)/rho_L-Gamma_L*b_L/(Gamma_L-1.0)*(p_L+B_L)+b_L*p_L)-(Gamma_L/(Gamma_L-1.0)*(p_Inf+B_L)/rho_Inf-Gamma_L*b_L/(Gamma_L-1.0)*(p_Inf+B_L)+b_L*p_Inf) #h_L-h_Inf
    H_dot_e=p_L_dot_e/rho_L-p_Inf_dot/rho_Inf
    c_L=(Gamma_L*(p_L+B_L)/(rho_L-b_L*rho_L*rho_L))**0.5
    Nom=((1.0+R_dot/c_L)*H-3.0/2.0*(1.0-R_dot/(3.0*c_L))*R_dot*R_dot)/((1.0-R_dot/c_L)*R)+H_dot_e/c_L
    Den=1.0+4.0*mu_L/(rho_L*R*c_L)
    return Nom, Den


"""________________________________NASA polynomials________________________________"""

# returns molar heat capacities, enthalpies and entropies
@njit(float64[:, :](float64))
def _thermodynamic(T):
    ret = np.zeros((4, par.K), dtype=np.float64)   # [C_p, H, S, C_v]
    for k in range(par.K):
    # get coefficients for T
        if T <= par.TempRange[k][2]: # T <= T_mid
            a = par.a_low[k]
        else:  # T_mid < T
            a = par.a_high[k]
    # calculate sums
        C_p = H = S = 0.0
        for n in range(par.N): # [0, 1, 2, 3, 4]
            T_pow = T**n
            C_p += a[n] * T_pow
            H += a[n] * T_pow / float64(n+1)
            if n != 0:
                S += a[n] * T_pow / float64(n)
    # calculations outside the sums
        # Molar heat capacities at constant pressure (isobaric) [erg/mol/K]
        ret[0][k] = par.R_erg * C_p
        # Enthalpies [erg/mol]
        ret[1][k] = par.R_erg * (T * H + a[par.N])
        # Entropies [erg/mol/K]
        ret[2][k] = par.R_erg * (a[0] * np.log(T) + S + a[par.N+1])
        # Molar heat capacities at constant volume (isochoric) [erg/mol/K]
        ret[3][k] = ret[0][k] - par.R_erg

    return ret


"""________________________________Evaporation________________________________"""

@njit(Tuple((float64, float64))(float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64))
def _evaporation(p, T, X_H2O, H_steam, alfa_M, sigma_evap, Gamma, T_inf, P_v, C_4_starred, C_p_water, J2erg):
# condensation and evaporation
    p_H2O = X_H2O * p
    #Old:
    
    n_eva_dot = 1.0e3 * alfa_M * P_v / (par.W[par.indexOfWater] * np.sqrt(2.0 * np.pi * par.R_v * T_inf))
    n_con_dot = 1.0e3 * Gamma * alfa_M * p_H2O / (par.W[par.indexOfWater] * np.sqrt(2.0 * np.pi * par.R_v * T))
    n_net_dot = sigma_evap * (n_eva_dot - n_con_dot)
    m_net_dot = n_net_dot * par.W[par.indexOfWater] / 1.0e3 #mol/s
    
    #m_net_dot =  p_H2O/(P_v  * (2.0*np.sqrt(np.pi) * (1.0-alfa_M)/alfa_M - C_4_starred)) * np.sqrt(2.0/par.R_v) * (P_v-p_H2O)/np.sqrt(T_inf) #kg/s
    #n_net_dot = m_net_dot * par.W[par.indexOfWater] / 1.0e3 #mol/s

    #Old:
# Molar heat capacity of water at constant volume (isochoric) [J/mol/K]
    # get coefficients for T
    #if T <= par.TempRange[par.indexOfWater][2]: # T <= T_mid
    #    a = par.a_low[par.indexOfWater]
    #else:  # T_mid < T
    #    a = par.a_high[par.indexOfWater]
    # calculate sum

    #C_V = 0.0
    #for n in range(par.N): # [0, 1, 2, 3, 4]
    #    C_V += a[n] * T**n
    #C_V = par.R_erg * (C_V - 1.0)
    
    # get coefficients for T
    #if T_inf <= par.TempRange[par.indexOfWater][2]: # T_inf <= T_mid
    #    a = par.a_low[par.indexOfWater]
    #else:  # T_mid < T_inf
    #    a = par.a_high[par.indexOfWater]
    
    # calculate sum
    #C_V_inf = 0.0
    #for n in range(par.N): # [0, 1, 2, 3, 4]
    #    C_V_inf += a[n] * T_inf**n
    #C_V_inf = par.R_erg * (C_V_inf - 1.0)
# Evaporation energy [J/mol]
    #Old:
    #e_eva = C_V_inf * T_inf * 1e-7
    #e_con = C_V * T * 1e-7
    #evap_energy = n_eva_dot * e_eva - n_con_dot * e_con    # [W/m^2]
    
    #print(m_net_dot)
    #print(H_steam)
    #print(C_p_water*(T_inf-273.15)*J2erg)
    evap_energy = m_net_dot * (H_steam - C_p_water/1000.0*par.W[par.indexOfWater]*J2erg*(T_inf-273.15)) #T_ref=273.15 K where H_water = 0.0
    #C_p_water: J/(kg*K) -> *1000.0/par.W[par.indexOfWater]: J/(mol*K) ->*J2erg: erg/(mol*K)
    #print(evap_energy)
    evap_energy=0.0
    return n_net_dot, evap_energy


"""________________________________Reaction rates________________________________"""

@njit(float64[:](float64, float64[:], float64, float64,float64))
def _forward_rate(T, M_eff, M, p, reaction_rate_threshold):
# Reaction rate
    k_forward = par.A * T ** par.b * np.exp(-par.E / (par.R_cal * T))
    
# Pressure dependent reactions
    Troe_Index=0
    SRI_index=0
    for j, i in enumerate(par.PressureDependentIndexes):    # i is the number of reaction, j is the index of i's place in par.PressureDependentIndexes
        k_inf = k_forward[i]    # par.A[i] * T ** par.b[i] * np.exp(-par.E[i] / (par.R_cal * T))
        k_0 = par.ReacConst[j][0] * T ** par.ReacConst[j][1] * np.exp(-par.ReacConst[j][2] / (par.R_cal * T))
        if i in par.ThirdBodyIndexes:
            thirdBodyIndex = 0
            while par.ThirdBodyIndexes[thirdBodyIndex] != i:
                thirdBodyIndex += 1
            M_eff_loc = M_eff[thirdBodyIndex]
            thirdBodyIndex += 1
        else:
            M_eff_loc = M
        P_r = k_0 / k_inf * M_eff_loc
        
        # Lindemann formalism
        if i in par.LindemannIndexes:
            F = 1.0

        # Troe formalism
        elif i in par.TroeIndexes:
            F_cent = (1.0 - par.Troe[Troe_Index][0]) * np.exp(-T / par.Troe[Troe_Index][1]) + par.Troe[Troe_Index][0] * np.exp(-T / par.Troe[Troe_Index][2]) + np.exp(-par.Troe[Troe_Index][3] / T)
            logF_cent = np.log10(F_cent)
            c2 = -0.4 - 0.67 * logF_cent
            n = 0.75 - 1.27 * logF_cent
            d = 0.14
            logP_r = np.log10(P_r)
            logF = 1.0 / (1.0 + ((logP_r + c2) / (n - d * (logP_r + c2))) ** 2) * logF_cent
            F = 10.0 ** logF
            Troe_Index = Troe_Index + 1
        
        # SRI formalism
        elif i in par.SRIIndexes: 
            X = 1.0 / (1.0 + np.log10(P_r)**2)
            F = par.SRI[SRI_index][3] * (par.SRI[SRI_index][0] * np.exp(-par.SRI[SRI_index][1] / T) + np.exp(-T / par.SRI[SRI_index][2]))**X * T ** par.SRI[SRI_index][4]
            SRI_index = SRI_index + 1
    # Pressure dependent reactions END
    
        k_forward[i] = k_inf * P_r / (1.0 + P_r) * F

# PLOG reactions
    for j, i in enumerate(par.PlogIndexes):
        # determne indexes of the lower and upper pressures
        lower = par.PlogStart[j]
        for k in range(par.PlogStart[j]+1, par.PlogStop[j]-1):  # smallest and largest pressure skipped
            if par.Plog[k][0] < p:
                lower = k
        upper = lower + 1

        # reaction rates at the lower and upper pressures
        k_lower = par.Plog[lower][1] * T ** par.Plog[lower][2] * np.exp(-par.Plog[lower][3] / (par.R_cal * T))
        k_upper = par.Plog[upper][1] * T ** par.Plog[upper][2] * np.exp(-par.Plog[upper][3] / (par.R_cal * T))

        # interpolation
        if par.Plog[par.PlogStart[j]][0] > p:   # p < smallest pressure level
            ln_k=np.log(k_lower)
        elif par.Plog[par.PlogStop[j]-1][0] < p:    # largest pressure level < p
            ln_k=np.log(k_upper)
        else:
            ln_k = np.log(k_lower) + (np.log(p) - np.log(par.Plog[lower][0])) / (np.log(par.Plog[upper][0]) - np.log(par.Plog[lower][0])) * (np.log(k_upper) - np.log(k_lower))
            
        k_forward[i] = np.exp(ln_k)
    return k_forward


@njit(Tuple((float64[:],float64[:]))(float64[:], float64[:], float64[:], float64, float64, int64[:,:]))
def _backward_rate(k_forward, S, H, T, reaction_rate_threshold, reaction_order):
    k_backward = np.zeros((par.I), dtype=np.float64)
    K_c = np.zeros((par.I), dtype=np.float64)
    
    #Forward rate thresholding:
    if(enable_reaction_rate_threshold):
        for i in range(par.I):
            if(abs(k_forward[i]) > reaction_rate_threshold * np.power(par.N_A,reaction_order[i, 0])):
                k_forward[i] = reaction_rate_threshold * np.sign(k_forward[i])
    
    for i in range(par.I):
        DeltaS = 0.0
        DeltaH = 0.0
        for k in range(par.K):
            DeltaS += par.nu[i][k] * S[k]
            DeltaH += par.nu[i][k] * H[k]
        K_p = np.exp(DeltaS / par.R_erg - DeltaH / (par.R_erg * T))
        K_c[i] = K_p * (par.atm2Pa * 10.0 / (par.R_erg * T)) ** np.sum(par.nu[i])
        K_c[i] += (K_c[i] == 0.0) * 1.0e-323  # MODIFIED
        k_backward[i] = k_forward[i] / K_c[i]
    for i in par.IrreversibleIndexes:
        k_backward[i] = 0.0
    
    if(enable_reaction_rate_threshold): 
        for i in range(par.I):
            if(abs(k_backward[i]) > reaction_rate_threshold * np.power(par.N_A,reaction_order[i, 0]) / K_c[i]):
                k_backward[i] = reaction_rate_threshold * np.sign(k_backward[i])
                k_forward[i] = K_c[i]*k_backward[i]
    return k_forward,k_backward

@njit(float64[:](float64, float64[:], float64[:], float64[:], float64, float64, float64))
def _production_rate(T, H, S, c, P_amb, p, M):
    
# Third body correction factors
    M_eff = np.zeros((par.ThirdBodyCount), dtype = np.float64)   # effective total concentration of the third-body 
    for j, i in enumerate(par.ThirdBodyIndexes):
        for k in range(par.K):
            M_eff[j] += par.alfa[j][k] * c[k]
# Forward and backward rates
    reaction_rate_threshold = par.k_B * T / par.h
    k_forward = _forward_rate(T=T, M_eff=M_eff, M=M, p=p, reaction_rate_threshold=reaction_rate_threshold)
    k_forward_limited,k_backward_limited = _backward_rate(k_forward=k_forward, S=S, H=H, T=T, reaction_rate_threshold=reaction_rate_threshold,reaction_order=par.reaction_order)

# Net rates
    q = np.zeros((par.I), dtype = np.float64)
    for i in range(par.I):
        forward = 1.0
        backward = 1.0
        for k in range(par.K):
            forward *= c[k] ** par.nu_forward[i][k]
            backward *= c[k] ** par.nu_backward[i][k]
                
        q[i] = k_forward_limited[i] * forward - k_backward_limited[i] * backward
# Third body reactions
    for j, i in enumerate(par.ThirdBodyIndexes):    # i is the number of reaction, j is the index of i in par.ThirdBodyIndexes
        if i not in par.PressureDependentIndexes:
            q[i] *= M_eff[j]
# Production rates
    omega_dot = np.zeros((par.K), dtype=np.float64)
    for k in range(par.K):
        for i in range(par.I):
            #Check negative concentrations:
            #if(par.nu[i,k] > 0 and c[k]<0):
                #print('i,k,par.nu[i,k],c')
                #print(i)
                #print(k)
                #print(par.nu[i,k])
                #print(c)
                #q[i]=0
            omega_dot[k] += par.nu[i, k] * q[i]
    return omega_dot


"""________________________________Differential equation________________________________"""

@njit(float64[:](float64, float64[:],float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64[:], int64))
def _f(t, x, R_E, P_amb, alfa_M, Gamma, sigma_evap, T_inf, surfactant, P_v, C_4_starred, mu_L, rho_L_ref, p_L_ref, B_L, b_L, Gamma_L, c_L_ref, r_hc, kappa, thermodynamicalcase, ex_args, extra_dims=0):
    """ODE function for the bubble model. Returns the derivative of the state vector x at time t. 
    Use extra_dims to plot extra variables (e.g. energy) during the simulation. Will impact the performance."""   
    R = x[0]      # bubble radius [m]
    R_dot = x[1]  # [m/s]
    if (thermodynamicalcase<2):
        R_dot = 0.0
    T = x[2]      # temperature [K]
    c = x[3:3+par.K]     # molar concentration [mol/cm^3]
    
    M = np.sum(c) # sum of concentration
    X = c / M     # mole fraction [-]
    p = 0.1 * M * par.R_erg * T #P_amb * ((R_E*R_E*R_E-r_hc*r_hc*r_hc)/(R*R*R-r_hc*r_hc*r_hc))**kappa # Partial pressure of the gases [Pa]
    dxdt = np.zeros(x.shape, dtype = np.float64)
    
# d/dt R
    dxdt[0] = R_dot
# Thermodynamics
    (C_p, H, S, C_v) = _thermodynamic(T=T)
    W_avg = C_p_avg = C_v_avg = lambda_avg = 0.0
    for k in range(par.K):
        W_avg += X[k] * par.W[k]
        C_p_avg += X[k] * C_p[k]
        C_v_avg += X[k] * C_v[k]
        lambda_avg += X[k] * par.lambdas[k]

    if enable_heat_transfer:
        rho_avg = W_avg * M # or np.sum(c * par.W)
        chi_avg = 10.0 * lambda_avg * W_avg / (C_p_avg * rho_avg)
        l_th = np.inf
        if R_dot != 0.0:
            l_th = np.sqrt(R * chi_avg / abs(R_dot))
        l_th = min(l_th, R / np.pi)
        Q_th_dot = lambda_avg * (T_inf - T) / l_th
    else:
        Q_th_dot = 0.0
# d/dt c
    if enable_reactions:
        omega_dot = _production_rate(T=T, H=H, S=S, c=c, P_amb=P_amb, p=p, M=M)
    else:
        omega_dot = np.zeros((par.K), dtype = np.float64)
    V = 4.0 * (R*R*R - r_hc*r_hc*r_hc) * np.pi/3.0
    V_dot = 4.0*R*R*R_dot*np.pi
    A = 4.0*R*R*np.pi
    c_dot = omega_dot - c * V_dot/V 
# Evaporation
    if enable_evaporation:
        n_net_dot, evap_energy = _evaporation(p=p, T=T, X_H2O=X[par.indexOfWater], H_steam=H[par.indexOfWater], alfa_M=alfa_M, Gamma=Gamma, sigma_evap=sigma_evap,T_inf=T_inf, P_v=P_v, C_4_starred=C_4_starred, C_p_water=par.C_p_water, J2erg=par.J2erg)
        c_dot[par.indexOfWater] += 1.0e-6 * n_net_dot * 3.0 / R    # water evaporation
    else:
        n_net_dot = evap_energy = 0.0
    
    #Check negative concentrations:
    #for k in range(len(c)):
    #    if c[k]<=0.0:
            #print('c_k lemegy!')
            #print('k,c,c[k], omega_dot,dt:')
            #print(k)
            #print(c)
            #print(c[k])
            #print(c_dot)
            #print(dt)
            #print(-c[k]/dt)
            #if(c_dot[k]<0.0):
                #c_dot[k] = 1.0e-300#0.0
            #print('\nc_dot,jav:\n')
            #print(c_dot)

    dxdt[3:3+par.K] = c_dot
    #for k in range(len(c)):
    #    if c[k]<0:
    #        print('k,c:')
    #        print(k)
    #        print(c)
    #        print('dxdt:')
    #        print(dxdt)
# d/dt T
    sum_omega_dot = np.sum(omega_dot)
    Q_r_dot = 0.0
    for k in range(par.K):
        Q_r_dot -= omega_dot[k] * H[k]
    Q_r_dot += sum_omega_dot * par.R_erg * T
    T_dot = (Q_r_dot + 10.0 * (-p*V_dot/V + Q_th_dot*A/V) + evap_energy/V)/ (M * C_v_avg)
    if(thermodynamicalcase==1):
        T_dot = 0.0
    p_dot = p * (sum_omega_dot/M + T_dot/T - V_dot/V) #-3.0 * p *kappa*R*R*R_dot/(R*R*R - r_hc*r_hc*r_hc)
    
    dxdt[2] = T_dot
    if(thermodynamicalcase==1):
        dxdt[2] = 0.0
# d/dt R_dot
    (Nom, Den) = _pressure(t=t,
        R=R, R_dot=R_dot, mu_L=mu_L, surfactant=surfactant, rho_L_ref=rho_L_ref, p_L_ref=p_L_ref, B_L=B_L, b_L=b_L, Gamma_L=Gamma_L,
        p=p, p_dot=p_dot, P_amb=P_amb, args=ex_args
    )   # delta = (p_L-P_amb) / rho_L_ref
    
    #Nom = (1.0 + R_dot / c_L_ref) * delta + R / c_L_ref * delta_dot - (1.5 - 0.5 * R_dot / c_L_ref) * R_dot ** 2
    #Den = (1.0 - R_dot / c_L_ref) * R + 4.0 * mu_L / (c_L_ref * rho_L_ref)
    
    dxdt[1] = Nom / Den
    if (thermodynamicalcase<2):
        dxdt[1] = 0.0
    
    if enable_dissipated_energy:
        V_dot=4.0 * R * R * R_dot * np.pi
        integrand_th = -(p * (1 + R_dot / c_L_ref) + R / c_L_ref * p_dot) * V_dot
        integrand_v = 16.0 * np.pi * mu_L * (R * R_dot*R_dot + R * R * R_dot * dxdt[1] / c_L_ref)
        integrand_r = 4.0 * np.pi / c_L_ref * R * R * R_dot * (R_dot * p + p_dot * R - 0.5 * rho_L_ref * R_dot * R_dot * R_dot - rho_L_ref * R * R_dot * dxdt[1])

        dxdt[3+par.K] = integrand_th + integrand_v + integrand_r
    else:
        dxdt[3+par.K] = 0.0

# PLOT EXTRA VARIABLES HERE
    if extra_dims > 0:  # You might change this to plot whatever your heart desires
        if enable_dissipated_energy and extra_dims==4:
            dxdt[3+par.K+1] = integrand_th
            dxdt[3+par.K+2] = integrand_v
            dxdt[3+par.K+3] = integrand_r
            dxdt[3+par.K+4] = integrand_th + integrand_v + integrand_r 
    return dxdt

"""________________________________Stop event________________________________"""

@njit(float64(float64, float64[:], float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64[:], int64))
def stop_event(t, x, R_E, P_amb, alfa_M, Gamma, sigma_evap, T_inf, surfactant, P_v, C_4_starred, mu_L, rho_L_ref, p_L_ref, B_L, b_L, Gamma_L, c_L_ref, r_hc, kappa, thermodynamicalcase, ex_args, extra_dims=0):
    # Ha minden derivált abszolút értéke kisebb, mint 10e-10, az esemény bekövetkezik
    if(t>1.0e-3):
        dxdt = _f(t, x, R_E, P_amb, alfa_M, Gamma, sigma_evap, T_inf, surfactant, P_v, C_4_starred, mu_L, rho_L_ref, p_L_ref, B_L, b_L, Gamma_L, c_L_ref, r_hc, kappa, thermodynamicalcase, ex_args, extra_dims=0)
        check=np.zeros(3, dtype=np.float64)
        check[0]=dxdt[0]/(10000.0 * R_E) #10000 Hz as minimal excitation
        check[1]=dxdt[1]/(10000.0**2.0 * R_E) #10000 Hz as minimal excitation
        check[2]=dxdt[2]/(10000.0 * T_inf) #10000 Hz as minimal excitation
        return np.max(np.abs(check)) - 1.0e-8
    else:
        return 1.0



"""________________________________Solving________________________________"""

def solve(cpar, t_int=np.array([0.0, 1.0]), LSODA_timeout=30.0, Radau_timeout=300.0, extra_dims=0, print_errors=False):
    """
    This funfction solves the differential equation, and returns the numerical solution.
    Parameters:
     * cpar: control parameters
     * t_int: time interval
     * LSODA_timeout: timeout for LSODA solver in seconds
     * Radau_timeout: timeout for Radau solver in seconds
     * extra_dims: add extra dimensions to the initial condition array (initial value: 0.0) | 
                   use it to plot extra variables (e.g. energy) during the simulation
     * print_errors: if True, LSODA and Radau errors will be printed during fatal failiures. | 
                     disable JIT to see the exact line the error occured

    Returns:
     * num_sol: numerical solution. Use num_sol.t and num_sol.y to get the time and the solution. Can be None
     * error_code: see de.error_codes: dict, de.get_errors()
     * elapsed_time: elapsed time
    """
    # The event must specify whether it should terminate and whether to search for the roots
    stop_event.terminal = True
    stop_event.direction = 0  # Search for zero-crossing in both directions (up or down)
    
    if type(cpar) == dict:
        cpar = dotdict(cpar)
        
    error_code = 0
    start = time.time()
    num_sol = None
    if not check_cpar(cpar):
        error_code += 300
        return None, error_code, 0.0
    
    ex_args = []
    for key in excitation_args:
        ex_args.append(cpar.get(key, 0.0))
    ex_args = np.array(ex_args, dtype=np.float64)

    IC, lowpressure_error, lowpressure_warning = _initial_condition(cpar, enable_evaporation, extra_dims)
    if lowpressure_error:
        error_code += 100
        return None, error_code, 0.0
    elif lowpressure_warning:
        error_code += 200

    # solving d/dt x=f(t, x, cpar)
    t_eval=np.array([t_int[0],t_int[1]])

    first_step=1.0e-7
    try: # try-catch block
        if(enable_time_evaluation_limit):
            num_sol = func_timeout( # timeout block
                LSODA_timeout, solve_ivp,
                kwargs=dict(fun=_f, t_span=t_int, y0=IC, t_eval=t_eval, method='LSODA', atol = 1e-10, rtol=1e-10, first_step=first_step, events=stop_event,# solve_ivp()'s arguments
                        args=(cpar.R_E, cpar.P_amb, cpar.alfa_M, cpar.Gamma, cpar.sigma_evap, cpar.T_inf, cpar.surfactant, cpar.P_v, cpar.C_4_starred, cpar.mu_L, cpar.rho_L_ref, cpar.p_L_ref, cpar.B_L, cpar.b_L, cpar.Gamma_L, cpar.c_L_ref, cpar.r_hc, cpar.kappa, cpar.thermodynamicalcase, ex_args, extra_dims) # _f()'s arguments
            ))
        else:
            num_sol = func_timeout( # timeout block
                LSODA_timeout, solve_ivp,
                kwargs=dict(fun=_f, t_span=t_int, y0=IC, method='LSODA', atol = 1e-10, rtol=1e-10, first_step=first_step, events=stop_event,# solve_ivp()'s arguments
                       args=(cpar.R_E, cpar.P_amb, cpar.alfa_M, cpar.Gamma, cpar.sigma_evap, cpar.T_inf, cpar.surfactant, cpar.P_v, cpar.C_4_starred, cpar.mu_L, cpar.rho_L_ref, cpar.p_L_ref, cpar.B_L, cpar.b_L, cpar.Gamma_L, cpar.c_L_ref, cpar.r_hc, cpar.kappa, cpar.thermodynamicalcase, ex_args, extra_dims) # _f()'s arguments   
            ))
        if num_sol.success == False:
            error_code += 1
            if print_errors:
                print(colored(f'Error in solve(): LSODE didn\'t converge: ', 'yellow'), num_sol.message)
    except FunctionTimedOut:
        error_code += 2
    except Exception as error:
        error_code += 3
        if print_errors:
            tb = error.__traceback__
            print(colored(f'Error in solve(): LSODE had a fatal error:', 'red'))
            print(''.join(traceback.format_exception(error, error, tb, limit=15)))
    if error_code % 10 != 0:
        try: # try-catch block
            if(enable_time_evaluation_limit):
                num_sol = func_timeout( # timeout block
                    Radau_timeout, solve_ivp, 
                    kwargs=dict(fun=_f, t_span=t_int, y0=IC, t_eval=t_eval, method='Radau', atol = 1e-10, rtol=1e-10, first_step=first_step, events=stop_event,# solve_ivp()'s arguments
                        args=(cpar.R_E, cpar.P_amb, cpar.alfa_M, cpar.Gamma, cpar.sigma_evap, cpar.T_inf, cpar.surfactant, cpar.P_v, cpar.C_4_starred, cpar.mu_L, cpar.rho_L_ref, cpar.p_L_ref, cpar.B_L, cpar.b_L, cpar.Gamma_L, cpar.c_L_ref, cpar.r_hc, cpar.kappa, cpar.thermodynamicalcase, ex_args, extra_dims) # _f()'s arguments
                ))
            else:
                num_sol = func_timeout( # timeout block
                    Radau_timeout, solve_ivp, 
                    kwargs=dict(fun=_f, t_span=t_int, y0=IC, method='Radau', atol = 1e-10, rtol=1e-10, first_step=first_step, events=stop_event,# solve_ivp()'s arguments
                        args=(cpar.R_E, cpar.P_amb, cpar.alfa_M, cpar.Gamma, cpar.sigma_evap, cpar.T_inf, cpar.surfactant, cpar.P_v, cpar.C_4_starred, cpar.mu_L, cpar.rho_L_ref, cpar.p_L_ref, cpar.B_L, cpar.b_L, cpar.Gamma_L, cpar.c_L_ref, cpar.r_hc, cpar.kappa, cpar.thermodynamicalcase, ex_args, extra_dims) # _f()'s arguments
                ))
            if num_sol.success == False:
                error_code += 40
                if print_errors:
                    print(colored(f'Error in solve(): Radau didn\'t converge: ', 'yellow'), num_sol.message)
        except FunctionTimedOut:
            error_code += 50
        except Exception as error:
            error_code += 60
            if print_errors:
                tb = error.__traceback__
                print(colored(f'Error in solve(): Radau had a fatal error:', 'red'))
                print(''.join(traceback.format_exception(error, error, tb, limit=15)))
    end = time.time()
    elapsed_time = (end - start)
    return num_sol, error_code, elapsed_time

# error codes description
error_codes = { # this is also a dictionary
    'xx0': dict(describtion='succecfully solved with LSODA solver', color='green'),
    'xx1': dict(describtion='LSODA solver didn\'t converge', color='yellow'),
    'xx2': dict(describtion='LSODA solver timed out', color='yellow'),
    'xx3': dict(describtion='LSODA solver had a fatal error', color='yellow'),
    'x0x': dict(describtion='succecfully solved with Radau solver', color='green'),
    'x4x': dict(describtion='Radau solver didn\'t converge (NO SOLUTION!)', color='red'),
    'x5x': dict(describtion='Radau solver timed out (NO SOLUTION!)', color='red'),
    'x6x': dict(describtion='Radau solver had a fatal error (NO SOLUTION!)', color='red'),
    '1xx': dict(describtion='Low pressure error: The pressure of the gas is negative', color='red'),
    '2xx': dict(describtion='Low pressure warning: The pressure during the expansion is lower, than the saturated water pressure', color='yellow'),
    '3xx': dict(describtion='Invalid control parameters', color='red'),
}

def get_errors(error_code, printit=False):
    """
    * Input: error_code (int) 
    * Output: list of error codes (str)
    * Also prints colored errors, if printit=True
    """

    # get digits of error_code
    first_digit = f'xx{error_code % 10}'
    second_digit = f'x{(error_code // 10) % 10}x'
    third_digit = f'{(error_code // 100) % 10}xx'

    # get only relevant errors
    errors = []
    errors.append(first_digit)
    if first_digit != 'xx0':
        errors.append(second_digit)
    if third_digit != '0xx':
        errors.append(third_digit)

    # determine if soultion was succesfull
    success = ('1xx' not in errors and '3xx' not in errors) and ('xx0' in errors or 'x0x' in errors)

    # print errors
    if printit:
        for error in errors:
            print(colored(error_codes[error]['describtion'], error_codes[error]['color']))
    return errors, success

"""________________________________Post processing________________________________"""

# This function gets the numerical solution and the control parameters, and returns some datas about the simulation
def get_data(cpar, num_sol, error_code, elapsed_time):
    """This function gets the numerical solution and the control parameters, and returns some datas about the simulation. Arguments:
     * cpar: control parameters (dict or dotdict)
     * num_sol: numerical solution from solve()
     * error_code: error code from solve()
     * elapsed_time: elapsed time from solve()

    Returns:
     * data: dotdict with the post processing data (e.g. collapse time, energy demand, etc.)
    """

    if type(cpar) == dict:
        cpar = dotdict(cpar)
    # copy cpar:
    data = dotdict(dict(
        ID=cpar.ID,
        R_E=cpar.R_E,
        ratio=cpar.ratio,
        P_amb=cpar.P_amb,
        alfa_M=cpar.alfa_M,
        Gamma=cpar.Gamma,
        sigma_evap=cpar.sigma_evap,
        T_inf=cpar.T_inf,
        P_v=cpar.P_v,
        mu_L=cpar.mu_L,
        rho_L_ref=cpar.rho_L_ref,
        c_L_ref=cpar.c_L_ref,
        surfactant=cpar.surfactant,
        gases=cpar.gases,
        fractions=cpar.fractions,
        kappa=cpar.kappa,
        r_hc=cpar.r_hc, 
        Gamma_L=cpar.Gamma_L,
        B_L=cpar.B_L, 
        b_L=cpar.b_L, 
        cV_L=cpar.cV_L, 
        p_L_ref=cpar.p_L_ref,
    ))
    for key in excitation_args:
        data[key] = cpar.get(key, 0.0)
    
    # runtime and error
    data.error_code = error_code
    data.elapsed_time = elapsed_time # [s]
    
    # default values
    data.steps = 0
    #data.collapse_time = 0.0
    #data.T_max = 0.0
    data.x_initial = np.zeros((4+par.K), dtype=np.float64)
    data.x_final = np.zeros((4+par.K), dtype=np.float64)
    data[f'n_{target_specie}'] = 0.0
    data.m_target = 0.0
    data.expansion_work = 0.0
    data.dissipated_acoustic_energy = 0.0
    data.energy_demand = 1.0e30
    data.energy_efficiency = 1.0e30 # legacy for data.energy_demand
    data.enable_heat_transfer = enable_heat_transfer
    data.enable_evaporation = enable_evaporation
    data.enable_reactions = enable_reactions
    data.enable_dissipated_energy = enable_dissipated_energy
    data.enable_reaction_rate_threshold=enable_reaction_rate_threshold
    data.enable_time_evaluation_limit = enable_time_evaluation_limit
    data.excitation_type = excitation_type
    data.target_specie = target_specie
    errors, success = get_errors(error_code)
    data.success = num_sol.success
    if num_sol is None:
        return data
    
    # normal functioning
    data.steps = len(num_sol.t)
    
    from scipy.interpolate import interp1d

    # 1. Loading data
    # Assumption (first three columns): t, R, R_dot
    data2 = np.loadtxt(cpar.file_name, delimiter=',', skiprows=1)
    t_data2 = data2[:, 0]
    R_data2 = data2[:, 1]
    #R_dot_data2 = data2[:, 2]
    # 2. Interpolation
    new_R = interp1d(t_data2, R_data2, kind='linear', fill_value="extrapolate") #np.power((3.0 * V_interp) / (4.0 * np.pi), 1.0/3.0)

    # R_dot = V_dot / (4 * pi * R^2)
    #new_R_dot = interp1d(t_data2, R_dot_data2, kind='linear', fill_value="extrapolate")

    #mask = new_R > 1e-15
    #new_R_dot[mask] = V_dot_interp[mask] / (4.0 * np.pi * np.power(new_R[mask], 2))

    num_sol.y[0, :] = new_R(num_sol.t)
    #num_sol.y[1, :] = new_R_dot(num_sol.t)
    
    data.x_initial = num_sol.y[:, 0] # initial values of [R, R_dot, T, c_1, ... c_K]
    # collapse time (first loc min of R)    TODO fix
    #loc_min = argrelmin(num_sol.y[:][0])
    #data.collapse_time = 0.0
    #if not len(loc_min[0]) == 0:
    #    data.collapse_time = num_sol.t[loc_min[0][0]]
        
    # Energy calculations
    #data.T_max = np.max(num_sol.y[2,:]) # maximum of temperature peaks [K]
    data.x_final = num_sol.y[:,-1] # final values of [R, R_dot, T, c_1, ... c_K]
    last_V = 4.0 / 3.0 * (100.0 * data.x_final[0]) ** 3#((100.0 * data.x_final[0]) ** 3 - (100.0 * cpar.r_hc) ** 3) * np.pi # [cm^3]
    data[f'n_{target_specie}'] = data.x_final[3+par.index[target_specie]] # [mol]
    m_target = 1.0e-3 * data[f'n_{target_specie}'] * par.W[par.index[target_specie]] # [kg]
    data.expansion_work = _work(cpar, enable_evaporation) # [J]
    data.dissipated_acoustic_energy = data.x_final[3+par.K]  # [J]
    all_work = data.expansion_work + data.dissipated_acoustic_energy
    data.energy_demand = 1.0e-6 * all_work / m_target if m_target > 0.0 else 1.0e30 # [MJ/kg]
    data.energy_efficiency = data.energy_demand # legacy for data.energy_demand
    data.target_specie = target_specie
    return data

# keys of data: (except x_final and x_initial)
keys = ['ID', 'R_E', 'ratio', 'P_amb', 'alfa_M', 'Gamma', 'sigma_evap', 'T_inf', 'P_v', 'mu_L', 'rho_L_ref', 'gases', 'fractions', 'surfactant', 'c_L_ref', 'kappa', 'r_hc', 'Gamma_L', 'B_L', 'b_L', 'cV_L', 'p_L_ref',
        'error_code', 'success', 'elapsed_time', 'steps', #'collapse_time', 'T_max', 
        f'n_{target_specie}', 'expansion_work', 'dissipated_acoustic_energy', 'energy_demand',
        'enable_heat_transfer', 'enable_evaporation', 'enable_reactions', 'enable_dissipated_energy', 'excitation_type', 'target_specie','excitation_cycles','ramp_up_cycles','R_and_R_dot_from_file','file_name','const_V'] + excitation_args

def _print_line(name, value, comment, print_it=False):
    """Prints a name=value pair in an organised way, with nicely formatted floats. Arguments:
     * name: name of the parameter
     * value: value of the parameter
     * comment: comment (e.g. the unit of the parameter)
     * print_it: if True, the function will print the text. If False, it will return with the text (sting)"""

    text = f'    {name} = '
    if type(value) == str:
        text += f'\'{value}\','
    elif type(value) == int:
        text += f'{value},'
    elif type(value) == float:
        if abs(value) < 1e-12:
            text += f'0.0,'
        elif abs(value) < 1e-4:
            text += f'{value: .6e},'
        elif abs(value) < 1e-2:
            text += f'{value: .6f},'
        elif abs(value) < 1.0:
            text += f'{value: .4f},'
        else:
            text += f'{value: .2f},'
    elif type(value) == np.ndarray or type(value) == list:
        if name == 'gases':
            gases= ''.join([f'par.index[\'{par.species[i]}\'], ' for i in value])
            text += f'[{gases[:-2]}],'
        if name == 'fractions':
            fractions = ''.join([f'{i}, ' for i in value])
            text += f'[{fractions[:-2]}],'
    else:
        text += f'{value},'

    if print_it:
        print(f'{text: <48} # {comment}')
    else:
        return f'{text: <48} # {comment}\n'

def print_cpar(cpar, without_code=False, print_it=True):
    """Prints the control parameters (cpar) in an organised way. Arguments:
     * cpar: control parameters (dict or dotdict)
     * without_code: if True, an easier to read version is printed. If False, then the result is a valid python code
     * print_it: if True, the function will print the text. If False, it will return with the text (sting)
    """

    text = ''
    if not without_code:
        text += f'cpar = de.dotdict(dict(\n'
    text += _print_line('ID',         int(cpar['ID']),       cpar_rules['ID']['comment'])
    text += f'  # Initial conditions:\n'
    text += _print_line('R_E',        float(cpar['R_E']),    cpar_rules['R_E']['comment'])
    text += _print_line('ratio',      float(cpar['ratio']),  cpar_rules['ratio']['comment'])
    text += _print_line('gases',      cpar['gases'],         cpar_rules['gases']['comment'])
    text += _print_line('fractions',  cpar['fractions'],     cpar_rules['fractions']['comment'])
    text += f'  # Ambient parameters:\n'
    text += _print_line('P_amb',      float(cpar['P_amb']),  cpar_rules['P_amb']['comment'])
    text += _print_line('T_inf',      float(cpar['T_inf']),  cpar_rules['T_inf']['comment'])
    text += f'  # Liquid parameters:\n'
    text += _print_line('alfa_M',     float(cpar['alfa_M']), cpar_rules['alfa_M']['comment'])
    text += _print_line('P_v',        float(cpar['P_v']),    cpar_rules['P_v']['comment'])
    text += _print_line('mu_L',       float(cpar['mu_L']),   cpar_rules['mu_L']['comment'])
    text += _print_line('rho_L_ref',      float(cpar['rho_L_ref']),  cpar_rules['rho_L_ref']['comment'])
    text += _print_line('c_L_ref',        float(cpar['c_L_ref']),    cpar_rules['c_L_ref']['comment'])
    text += _print_line('surfactant', float(cpar['surfactant']), cpar_rules['surfactant']['comment'])
    text += f'  # Excitation parameters: (excitation_type = {excitation_type})\n'
    for arg, unit in zip(excitation_args, excitation_units):
        text += _print_line(arg,      float(cpar[arg]),      f'[{unit}]')
    if not without_code:
        text += f'))\n\n# Calculate pressure/temperature dependent parameters:\n'
        text += f'cpar.mu_L = de.viscosity(cpar.T_inf)\n'
        text += f'cpar.P_v = de.vapour_pressure(cpar.T_inf)\n'
    if print_it:
        print(text)
    else:
        return text

def print_data(cpar, data, print_it=True):
    """Prints the data dictionary from get_data() in an organised way. Arguments:
     * cpar: control parameters
     * data: data dictionary from get_data()
     * print_it: if True, the function will print the text. If False, it will return with the text (sting)
    """

    text = f'Control parameters:\n'
    text += print_cpar(cpar, without_code=True, print_it=False)
    text += f'''\nSimulation info:
    error_code ={data.error_code: .0f} (success = {data.success})
    elapsed_time ={data.elapsed_time: .2f} [s]
    steps ={data.steps: .0f} [-]'''
    
    text += f'''\nFinal state:
    R_final ={1e6*data.x_final[0]: .2f} [um];   R_dot_final ={data.x_final[1]} [m/s];   T_final ={data.x_final[2]: .2f} [K]
    n_{target_specie}_final ={data[f'n_{target_specie}']: .6e} [mol]
    Final molar concentrations: [mol/cm^3]\n        '''
    
    for k, specie in enumerate(par.species):
        text += f'{specie: <6}: {data.x_final[3+k]: 12.6e};    '
        if (k+1) % 4 == 0: text += f'\n        '
    
    text += f'''\nResults:'''
    #collapse_time = {data.collapse_time: .6e} [s]
    #T_max ={data.T_max: .2f} [K]
    text +=f'''
    expansion work = {data.expansion_work: .6e} [J]
    dissipated acoustic energy = {data.dissipated_acoustic_energy: .6e} [J]
    energy demand = {data.energy_demand: .2f} [MJ/kg of {target_specie}]'''
    
    if print_it:
        print(text)
    else:
        return text
    
def simulate(kwargs):
    """This function runs solve() and get_data(), then return with data. 
    Input and output is (or can be) normal dictionary. 
    It is used for multithreading (e.g. in Bruteforce_parameter_sweep.ipynb). 
    The input (kwargs) is a dictionary with the keyword-argument pairs of solve().  
    """

    args = dict(t_int=np.array([0.0, 1.0]), LSODA_timeout=30.0, Radau_timeout=300.0, extra_dims=0)
    for key in kwargs:
        args[key] = kwargs[key]
    args = dotdict(args)
    cpar = dotdict(args.cpar)
    num_sol, error_code, elapsed_time = solve(cpar, args.t_int, LSODA_timeout=args.LSODA_timeout, Radau_timeout=args.Radau_timeout, extra_dims=args.extra_dims)
    data = get_data(cpar, num_sol, error_code, elapsed_time)
    return dict(data)


"""________________________________Plotting________________________________"""

def plot(cpar, t_int=np.array([0.0, 1.0]), n=5.0, base_name='', format='png', LSODA_timeout=30.0, Radau_timeout=300.0,
         presentation_mode=False, plot_pressure=False, plot_extra=False, show_legend=False, show_cpar=True):
    """
    This funfction runs solve(), get_data(), print_data(), and plots the numerical solution. 
    By default, R(t) and T(t) are plotted on the first plot. 
    The amount of chemical species are plotted on the second plot. 
    Optionally, the internal and external pressure can be plotted on the third plot. 
    Parameters:
     * cpar: control parameters in a dictionary
     * t_int: time interval to solve the diffeq in (default: [0, 1] [s]) |   
           graphs will be plotted in this intervall, if not default
     * n: how long should the plotted time interval be compared to the collapse time (default: 5 [-])
     * base_name: save plots as image (default: '' alias do not save) | 
               use base_name='plot' --> plot_1.png, plot_2.png |  
               use base_name='images/plot' to save into images folder |  
               using a folder for images is recommend |  
               this folder have to be created manually
     * format: format of the saved images (available: png, pdf, ps, eps, svg)
     * LSODA_timeout, Radau_timeout: timeout (maximum runtime) for different solvers in solve() in seconds
     * presentation_mode: if True, the plot will be in presentation mode (default: False)
     * plot_pressure: if True, the internal and external pressure will be plotted (default: False)
     * plot_extra: if True, extra dimensions will be plotted (default: False) | 
                   To use this, you have to modify the end of the _f() function and set extra_dims and extra_dim_labels below. 
     * show_legend: if True, the legend will be visible with every single species (default: False)
     * show_cpar: if True, the control parameters will be printed on the plot (default: False)
    """

# Extra dimensions setting. You may change this
    if plot_extra:
        extra_dims = 4
    else:
        extra_dims = 0
    extra_dim_labels = ['Dissipated energy thermal term [J]', 'Dissipated energy viscous term [J]', 'Dissipated energy acoustic term [J]', 'Total dissipated energy [J]']

# Solve
    if type(cpar) == dict:
        cpar = dotdict(cpar)

    num_sol, error_code, elapsed_time = solve(cpar, t_int, LSODA_timeout, Radau_timeout, extra_dims, print_errors=True)
    data = get_data(cpar, num_sol, error_code, elapsed_time)
    
# Print errors
    errors, success = get_errors(error_code, printit=True)
    if num_sol is None:
        print_data(cpar, data)
        return None
    
# Calculations
    #if t_int[1] != 1.0 or not success:
    if(True):
        end_index = -1
    #else:
    #    end_index = np.where(num_sol.t > n * data.collapse_time)[0][0]

    if num_sol.t[end_index] < 1e-3:
        t = num_sol.t[:end_index] * 1e6 # [us]
    else:
        t = num_sol.t[:end_index] * 1e3 # [ms]
    R = num_sol.y[0, :end_index] # [m]
    R_dot = num_sol.y[1, :end_index] # [m/s]
    T = num_sol.y[2, :end_index] # [K]
    #NEW!
    #c = num_sol.y[3:3+par.K, :end_index] # [mol/cm^3]

    V = 4.0 / 3.0 * ((100.0 * R) ** 3 -(100.0 * cpar.r_hc) ** 3) * np.pi # [cm^3]
    #NEW!
    n = num_sol.y[3:3+par.K, :end_index] #c * V
    if plot_pressure:
        internal_pressure = np.sum(n, axis=0) * par.R_g * T / V # [MPa]
    if plot_extra:
        if len(num_sol.y[:, 0]) != 4+par.K+extra_dims or extra_dims != len(extra_dim_labels):
            print(colored('Error! The number of extra dimensions is incorrect. ', 'red'))
            #print(f'Number of dimensions: {len(num_sol.y[:, 0])=}')
            #print(f'Number of dimensions should be: {4+par.K+extra_dims=}')
            #print(f'Number of extra dimensions: {extra_dims=}')
            #print(f'Number of extra dimension labels: {len(extra_dim_labels)=}')
            print_data(cpar, data)
            return None
        extra_plots = num_sol.y[3+par.K : 4+par.K+extra_dims, :end_index]

# plot R and T
    linewidth = 2.0 if presentation_mode else 1.0
    plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
    fig1 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 6))
    ax1 = fig1.add_subplot(axisbelow=True)
    ax2 = ax1.twinx() 
    ax1.plot(t, R / cpar.R_E, color = 'b', linewidth = linewidth)
    ax2.plot(t, T, color = 'r', linewidth = linewidth, linestyle = '-.')

    if num_sol.t[end_index] < 1e-3:
        ax1.set_xlabel('$t$ [μs]')
    else:
        ax1.set_xlabel('$t$ [ms]')
    ax1.set_ylabel('$R/R_E$ [-]', color = 'b')
    ax2.set_ylabel('$T$ [K]', color = 'r')
    if not presentation_mode: ax1.grid()
    
# textbox with initial conditions
    text = f'Initial conditions:\n'
    text += f'    $R_E$ = {1e6*cpar.R_E: .2f} $[\mu m]$\n'
    if cpar.ratio != 1.0:
        text += f'    $R_0/R_E$ = {cpar.ratio: .2f} $[-]$\n'
    text += f'    $P_{{amb}}$ = {1e-5*cpar.P_amb: .2f} $[bar]$\n'
    text += f'    $T_{{inf}}$ = {cpar.T_inf-273.15: .2f} $[°C]$\n'
    text += f'    $P_{{vapour}}$ = {cpar.P_v: .1f} $[Pa]$\n'
    text += f'Initial content:\n    '
    for gas, fraction in zip(cpar.gases, cpar.fractions):
        text += f'{int(100*fraction)}% {par.species[gas]}, ' 
    text = text[:-2] + f'\nExcitation = {excitation_type}:\n'
    for name, unit in zip(excitation_args, excitation_units):
        text += f'    {name} = {cpar[name]: .2f} [{unit}]\n'
    text = text[:-1]

    if show_cpar and not presentation_mode:
        ax2.text(
            0.98, 0.95, # coordinates
            text, transform=ax1.transAxes,
            horizontalalignment='right', verticalalignment='top', multialignment='left',
            fontsize=14, fontstyle='oblique',
            bbox={'facecolor': 'white', 'alpha': 1.0, 'pad': 10},
        )
    
    plt.show()

# plot reactions
    plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
    fig2 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 9))
    ax = fig2.add_subplot(axisbelow=True)

    # plot the lines
        # use this to generate colors:
            # import seaborn as sns
            # colors = sns.color_palette('Set1', n_colors=10)
            # print(colors.as_hex()); colors
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#d48044', '#33adff', '#a65628', '#f781bf', '#d444ca', '#d4ae44']
    color_index = 0
    texts = []
    max_mol = np.max(n, axis=1) # maximum amounts of species [mol]
    indexes_to_plot = []
    for i, specie in enumerate(par.species):
        if n[i, -1] > 1e-24:
            indexes_to_plot.append(i)
    for i, specie in enumerate(par.species):
        name = specie
        for digit in range(10): # turns 'H2O2' into 'H_2O_2'
            name = name.replace(str(digit), '_' + str(digit))
        if i in indexes_to_plot:
            color = colors[color_index]
            color_index = color_index + 1 if color_index < len(colors) - 1 else 0
            linewidth = 2.0 if presentation_mode and n[i, -1] > 1e-24 else 1.0
            ax.plot(t, n[i], linewidth = linewidth, color=color, label = '$' + name + '$') # PLOT HERE
            texts.append((color, name, n[i, -1]))            
        elif not presentation_mode:
            linewidth = 2.0 if presentation_mode else 1.0
            ax.plot(t, n[i], linewidth = linewidth, label = '$' + name + '$')  # PLOT HERE

    # make legend
    texts.sort(key=lambda x: x[-1], reverse=True)
    last_n_final = 1.0e100
    for text in texts:
        color, name, n_final = text
        # spaceing
        if n_final < 1e-24: continue
        limit = 5.5 if presentation_mode else 3.5
        if last_n_final / n_final < limit:
            n_final = last_n_final / limit
        last_n_final = n_final
        # place text
        ax.text(
            t[-1],
            n_final,
            '$' + name + '$',
            color=color,
            fontsize=24 if presentation_mode else 18,
            verticalalignment='center',
            bbox={'facecolor': 'white', 'pad': 0, 'linewidth': 0.0},
        )

    # plot settings
    ax.set_ylim([1e-24, 5.0*max(max_mol)])
    ax.set_yscale('log')
    if num_sol.t[end_index] < 1e-3:
        ax.set_xlabel('$t$ [μs]')
    else:
        ax.set_xlabel('$t$ [ms]')
    ax.set_ylabel('$n_k$ [mol]')
    if not presentation_mode: ax.grid()
    if show_legend: ax.legend()

    plt.show()

# plot pressure excitation and internal pressure
    if plot_pressure:
        plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
        linewidth = 2.0 if presentation_mode else 1.0
        fig3 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 6))
        ax1 = fig3.add_subplot(axisbelow=True)
        ax2 = ax1.twinx()
        ex_args = np.array([cpar[name] for name in excitation_args], dtype=np.float64)
        external_pressure = [1e-3*Excitation(time, cpar.P_amb, ex_args)[0] for time in num_sol.t[:end_index]] # [MPa]
        ax1.plot(t, external_pressure, color='darkorange', label='external pressure', linewidth=linewidth)
        ax2.plot(t, internal_pressure, color='g', label='internal pressure', linewidth=linewidth)

        if num_sol.t[end_index] < 1e-3:
            ax1.set_xlabel('$t$ [μs]')
        else:
            ax1.set_xlabel('$t$ [ms]')
        ax1.set_ylabel('Pressure excitation [kPa]', color='darkorange')
        ax2.set_ylabel('Internal pressure [MPa]', color='g')
        ax2.set_ylim([5e-7*cpar['P_amb'], 2.0*max(internal_pressure)])
        ax2.set_yscale('log')
        if not presentation_mode: ax1.grid()

        plt.show()

# plot extra dimensions
    if plot_extra:
        plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
        linewidth = 2.0 if presentation_mode else 1.0
        fig4 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 6))
        ax = fig4.add_subplot(axisbelow=True)
        for extra_dim_label, extra_plot in zip(extra_dim_labels, extra_plots):
            ax.plot(t, extra_plot, label=extra_dim_label, linewidth=linewidth)

        if num_sol.t[end_index] < 1e-3:
            ax.set_xlabel('$t$ [μs]')
        else:
            ax.set_xlabel('$t$ [ms]')
        ax.set_ylabel('Extra dimensions')
        if not presentation_mode: ax.grid()
        ax.legend()

        plt.show()
    
# saving the plots
    if base_name != '':
        if format not in ['png', 'pdf', 'ps', 'eps', 'svg']:
            print(colored(f'Invalid image format {format}, png is used instead. ','yellow'))
            format = 'png'
        name = base_name + '_radial.' + format
        if os.path.exists(name):
            print(colored(f'Error: {name} already exists. ', 'red'))
            return None
        try:
            if format == 'png':
                metadata = {key: str(data[key]) for key in data.keys()}
            else:
                metadata = {}
            fig1.savefig(base_name+'_radial.'+format, format=format, metadata=metadata, bbox_inches='tight')
            fig2.savefig(base_name+'_molar.'+format, format=format, metadata=metadata, bbox_inches='tight')
            if plot_pressure:
                fig3.savefig(base_name+'_pressure.'+format, format=format, metadata=metadata, bbox_inches='tight')
            if plot_extra:
                fig4.savefig(base_name+'_extra.'+format, format=format, metadata=metadata, bbox_inches='tight')
        except Exception as error:
            print(print(colored(f'Error in saving {base_name}_1.png','red')))
            print(error)

# print data
    print_data(cpar, data)
    return None
           

"""________________________________Save to CSV________________________________"""

def get_settings_and_info():
    """Returns a string with information about the computer, and the settings of the simulation."""

    return f'''General information:
    Created: {datetime.now().strftime("%Y.%m.%d %H:%M:%S (YYYY.MM.DD HH:MM:SS)")}
    Computer name: {socket.gethostname()}
    User name: {os.getlogin()}
    Number of cores: {psutil.cpu_count(logical=False)}
    Number of logical threads: {psutil.cpu_count(logical=True)}
    CPU frequency: {psutil.cpu_freq().max: .2f} MHz
    RAM size: {psutil.virtual_memory().total / 1024**2: .2f} MB

parameters settings:
    model = {par.model}
    species = {par.species}
    number of species = {par.K}
    number of reactions = {par.I}

full_bubble_model settings:
    enable_heat_transfer = {enable_heat_transfer}
    enable_evaporation = {enable_evaporation} 
    enable_reactions = {enable_reactions}
    enable_dissipated_energy = {enable_dissipated_energy}
    enable_reaction_rate_threshold = {enable_reaction_rate_threshold}
    enable_time_evaluation_limit = {enable_time_evaluation_limit}
    target_specie = \'{target_specie}\' # Specie to calculate energy effiqiency
    excitation_type = \'{excitation_type}\' # function to calculate pressure excitation
'''

class Make_dir:
    """Class for saving simulation data to csv files. Constructor arguments:
     * folder_name: name of the folder to save the csv files into (e.g. 'folder', 'folder/subfolder')
     * file_base_name: base name of the csv files (default: 'output_' --> 'output_1.csv', 'output_2.csv', ...)
     * separator: separator character in the csv file (default: ',')"""

    # constructor
    def __init__(self, folder_name, file_base_name='output_', separator=','):
        self.folder_name = folder_name
        self.file_base_name = file_base_name
        self.separator = separator
        self.parent_dir = os.getcwd()
        self.save_dir = os.path.join(self.parent_dir, folder_name)
        self.number = 1    # uniqe ID number for csv files
        self.lines = 0    # number of data lines in currently opened CCSV
        self.is_opened = False    # True, if a CSV is opened
        
        if os.path.exists(self.save_dir):
            self.new = False
            self.number = len([1 for file in os.listdir(self.save_dir) if file[-4:] == '.csv']) + 1
            print(f'Folder already exists with {self.number-1} csv in it')
        else:
            self.new = True
            os.mkdir(self.save_dir)
    
    def _list_to_string(self, array):
        """ Returns a string from a list e.g. [1, 2, 3] -> '1,2,3'"""

        line = ''
        for element in array:
            element = str(element).replace(',', ' ').replace('[', '').replace(']', '')
            if isinstance(element, float):
                line += f'{float(element): e}' + self.separator
            else:     
                line += element + self.separator
        return line[:-1]
        
    def write_line(self, data):
        """Writes the data dict into the currently opened csv as a new line. The line contains the values of the keys in the data dict."""

        line = self._list_to_string([data[key] for key in keys])
        line += self.separator + self._list_to_string([x for x in data['x_initial'][:3+par.K]] + [x for x in data['x_final'][:3+par.K]])
        self.file.write(line + '\n')
        self.lines += 1
        
    def write_solution(self, data, num_sol, file_base_name):
        """Writes the data dict and the numerical solution into two new csv files. Arguments:
         * data: dictionary containing the simulation data from get_data()
         * num_sol: solution of the differential equation from solve()
         * file_base_name: base name of the csv files (e.g. 'name' --> 'name_data.csv', 'name_num_sol.csv')"""

    # create file containing data
        file = os.path.join(self.save_dir, file_base_name + '_data.csv')
        file = open(file, 'w')
        # write header line
        line = self._list_to_string(keys + ['R_0', 'R_dot_0', 'T_0'] + ['c_' + specie + '_0' for specie in par.species] + ['R_last', 'R_dot_last', 'T_last'] + ['c_' + specie + '_last' for specie in par.species])
        file.write(line + '\n')
        # write data
        line = self._list_to_string([data[key]] for key in keys)
        line += self.separator + self._list_to_string([x for x in data.x_initial[:3+par.K]] + [x for x in data.x_final[:3+par.K]])
        file.write(line + '\n')
        file.close()

    # create file containing num_sol
        file = os.path.join(self.save_dir, file_base_name + '_num_sol.csv')
        file = open(file, 'w')
        # write header line
        line = self._list_to_string(['t', 'R', 'R_dot', 'T'] + ['c_' + specie for specie in par.species] + ['dissipated_acoustic_energy']) 
        file.write(line + '\n')
        # write data
        for i in range(len(num_sol.t)):
            line = self._list_to_string([num_sol.t[i]] + list(num_sol.y[:, i]))
            file.write(line + '\n')
        file.close()
    
    def write_string(self, string, file_base_name):
        """Writes the string into a new txt file. Also saves header with get_settings_and_info(). Arguments:
         * string: arbitrary string to write into the file
         * file_base_name: base name of the txt file (e.g. 'name' --> 'name.txt')"""
        
        # create file
        file = os.path.join(self.save_dir, file_base_name + '.txt')
        try:
            file = open(file, 'x')
        except FileExistsError:
            print(colored(f'Error, file \'{file}\' already exists. ', 'red'))
            return None
        # create header
        #file.write(get_settings_and_info())
        # write string
        file.write(string)
        file.close()

    def new_file(self):
        """Creates a new csv file from file_base_name and an unique number. (e.g. 'output_1.csv', 'output_2.csv', ...)
         Automatically closes opened file. Use write_line() to write data into the file."""

        if self.is_opened:
            return None
        self.number = len([1 for file in os.listdir(self.save_dir) if file[-4:] == '.csv']) + 1
        file = os.path.join(self.save_dir, self.file_base_name + str(self.number) + '.csv')
        self.file = open(file, 'w')
        self.is_opened = True
        self.number += 1
        self.lines = 0
        # write header line:
        line = self._list_to_string(keys + ['R_0', 'R_dot_0', 'T_0'] + ['c_' + specie + '_0' for specie in par.species] + ['R_last', 'R_dot_last', 'T_last'] + ['c_' + specie + '_last' for specie in par.species])
        self.file.write(line + '\n')
    
    def close(self):
        """Closes the currently opened csv file."""

        if self.is_opened:
            self.file.close()
            self.is_opened = False
        
    # destructor
    def __del__(self):
        if self.is_opened:
            self.file.close()