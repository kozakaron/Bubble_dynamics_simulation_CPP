"""
This file contains different pressure excitation functions. 
Specify which one to use under "___Settings___" in full_bubble_model.py 

Use getExcitation(excitation_type) to get the excitation function and its arguments. 
Use plotExcitation(excitation_type) to plot the excitation function. 
"""

import numpy as np   # matrices, math
from termcolor import colored   # colored error messages
from numba import njit   # Just In Time compiler
from numba.types import Tuple, float64   # JIT types
from matplotlib import pyplot as plt   # for plotting
import os    # file management

def getExcitation(excitation_type='no_excitation'):
    """
    Returns the pressure excitation function and the arguments it takes. Use plotExcitation() to plot the excitation function.
    Available excitation types:
     * 'no_excitation': constant ambient pressure
     * 'two_sinusoids': two sinusoids with different frequencies and amplitudes, and a phase shift between them
     * 'sin_sqr': sinusoid squared with only n amplitude cycles
     * 'slow_expansion': decrease from abient pressure to min_pressure (decay_time), followed by a slow increase back to ambient pressure (increase_time)
     * 'sin_impulse_flat_ends': sinusoid with only n amplitude cycles, the ends are smoothed out
     * 'sin_impulse': sinusoid with only n amplitude cycles, the ends are not smoothed out
     * 'sin_impulse_logf': same as 'sin_impulse', but log10(freq) is used instead of freq
     * 'double_sin_impulse': n cycle of two sinusoids with different frequencies and amplitudes and no phase shift; the 2nd frequency is freq_ratio times the 1st frequency
     * 'double_sin_impulse_with_phase_shift': n cycle of two sinusoids with different frequencies and amplitudes and phase shift; the 2nd frequency is freq_ratio times the 1st frequency; phase shift is measured using the 1st frequency (which is LOWER)
     * 'multi_sin_impulse': 5 cycles of sinusoids with different frequencies and amplitudes for each cycle
     * 'double_multi_sin_impulse': 5 cycles of two sinusoids with different frequencies and amplitudes for each cycle; 2 sine waves are used each time, similar to 'double_sin_impulse'
     
    Returns:
     * Excitation: jitted excitation function
        \t* inputs: t (time, sec), P_amb (ambient pressure, Pa), args (list of control parameters)
        \t* outputs: p_Inf (excitation pressure, Pa), p_Inf_dot (excitation pressure derivative, Pa/s)
     * args: list of argument names the Excitation function takes
     * units: list of units of the arguments
     * defaults: list of default values for the arguments
    """

    if excitation_type == 'no_excitation':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            return P_amb, 0.0
        
        args = []
        units = []
        defaults = []
        return Excitation, args, units, defaults
    
    elif excitation_type == 'two_sinusoids':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A1, p_A2, freq1, freq2, theta_phase = args
            p_Inf = P_amb + p_A1*np.sin(2.0*np.pi*freq1*t) + p_A2*np.sin(2.0*np.pi*freq2*t + theta_phase) 
            p_Inf_dot = p_A1*2.0*np.pi*freq1*np.cos(2.0*np.pi*freq1*t) + p_A2*2.0*np.pi*freq2*np.cos(2.0*np.pi*freq2*t+theta_phase)
            return p_Inf, p_Inf_dot
        
        args = ['p_A1', 'p_A2', 'freq1', 'freq2', 'theta_phase']
        units = ['Pa', 'Pa', 'Hz', 'Hz', 'rad']
        defaults = [2e5, 1.5e5, 30e3, 60e3, np.pi/3]
        return Excitation, args, units, defaults
    
    elif excitation_type == 'sin_sqr':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A, freq, n = args
            if t < n / freq:
                p_Inf = P_amb + p_A * np.sin(2.0*np.pi*freq*t)**2
                p_Inf_dot = 4.0*p_A*np.pi*freq * np.sin(4.0*np.pi*freq*t)
            else:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            return p_Inf, p_Inf_dot
        
        args = ['p_A', 'freq', 'n']
        units = ['Pa', 'Hz', '-']
        defaults = [-2e5, 30e3, 0.5]
        return Excitation, args, units, defaults
    
    elif excitation_type == 'slow_expansion':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            decay_time, increase_time, min_pressure = args
            if t < 0.0:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t > decay_time+increase_time:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t < decay_time:
                p_Inf = P_amb - (P_amb-min_pressure) / decay_time * t
                p_Inf_dot = - (P_amb-min_pressure) / decay_time
            else:
                p_Inf = min_pressure + (P_amb-min_pressure) / (increase_time) * (t-decay_time)
                p_Inf_dot = (P_amb-min_pressure) / (increase_time)

            return p_Inf, p_Inf_dot
        
        args = ['decay_time', 'increase_time', 'min_pressure']
        units = ['s', 's', 'Pa']
        defaults = [5e-5, 5e-6, 1e4]
        return Excitation, args, units, defaults
    
    elif excitation_type == 'sin_impulse_flat_ends':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A, freq, n = args
            if t < 0.0:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t > n / freq:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t < 0.25 / freq:
                p_Inf = P_amb + p_A*np.sin(2.0*np.pi*freq*t)**2
                p_Inf_dot = 4.0*p_A*np.pi*freq * np.sin(4.0*np.pi*freq*t)
            elif t > (n - 0.25) / freq:
                sign = np.sign( np.sin(2.0*np.pi*freq*t) )
                p_Inf = P_amb + p_A*np.sin(2.0*np.pi*freq*t)**2 * sign
                p_Inf_dot = 4.0*p_A*np.pi*freq * np.sin(4.0*np.pi*freq*t) * sign
            else:
                p_Inf = P_amb + p_A*np.sin(2.0*np.pi*freq*t)
                p_Inf_dot = p_A*2.0*np.pi*freq*np.cos(2.0*np.pi*freq*t)

            return p_Inf, p_Inf_dot
        
        args = ['p_A', 'freq', 'n']
        units = ['Pa', 'Hz', '-']
        defaults = [-2e5, 30e3, 1.0]
        return Excitation, args, units, defaults
    
    elif excitation_type == 'sin_impulse':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A, freq, n = args
            if t < 0.0:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t > n / freq:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            else:
                insin = 2.0*np.pi*freq
                p_Inf = P_amb + p_A*np.sin(insin*t)
                p_Inf_dot = p_A*insin*np.cos(insin*t)

            return p_Inf, p_Inf_dot
        
        args = ['p_A', 'freq', 'n']
        units = ['Pa', 'Hz', '-']
        defaults = [-2e5, 30e3, 1.0]
        return Excitation, args, units, defaults

    elif excitation_type == 'sin_impulse_logf':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A, logf, n = args
            freq = 10**logf
            if t < 0.0:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t > n / freq:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            else:
                insin = 2.0*np.pi*freq
                p_Inf = P_amb + p_A*np.sin(insin*t)
                p_Inf_dot = p_A*insin*np.cos(insin*t)

            return p_Inf, p_Inf_dot
        
        args = ['p_A', 'logf', 'n']
        units = ['Pa', '-', '-']
        defaults = [-2e5, 4.5, 1.0]
        return Excitation, args, units, defaults
    
    elif excitation_type == 'double_sin_impulse':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A1, p_A2, freq, freq_ratio, n = args
            if t < 0.0:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t > n / freq:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            else:
                insin = 2.0*np.pi*freq
                p_Inf = P_amb + p_A1*np.sin(insin*t) + p_A2*np.sin(freq_ratio*insin*t)
                p_Inf_dot = p_A1*insin*np.cos(insin*t) + p_A2*freq_ratio*insin*np.cos(freq_ratio*insin*t)

            return p_Inf, p_Inf_dot
        
        args = ['p_A1', 'p_A2', 'freq', 'freq_ratio', 'n']
        units = ['Pa', 'Pa', 'Hz', '-', '-']
        defaults = [-2e5, -1e5, 30e3, 3.0, 2.0]
        return Excitation, args, units, defaults
    
    elif excitation_type == 'double_sin_impulse_with_phase_shift':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A1, p_A2, freq, freq_ratio, theta_phase, n = args
            if t < 0.0:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t > n / freq:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            else:
                insin = 2.0*np.pi*freq
                if(t<-theta_phase/(2.0*np.pi*freq*freq_ratio) or t>-theta_phase/(2.0*np.pi*freq*freq_ratio) + n/(freq*freq_ratio)):
                    p_Inf = P_amb + p_A1*np.sin(insin*t)
                    p_Inf_dot = p_A1*insin*np.cos(insin*t)
                else:
                    p_Inf = P_amb + p_A1*np.sin(insin*t) + p_A2*np.sin(freq_ratio*insin*t+theta_phase)
                    p_Inf_dot = p_A1*insin*np.cos(insin*t) + p_A2*freq_ratio*insin*np.cos(freq_ratio*insin*t+theta_phase)
            return p_Inf, p_Inf_dot
        
        args = ['p_A1', 'p_A2', 'freq', 'freq_ratio', 'theta_phase', 'n']
        units = ['Pa', 'Pa', 'Hz', '-', '-', '-']
        defaults = [-2e5, -1e5, 30e3, 3.0, 0.0, 2.0]
        return Excitation, args, units, defaults
    
    
    elif excitation_type == 'multi_sin_impulse':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A1, p_A2, p_A3, p_A4, p_A5, freq1, freq2, freq3, freq4, freq5 = args
            if t < 0.0: # 1. wave
                return P_amb, 0.0
            T = 1.0 / freq1
            if t < T:
                insin = 2.0*np.pi*freq1
                p_Inf = P_amb + p_A1*np.sin(insin*t)
                p_Inf_dot = p_A1*insin*np.cos(insin*t)
                return p_Inf, p_Inf_dot
            t -= T
            T = 1.0 / freq2
            if t < T: # 2. wave
                insin = 2.0*np.pi*freq2
                p_Inf = P_amb + p_A2*np.sin(insin*t)
                p_Inf_dot = p_A2*insin*np.cos(insin*t)
                return p_Inf, p_Inf_dot
            t -= T
            T = 1.0 / freq3
            if t < T: # 3. wave
                insin = 2.0*np.pi*freq3
                p_Inf = P_amb + p_A3*np.sin(insin*t)
                p_Inf_dot = p_A3*insin*np.cos(insin*t)
                return p_Inf, p_Inf_dot
            t -= T
            T = 1.0 / freq4
            if t < T: # 4. wave
                insin = 2.0*np.pi*freq4
                p_Inf = P_amb + p_A4*np.sin(insin*t)
                p_Inf_dot = p_A4*insin*np.cos(insin*t)
                return p_Inf, p_Inf_dot
            t -= T
            T = 1.0 / freq5
            if t < T: # 5. wave
                insin = 2.0*np.pi*freq5
                p_Inf = P_amb + p_A5*np.sin(insin*t)
                p_Inf_dot = p_A5*insin*np.cos(insin*t)
                return p_Inf, p_Inf_dot
            return P_amb, 0.0
        
        args = ['p_A1', 'p_A2', 'p_A3', 'p_A4', 'p_A5', 'freq1', 'freq2', 'freq3', 'freq4', 'freq5']
        units = ['Pa', 'Pa', 'Pa', 'Pa', 'Pa', 'Hz', 'Hz', 'Hz', 'Hz', 'Hz']
        defaults = [-1e5, -2e5, -1.5e5, 0.0, 0.0, 50e3, 100e3, 25e3, 10e3, 10e3]
        return Excitation, args, units, defaults
        
    elif excitation_type == 'double_multi_sin_impulse':
        @njit(Tuple((float64, float64))(float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A11, p_A21, p_A31, p_A41, p_A51, p_A12, p_A22, p_A32, p_A42, p_A52, \
            freq1, freq2, freq3, freq4, freq5, freq_ratio1, freq_ratio2, freq_ratio3, freq_ratio4, freq_ratio5 = args
            if t < 0.0: # 1. wave
                return P_amb, 0.0
            T = 1.0 / freq1
            if t < T:
                insin = 2.0*np.pi*freq1
                p_Inf = P_amb + p_A11*np.sin(insin*t) + p_A12*np.sin(freq_ratio1*insin*t)
                p_Inf_dot = p_A11*insin*np.cos(insin*t) + p_A12*freq_ratio1*insin*np.cos(freq_ratio1*insin*t)
                return p_Inf, p_Inf_dot
            t -= T
            T = 1.0 / freq2
            if t < T: # 2. wave
                insin = 2.0*np.pi*freq2
                p_Inf = P_amb + p_A21*np.sin(insin*t) + p_A22*np.sin(freq_ratio2*insin*t)
                p_Inf_dot = p_A21*insin*np.cos(insin*t) + p_A22*freq_ratio2*insin*np.cos(freq_ratio2*insin*t)
                return p_Inf, p_Inf_dot
            t -= T
            T = 1.0 / freq3
            if t < T: # 3. wave
                insin = 2.0*np.pi*freq3
                p_Inf = P_amb + p_A31*np.sin(insin*t) + p_A32*np.sin(freq_ratio3*insin*t)
                p_Inf_dot = p_A31*insin*np.cos(insin*t) + p_A32*freq_ratio3*insin*np.cos(freq_ratio3*insin*t)
                return p_Inf, p_Inf_dot
            t -= T
            T = 1.0 / freq4
            if t < T: # 4. wave
                insin = 2.0*np.pi*freq4
                p_Inf = P_amb + p_A41*np.sin(insin*t) + p_A42*np.sin(freq_ratio4*insin*t)
                p_Inf_dot = p_A41*insin*np.cos(insin*t) + p_A42*freq_ratio4*insin*np.cos(freq_ratio4*insin*t)
                return p_Inf, p_Inf_dot
            t -= T
            T = 1.0 / freq5
            if t < T: # 5. wave
                insin = 2.0*np.pi*freq5
                p_Inf = P_amb + p_A51*np.sin(insin*t) + p_A52*np.sin(freq_ratio5*insin*t)
                p_Inf_dot = p_A51*insin*np.cos(insin*t) + p_A52*freq_ratio5*insin*np.cos(freq_ratio5*insin*t)
                return p_Inf, p_Inf_dot
            return P_amb, 0.0
        
        args = ['p_A11', 'p_A21', 'p_A31', 'p_A41', 'p_A51', 'p_A12', 'p_A22', 'p_A32', 'p_A42', 'p_A52', \
                'freq1', 'freq2', 'freq3', 'freq4', 'freq5', 'freq_ratio1', 'freq_ratio2', 'freq_ratio3', 'freq_ratio4', 'freq_ratio5']
        units = ['Pa', 'Pa', 'Pa', 'Pa', 'Pa', 'Pa', 'Pa', 'Pa', 'Pa', 'Pa', 'Hz', 'Hz', 'Hz', 'Hz', 'Hz', '-', '-', '-', '-', '-']
        defaults = [-2e5, -2.5e5, -2e5, 0.0, 0.0, -1e5, -1.5e5, -1e5, 0.0, 0.0, 50e3, 100e3, 25e3, 10e3, 10e3, 3.0, 2.0, 4.0, 1.0, 1.0]
        return Excitation, args, units, defaults
    
    else:
        print(colored(f'Warning: Excitation excitation_type \'{excitation_type}\' not recognized. Using \'no_excitation\' instead. ', 'yellow'))
        return getExcitation('no_excitation')
    



def plotExcitation(excitation_type, ex_args=None, t_int=[0.0, 1e-4], base_name='', format='png', presentation_mode=False):
    """
    Plots the excitation function.
    Arguments:
     * excitation_type: string, excitation type
     * ex_args: list of control parameters
     * t_int: time interval to plot excitation (default: [0, 0.0001] [s])
     * base_name: save plots as image (default: '' alias do not save) | 
               use base_name='plot' --> plot_1.png, plot_2.png |  
               use base_name='images/plot' to save into images folder |  
               using a folder for images is recommend |  
               this folder have to be created manually
     * format: format of the saved images (available: png, pdf, ps, eps, svg)
     * presentation_mode: if True, the plot will be in presentation mode (default: False)
    """

# Check inputs
    Excitation, args, units, defaults = getExcitation(excitation_type)
    if ex_args is None:
        ex_args = defaults
        print(colored(f'Using default excitation parameters. ', 'white'))
    if type(ex_args) is not np.ndarray:
        try:
            ex_args = np.array(ex_args, dtype=np.float64)
            if len(ex_args) != len(args):
                print(colored(f'Error: excitation parameters must contain {args}. ', 'red'))
                return None
        except:
            print(colored(f'Error: ex_args must be a list or numpy array of floats. ', 'red'))
            return None
    
# Calculations
    P_inf = 1e5
    time = np.linspace(t_int[0], t_int[1], 1000)
    p_inf = np.zeros(time.shape, dtype=np.float64)
    p_inf_dot = np.zeros(time.shape, dtype=np.float64)
    for i, t in enumerate(time):
        p_inf[i], p_inf_dot[i] = Excitation(t, P_inf, ex_args)
    if t_int[1] < 1e-3:
        time *= 1e6
    else:
        time *= 1e3
        
# Plot
    linewidth = 2.0 if presentation_mode else 1.0
    plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
    fig = plt.figure(figsize=(16, 9) if presentation_mode else (20, 6))
    ax1 = fig.add_subplot(axisbelow=True)
    ax2 = ax1.twinx() 
    ax1.plot(time, 1e-5*p_inf, color = 'blue', linewidth = 2.0*linewidth)
    ax2.plot(time, 1e-8*p_inf_dot, color = 'grey', linewidth = linewidth, linestyle = '--')

    if t_int[1] < 1e-3:
        ax1.set_xlabel('$t$ [μs]')
    else:
        ax1.set_xlabel('$t$ [ms]')
    ax1.set_ylabel('$p_\infty \ [bar]$', color = 'blue')
    ax2.set_ylabel('$\\frac{d}{dt}\ p_\infty \ \left[\\frac{kbar}{s}\\right]$', color = 'black')
    ax1.set_title(excitation_type.replace('_', ' ').capitalize())
    if not presentation_mode: ax1.grid()

    plt.show()
    
# Text with excitation parameters
    text = colored(f'Excitation type: {excitation_type}\n', attrs=['bold'])
    for i, arg in enumerate(args):
        if abs(ex_args[i]) < 1e-12:
            arg_text = f'0.0'
        elif abs(ex_args[i]) < 1e-6:
            arg_text = f'{ex_args[i]: .6e}'
        elif abs(ex_args[i]) < 1e-3:
            arg_text = f'{ex_args[i]: .8f}'
        elif abs(ex_args[i]) < 1.0:
            arg_text = f'{ex_args[i]: .4f}'
        else:
            arg_text = f'{ex_args[i]: .1f}'

        text += f'    {arg: <12} = {arg_text: <12} [{units[i]}]'

        if 'log' in arg:
            text += f'    (10^{arg_text} = {10**ex_args[i]: .1f})'
        text += '\n'
    if len(args) > 0:
        text = text[:-1]

    print(text)

# Save as image
    if base_name != '':
        if format not in ['png', 'pdf', 'ps', 'eps', 'svg']:
            print(colored(f'Invalid image format {format}, png is used instead. ','yellow'))
            format = 'png'
        name = base_name + '_pressure_excitation.' + format
        if os.path.exists(name):
            print(colored(f'Error: {name} already exists. ', 'red'))
            return None
        try:
            if format == 'png':
                metadata = {args[i]: str(ex_args[i]) for i in range(len(args))}
            else:
                metadata = {}
            fig.savefig(name, format=format, metadata=metadata, bbox_inches='tight')
        except:
            print(print(colored(f'Error in saving {base_name}_1.png','red')))
