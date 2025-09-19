"""
This program extracts data from .inp files and creates parameters.py
Recommended usage:
    importing: import inp_data_extractor as inp
    usage: inp.extract(path)
"""

"""________________________________Libraries________________________________"""

import numpy as np   # matrices, math
import os    # file management
from termcolor import colored   # colored error messages
import argparse   # command line arguments
import molecular_data as data

comment = '!'
invalid_index = 65535

"""________________________________Functions________________________________"""

def  _separate(string: str, separator: str=' ') -> list:
    """Splits string by separator, and removes empty elements from the resulting list. Arguments:
     * string (str): string to be split
     * separator (str): separator for splitting"""
    
    ret = string.split(separator)
    ret.append('')
    while len(ret)>0 and ret[-1] == '':
        ret.remove('')
    return ret

def _find(looking_for: str, lines: list) -> int:
    """Finds the first line containing looking_for. Arguments:
     * looking_for (str): string to be found
     * lines (list): list of strings to be searched in"""

    i = 0
    while i < len(lines) and not looking_for in lines[i]:
        i += 1
    if i >= len(lines):
        print(colored(f'Error in _find(), \'{looking_for}\' not found', 'red'))
        return -1
    return i

def _rearrange(array: list, original: list, new: list) -> list:
    """Changes order of lines in array, so that it will be ordered as new instead of original. Arguments:
     * array (list): array to be rearranged
     * original (list): original order
     * new (list): new order"""

    new_array = []
    if len(array) != len(original) or len(array) != len(new):
        print(colored(f'Warning in _rearrange(), lists have different lengths', 'yellow'))
    for orig in original:
        if not orig in new:
            print(colored(f'Warning in _rearrange(), \'{orig}\' is in original, but not in new', 'yellow'))
    for i in range(len(new)):
        index = 0
        if new[i] in original:
            index = original.index(new[i])
        else:
            print(colored(f'Warning in _rearrange(), \'{new[i]}\' is in new, but not in original', 'yellow'))
        new_array.append(array[index])
    return new_array

def print_array(array, width=0, comments=[], columns=[], max_len=0):
    """Prints a 1D or 2D array or list. Can leave comments above columns and after lines. Arguments:
     * array (list or np.ndarray): array to be printed
     * width (int): width of each element in characters
     * comments (list): list of comments after each line
     * columns (list): list of comments above each column
     * max_len (int): maximum number of elements in one line, if 0, all elements are printed in one line"""
    
    separator = ','
    arr_opener = '{'
    arr_closer = '}'
    line_opener = '{'
    line_closer = '},'
    remark = '//'
    text = ''

    # empty array
    if len(array) == 0:
        return '/*TO' + 'DO: empty array, use nullptr instead.*/'
    if type(array) == np.ndarray:
        array = list(array)
    if type(array[0]) == np.ndarray:
        for i, x in enumerate(array):
            array[i] = list(x)
    if isinstance(array[0], list) and len(array[0]) == 0:
        return '/*TO' + 'DO: empty array, use nullptr instead.*/'

    if isinstance(array[0], list):
    # 2D array
        text += arr_opener + '\n'
        if columns != []:
            text += f'\t{remark}'
            for col in columns:
                text += f'{col: >{width}} '
            text += '\n'
        for i, x in enumerate(array):
            text += '\t' + line_opener
            for y in x:
                if type(y) == str:
                    y = '"' + y + '"'
                if type(y) == bool:
                    y = str(y).lower()
                text += f'{y: >{width}}' + separator
            if i != len(array)-1:
                text = text[:-1]
                text += line_closer
            else:
                if array[0] != []: text = text[:-1]
                text += line_closer
                text = text[:-1]
                text += ' '
            if len(comments) == len(array):
                text += f'    {remark} {comments[i]}'
            text += '\n'
        text += arr_closer
    # 1D array in 1 line
    elif max_len==0:
        if len(comments) == len(array):
            text += remark
            for i, comment in enumerate(comments):
                text += f'{comment: >{width}} '
            text += '\n'
        text += arr_opener
        for i, x in enumerate(array):
            if type(x) == str:
                x = '"' + x + '"'
            if type(x) == bool:
                x = str(x).lower()
            text += f'{x: >{width}}' + separator
        text = text[:-1]
        text += arr_closer
    # 1D array in multiple lines
    else:
        text += arr_opener + '\n\t'
        i = 0
        while i < len(array):
            if len(comments) == len(array):
                text += remark
                for j in range(i, i+max_len):
                    if j >= len(array): break
                    text += f'{comments[j]: >{width}} '
                text += '\n\t '
            for j in range(i, i+max_len):
                if j >= len(array): break
                x = array[j]
                if type(x) == str:
                    x = '\'' + x + '\''
                if type(x) == bool:
                    x = str(x).lower()
                text += f'{x: >{width}}' + separator
            
            i += max_len
            if i < len(array):
                text += '\n\t'
            else:
                text = text[:-1]
            
        text += '\n' + arr_closer
        
    return text
        
"""________________________________Line seperation, remove comments________________________________"""

def _get_lines(text):  
    """Separates text into lines, removes comments and empty lines. """

    text = text.upper()
    text = text.replace('\\\\', '\\')
    text = text.replace('\t', ' ')
    text = text.replace('\r', ' ')
    text = text.replace('\v', '')
    text = text.replace('\b', '')
    text = text.replace('\f', '')
    text = text.replace('\a', '')
    text = text.replace('END', '\nEND')

    all_lines =  _separate(text, '\n')
    if len(all_lines) < 1:
        print(colored(f'Error, can not separate lines. Check if line ends are \'\\n\'', 'red'))
    all_lines.append('END')
    lines = []
    for line in all_lines:
        if line[0]==comment:
            continue
        if comment in line:
            line = line[:line.find(comment)]
        if line.replace(' ', '') != '':
            lines.append(line)
            
    return lines

"""________________________________Elements________________________________"""

def _get_elements(lines):
    """Gets the list of elements. e.g. ['O','H','N','HE','AR']"""

    i = _find('ELEM', lines)
    elements = []
    while not 'END' in lines[i]:
        line = lines[i]
        line = line.replace('ELEMENTS', '')
        line = line.replace('ELEM', '')
        line = line.replace('/', ' ')
        all_elements =  _separate(line, ' ')
        for element in all_elements:
            if element in data.W:
                if element in elements:
                    print(colored(f'Warning, element \'{element}\' is duplacated', 'yellow'))
                else:
                    elements.append(element)
            else:
                if element == '/': continue
                if element.replace('/', '').replace('+', '', 2).replace('-', '', 2).replace('.', '', 1).replace('E', '', 1).isnumeric(): continue
                print(colored(f'Warning, element \'{element}\' is not recognised, it won\'t be, included', 'yellow')) 
        i += 1
    
    return elements

"""________________________________Species data________________________________"""

def _get_species(lines, elements):
    """Gets species, their molar mass and thermal. e.g. ['NH3', 'H2', 'H', 'NH2', 'NH', 'N', 'NNH', 'N2H4', 'N2H3', 'N2H2', 'H2NN', 'N2']"""

  # Get species
    i = _find('SPEC', lines)
    species = []
    while not 'END' in lines[i]:
        line = lines[i]
        line = line.replace('SPECIES', '')
        line = line.replace('SPEC', '')
        species +=  _separate(line, ' ')
        i += 1

# Get W and lambda for species
    W = []
    lambdas = []
    components = np.zeros((len(species), len(elements)), dtype=np.int32)
    digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    for i, specie in enumerate(species):
        if species.count(specie) > 1:
            print(colored(f'Error, specie \'{specie}\' is duplicated', 'red'))
        if specie in data.lambdas:
            lambdas.append(round(data.lambdas[specie], 5))
        else:
            lambdas.append(0.0)
            print(colored(f'Note, lambda value for specie \'{specie}\' is not in data.py: 0.0 is used', 'blue'))

        specie += '  '
        specie = specie.replace('EX', '')
        while len(specie) > 0 and specie[0] != ' ':
            if specie[0:2] in elements:
                if specie[2] in digits:
                    components[i][elements.index(specie[0:2])] += int(specie[2])
                    specie = specie[3:]
                else:
                    components[i][elements.index(specie[0:2])] += 1
                    specie = specie[2:]
            elif specie[0] in elements:
                if specie[1] in digits:
                    components[i][elements.index(specie[0])] += int(specie[1])
                    specie = specie[2:]
                else:
                    components[i][elements.index(specie[0])] += 1
                    specie = specie[1:]
            else:
                print(colored(f'Warning, {specie[0]} is not recognised in specie \'{species[i]}\', it will be neglected', 'yellow'))
                specie = specie[2:]

        w = 0.0
        for element, num in zip(elements, components[i]):
            w += num * data.W[element]
        W.append(round(w, 5))
        
    return species, W, lambdas

"""________________________________Thermodynamic data________________________________"""

def _get_thermo(lines, species):
    """Gets thermodynamic data for species. (temperature range, low and high NASA coefficients)"""

    TempRange = []
    a_low = []
    a_high = []
    materials = []

    i = _find('THER', lines)
    if 'THER' in lines[i]:
        i += 2

    while not 'END' in lines[i]:
        first_line =  _separate(lines[i], ' ')
        if 'TEMP' in first_line:
            i += 1
            while  _separate(lines[i], ' ')[0].replace('+', '', 2).replace('-', '', 2).replace('.', '', 1).replace('E', '', 1).isnumeric():
                i+=1        

        if not first_line[0] in materials:
            materials.append(first_line[0])
            if first_line[-2] == '0':
                first_line = first_line[:-2]
            else:
                first_line = first_line[:-1]
            TempRange.append( [float(first_line[-3]), float(first_line[-2]), float(first_line[-1])] )

            other_lines = [] 
            for j in range(1, 4): 
                line = lines[i+j]
                line = line[:-2]
                line = line.replace('E-', '_')
                line = line.replace('-', ' -')
                line = line.replace('_', 'E-')
                other_lines +=  _separate(line, ' ')

            array_line = [float(num) for num in other_lines]
            a_high.append(array_line[0:7])
            a_low.append(array_line[7:14])

        i += 4

    TempRange = _rearrange(TempRange, materials, species)
    a_low = _rearrange(a_low, materials, species)
    a_high = _rearrange(a_high, materials, species)

    return TempRange, a_low, a_high

"""________________________________Reactions________________________________"""

def _get_reactions(lines, species):
    """Gets reaction data. (reaction equations, Arrhenius parameters, third body constants, pressure dependent reactions: Lindemann, Troe, SRI, PLOG parameters)"""

  # Declare lists
    reactions = []
    numbers = []
    A = []
    B = []
    E = []

    ThirdBodyIndexes = []
    alfa = []

    PressureDependentIndexes = []
    LindemannIndexes = []
    ReacConst = []
    TroeIndexes = []
    Troe = []
    SRIIndexes = []
    SRI = []

    PlogIndexes = []
    Plog = []

  # Get reaction datas
    keywords = ['LOW', 'TROE', 'SRI', 'HIGH', 'REV', 'DUP', 'LT', 'TDEP', 'XSMI', 'PLOG', 'FORD', 'RORD', 'MOME', 'EXCI', 'JAN', '/']
    i = _find('REAC', lines) + 1

    while not 'END' in lines[i]: 
        reaction_line = lines[i]
        reaction_line = reaction_line.replace('<=>', '=')
        reaction_line = reaction_line.replace('=>', '>')
        reaction_line =  _separate(reaction_line, ' ')
        reactions.append(''.join(reaction_line[:-3]).replace('>', '=>'))
        numbers.append(len(reactions)-1)
        A.append(float(reaction_line[-3]))
        B.append(float(reaction_line[-2]))
        E.append(float(reaction_line[-1]))
        
        line = lines[i]
        isThirdBody = isPressureDependent = isTroe = isSRI = isPLOG = False
        isManualThirdBodyCoefficients = False
        
        if '(+M)' in line.replace(' ','').replace('<=>', '=').replace('=>', '>').replace('>', '=>'):
            isPressureDependent = True
            isThirdBody = True
        elif '+M' in line.replace(' ','').replace('<=>', '=').replace('=>', '>').replace('>', '=>'): 
            isThirdBody = True
        
        if not any([keyword in lines[i] for keyword in keywords+['END']]):
            i += 1
            line = lines[i]
        
        while any([keyword in line for keyword in keywords+['END']]): #This cycle steps one reaction.
            if 'END' in line:
                break
            elif 'LOW' in line:
                isPressureDependent = True
                line = line.replace('LOW', '')
                line = line.replace('/', '')
                line =  _separate(line, ' ')
                ReacConst.append([float(line[-3]), float(line[-2]), float(line[-1])])
            elif 'TROE' in line:
                isTroe = True
                line = line.replace('TROE', '')
                line = line.replace('/', '')
                line =  _separate(line, ' ')
                if len(line) >= 4:
                    Troe.append([float(line[-4]), float(line[-3]), float(line[-2]), float(line[-1])])
                else:
                    Troe.append([float(line[-3]), float(line[-2]), float(line[-1]), 1e300]) # last coeff is inf
                    print(colored(f'Note, no T*** in line {i} (\'{lines[i]}\') in reaction \'{reactions[-1]}\': T***=inf is used', 'blue'))
            elif 'SRI' in line:
                isSRI = True
                line = line.replace('SRI', '')
                line = line.replace('/', '')
                line =  _separate(line, ' ')
                if len(line) == 5:
                    SRI.append([float(line[-5]), float(line[-4]), float(line[-3]), float(line[-2]), float(line[-1])])
                else:
                    SRI.append([float(line[-3]), float(line[-2]), float(line[-1]), 1.0, 0.0]) # d=1, e=0
                    print(colored(f'Note, no d or e in line {i} (\'{lines[i]}\') in reaction \'{reactions[-1]}\': d=1, e=0 is used', 'blue'))
            elif 'PLOG' in line:
                if not isPLOG:
                    PlogIndexes.append(numbers[-1])
                    isPLOG = True
                    Plog.append([])
                line = line.replace('PLOG', '')
                line = line.replace('MX', '')
                line = line.replace('SP', '')
                line = line.replace('/', '')
                line =  _separate(line, ' ')
            
                Plog[-1].append([float(line[-4]), float(line[-3]), float(line[-2]), float(line[-1])])
            elif '/' in line:
                line = line.replace('/ ', ' ')
                line = line.replace('/', ' ')
                line =  _separate(line, ' ')
                alfa_line = np.ones((len(species)), dtype=np.float64)
                isManualThirdBodyCoefficients = True
                j = 0
                while j < len(line):
                    if line[j] in species:
                        alfa_line[species.index(line[j])] = round(float(line[j+1]), 5)
                    else:
                        print(colored(f'Warning, third body \'{line[j]}\' is not in species in line {i} (\'{lines[i]}\') in reaction \'{reactions[-1]}\'', 'yellow'))
                    j += 2
                alfa.append(alfa_line)
            elif 'DUP' in line:
                line = line.replace('DUP','')
            else:
                keyword = keywords[[keyword in lines[i] for keyword in keywords].index(True)]
                print(colored(f'Warning, keyword \'{keyword}\' is not supported in line {i} (\'{line}\')', 'yellow'))
                i += 1
                line = lines[i]
            if not any([keyword in line for keyword in keywords+['END']]):
                i += 1
                line = lines[i]

        if isPressureDependent:
            PressureDependentIndexes.append(numbers[-1])
            if isTroe:
                TroeIndexes.append(numbers[-1])
            elif isSRI:
                SRIIndexes.append(numbers[-1])
            else:
                LindemannIndexes.append(numbers[-1])
        if isThirdBody:
            ThirdBodyIndexes.append(numbers[-1])
            if not isManualThirdBodyCoefficients:
                alfa_line = np.ones((len(species)), dtype=np.float64)
                alfa.append(alfa_line)

    for i in LindemannIndexes:
        if i not in PressureDependentIndexes:
            print(colored(f'Error, Lindemann reaction {i} (\'{reactions[i]}\') is not in PressureDependentIndexes', 'red'))
    for i in TroeIndexes:
        if i not in PressureDependentIndexes:
            print(colored(f'Error, Troe reaction {i} (\'{reactions[i]}\') is not in PressureDependentIndexes', 'red'))
    for i in SRIIndexes:
        if i not in PressureDependentIndexes:
            print(colored(f'Error, SRI reaction {i} (\'{reactions[i]}\') is not in PressureDependentIndexes', 'red'))
    for i in PressureDependentIndexes:
        if not (i in LindemannIndexes or i in TroeIndexes or i in SRIIndexes):
            print(colored(f'Error, reaction {i} (\'{reactions[i]}\') is in PressureDependentIndexes but not in LindemannIndexes or TroeIndexes or SRIIndexes', 'red'))

    for i in range(len(Plog)):
        Plog[i] = np.array(Plog[i], dtype=np.float64)

    if alfa == []: alfa = [[]]
    if ReacConst == []: ReacConst = [[]]
    if Troe == []: Troe = [[]]
    if SRI == []: SRI = [[]]
    if Plog == []: Plog = [[[]]]

    A = np.array(A)
    B = np.array(B)
    E = np.array(E)
    ReacConst = np.array(ReacConst)
    
  # Correct units
    i = _find('REAC', lines)
    if len(PlogIndexes) != 0:
        for j in range(len(Plog)):
            Plog[j][:, 0] *= 1.0e5   # convert bar to Pa
    if 'MOLEC' in lines[i]: # MOLECULE
        print(colored(f'Note, pre-exponential factor is modified from units of [cm^3/molecule/s] to [cm^3/mol/s]', 'blue'))
        A /= 6.02214e23
        if len(PressureDependentIndexes) != 0: ReacConst[:, 0] /= 6.02214e23
        if len(PlogIndexes) != 0:
            for j in range(len(Plog)):
                Plog[j][:, 1]  /= 6.02214e23
        # Avogadro's number: N_A = 6.02214e23 [-]
    elif 'MOL' in lines[i]:
        print(colored(f'Note, pre-exponential factor is modified from units of [cm^3/mol/s] to [cm^3/mol/s]', 'blue'))
    if 'CAL' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [cal/mol] to [cal/mol]', 'blue'))
    elif 'KCAL' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [kcal/mol] to [cal/mol]', 'blue'))
        E /= 1000.0 # [kcal/mol -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] /= 1000.0
        if len(PlogIndexes) != 0:
            for j in range(len(Plog)): 
                Plog[j][:, 3]  /= 1000.0
    elif 'JOU' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [J/mol] to [cal/mol]', 'blue'))
        E /= 4.184 # [J/mol -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] /= 4.184
        if len(PlogIndexes) != 0:
            for j in range(len(Plog)):
                Plog[j][:, 3]  /= 4.184
    elif 'KJOU' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [kJ/mol] to [cal/mol]', 'blue'))
        E *= 1000.0 # [kJ/mol -> J/mol]
        E /= 4.184 # [J/mol -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] *= 1000.0 / 4.184
        if len(PlogIndexes) != 0:
            for j in range(len(Plog)):
                Plog[j][:, 3]  *= 1000.0 / 4.184
    elif 'KELV' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [K] to [cal/mol]', 'blue'))
        E *= 1.9872 # [K -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] *= 1.9872
        if len(PlogIndexes) != 0:
            for j in range(len(Plog)):
                Plog[j][:, 3]  *= 1.9872
        # Universal gas constant: impal = 1.987 [cal/mol/K]
    elif 'EVOL' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [eV/mol] to [cal/mol]', 'blue'))
        E *= 1.602176634e-19 # [eV/mol -> J/mol]
        E /= 4.184 # [J/mol -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] *= 1.602176634e-19 / 4.184
        if len(PlogIndexes) != 0:
            for j in range(len(Plog)):
                Plog[j][:, 3]  *= 1.602176634e-19 / 4.184
    
    for i in range(len(E)):
        E[i] = round(E[i], 5)
    
    return (reactions, A, B, E,
            ThirdBodyIndexes, alfa,
            PressureDependentIndexes,
                LindemannIndexes, ReacConst,
                TroeIndexes, Troe,
                SRIIndexes, SRI,
            PlogIndexes, Plog)

"""________________________________Reaction matrixes (nu)________________________________"""

def _get_nu(reactions, species, W):
    """Gets reaction matrixes (nu_forward, nu_backward) and indexes of irreversible reactions"""

    nu_forward = np.zeros((len(reactions), len(species)), dtype=int)
    nu_backward = np.zeros((len(reactions), len(species)), dtype=int)
    reaction_order = np.zeros((len(reactions)), dtype=int)
    IrreversibleIndexes = []
    digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

    for x, reaction in enumerate(reactions):
        if '>' in reaction:
            IrreversibleIndexes.append(x)
        reaction = reaction.replace('+M', '')
        reaction = reaction.replace('(+M)', '')
        reaction = reaction.replace('()', '')
        reaction = reaction.replace('>', '=')
        forward =  _separate(reaction, '=')[0]
        backward =  _separate(reaction, '=')[1]
        forward =  _separate(forward, '+')
        backward =  _separate(backward, '+')
        for f in forward:
            num = 1
            if f[0] in digits:
                num = int(f[0])
                f = f[1:]
            if f in species:
                nu_forward[x][species.index(f)] += num
            else:
                if f =='E':
                    print(colored(f'Warning, electron (E) in reaction {x} (\'{reactions[x]}\') is ignored', 'yellow'))
                elif f == 'HV':
                    print(colored(f'Warning, photon (HV) in reaction {x} (\'{reactions[x]}\') is ignored', 'yellow'))
                else:
                    print(colored(f'Warning, \'{f}\' in reaction {x} (\'{reactions[x]}\') is not in species, it is ignored', 'yellow'))
        reaction_order[x] = int(np.sum(nu_forward[x]))

        for b in backward:
            num = 1
            if b[0] in digits:
                num = int(b[0])
                b = b[1:]
            if b in species:
                nu_backward[x][species.index(b)] += num
            else:
                if b =='E':
                    print(colored(f'Warning, electron (E) in reaction {x} (\'{reactions[x]}\') is ignored', 'yellow'))
                elif b == 'HV':
                    print(colored(f'Warning, photon (HV) in reaction {x} (\'{reactions[x]}\') is ignored', 'yellow'))
                else:
                    print(colored(f'Warning, \'{b}\' in reaction {x} (\'{reactions[x]}\') is not in species, it is ignored', 'yellow'))

    nu = nu_backward - nu_forward
    for i in range(0, len(reactions)):
        if abs(sum(nu[i] * W)) > 1e-5:
            print(colored(f'Warning, nonconsistent reaction {i} (\'{reactions[i]}\')', 'yellow'))

    max_participants = max([sum(abs(nu[i]) + abs(nu_forward[i]) + abs(nu_backward[i]) != 0) for i in range(len(reactions))])
    nu_indexes = np.zeros((len(reactions), max_participants), dtype=np.int32)
    nu_forward_small = np.zeros((len(reactions), max_participants), dtype=np.int32)
    nu_backward_small = np.zeros((len(reactions), max_participants), dtype=np.int32)
    nu_small = np.zeros((len(reactions), max_participants), dtype=np.int32)

    for i in range(len(reactions)):
        loc_nu = []
        loc_nu_forward = []
        loc_nu_backward = []
        indexes = []
        for j in range(len(species)):
            if nu[i][j] != 0 or nu_forward[i][j] != 0 or nu_backward[i][j] != 0:
                loc_nu.append(nu[i][j])
                loc_nu_forward.append(nu_forward[i][j])
                loc_nu_backward.append(nu_backward[i][j])
                indexes.append(j)
        loc_nu =          np.pad(np.array(loc_nu,          dtype=np.int32), (0, max_participants - len(loc_nu)),          'constant', constant_values=(0))
        loc_nu_forward =  np.pad(np.array(loc_nu_forward,  dtype=np.int32), (0, max_participants - len(loc_nu_forward)),  'constant', constant_values=(0))
        loc_nu_backward = np.pad(np.array(loc_nu_backward, dtype=np.int32), (0, max_participants - len(loc_nu_backward)), 'constant', constant_values=(0))
        indexes =         np.pad(np.array(indexes,         dtype=np.int32), (0, max_participants - len(indexes)),         'constant', constant_values=(invalid_index))
        nu_indexes[i] = indexes
        nu_forward_small[i] = loc_nu_forward
        nu_backward_small[i] = loc_nu_backward
        nu_small[i] = loc_nu
    
    return nu_forward, nu_backward, nu, IrreversibleIndexes, max_participants, nu_indexes, nu_forward_small, nu_backward_small, nu_small, reaction_order

"""________________________________Printing to file________________________________"""

def extract(path, name=''):
    """Extracts data from .inp file and creates parameters.py in the current working directory. Arguments:
     * path (str): path to .inp file"""

  # Open file
    print(f'path={path}')
    try:
        file = open(path, 'r', encoding='utf-8', errors='replace')
        text = file.read()
        model = os.path.basename(path)[:-4]
    except:
        print(colored(f'Error, \'{path}\' not found', 'red'))
    if name == '':
        name = model
    name = name.lower()
    define = name.upper() + '_H'
    
  # Extract data
    lines = _get_lines(text)
    elements = _get_elements(lines)
    species, W, lambdas = _get_species(lines, elements)
    TempRange, a_low, a_high = _get_thermo(lines, species)
    (reactions, A, B, E,
        ThirdBodyIndexes, alfa,
        PressureDependentIndexes,
            LindemannIndexes, ReacConst,
            TroeIndexes, Troe,
            SRIIndexes, SRI,
        PlogIndexes, Plog) = _get_reactions(lines, species)
    (nu_forward, nu_backward, nu, IrreversibleIndexes,
      max_participants, nu_indexes, nu_forward_small,
      nu_backward_small, nu_small, reaction_order) = _get_nu(reactions, species, W)

    
  # Create parameters.py
    line_start = '\n// '
    line_end = '\n\n'
    text = ''
    header = ''
    
    # Physical constants, Elements and Species data
    header += f'#ifndef {define}\n#define {define}\n\n#include "parameters.h"\n\n'
    header += f'struct {name}_struct' + '{\n\n'
    text += f'static constexpr char model[] = "{name}";\n'
    text += f'static constexpr char input_file[] = "{model}.inp";\n'
    
    # Species and elements
    text += line_start + 'SPECIES' + line_end
    text += f'static constexpr index_t num_elements = {len(elements)};\n'
    text += f'static constexpr index_t num_species = {len(species)};\n'
    text += f'static constexpr index_t index_of_water = ' + ( str(invalid_index) if not 'H2O' in species else str(species.index('H2O')) ) + ';\n'
    text += f'static constexpr index_t invalid_index = {invalid_index};\n'
    text += f'static constexpr std::pair<const char*, index_t> elements[{len(elements)}] = ' + '{' + ''.join(['{'+f'"{element}", {i}'+'}, ' for i, element in enumerate(elements)])[:-2] + '};\n'
    text += f'static constexpr std::pair<const char*, index_t> species[{len(species)}] = ' + '{' + ''.join(['{'+f'"{specie}", {i}'+'}, ' for i, specie in enumerate(species)])[:-2] + '};\n'
    text += f'static constexpr const char *species_names[] = ' + '{' + ''.join(['"'+f'{specie}'+'", ' for specie in species])[:-2] + '};\n'
    text += f'    //                                           ' + (print_array([species[i] for i in range(len(species))], 12)[1:-1]).replace('",', ' ').replace('"', '  ') + '\n'
    text += f'static constexpr double W[num_species] =        {print_array(W, 12, max_len=0)};\n'
    text += f'static constexpr double lambdas[num_species] =  {print_array(lambdas, 12, max_len=0)};\n\n'
    
    # NASA polynomials
    text += line_start + 'NASA POLYNOMIALS' + line_end
    text += f'static constexpr index_t NASA_order = 5;\n'
    text += f'static constexpr double temp_range[num_species][3] = '+ print_array(TempRange, 8, species, ['T_low', 'T_high', 'T_mid']) + ';\n\n'
    text += f'static constexpr double a_low[num_species][NASA_order+2] = '+ print_array(a_low, 16, species, ['a_1', 'a_2', 'a_3', 'a_4', 'a_5', 'a_6', 'a_7']) + ';\n\n'
    text += f'static constexpr double a_high[num_species][NASA_order+2] = '+ print_array(a_high, 16, species, ['a_1', 'a_2', 'a_3', 'a_4', 'a_5', 'a_6', 'a_7']) + ';\n\n'
    
    # Reaction constants
    text += line_start + 'REACTION CONSTANTS' + line_end
    text += f'static constexpr index_t num_reactions = {len(reactions)};\n'
    text += f'static constexpr double A[num_reactions] = '+ print_array(A, 20, max_len=5) + ';\n\n'
    text += f'static constexpr double b[num_reactions] = '+ print_array(B, 20, max_len=5) + ';\n\n'
    text += f'static constexpr double E[num_reactions] = '+ print_array(E, 20, max_len=5) + ';\n\n'
    text += f'static constexpr index_t reaction_order[num_reactions] = '+ print_array(reaction_order, 4, max_len=10) + ';\n\n'
    
    # Reaction matrixes
    text += line_start + 'REACTION MATRIXES' + line_end
    text += f'static constexpr index_t num_max_specie_per_reaction = {max_participants};\n'
    text += f'static constexpr index_t nu_indexes[num_reactions][num_max_specie_per_reaction] = '+ print_array(nu_indexes, 6, [f'{x:>2}. {reaction}' for x, reaction in enumerate(reactions)]) + ';\n\n'
    text += f'static constexpr stoich_t nu[num_reactions][3][num_max_specie_per_reaction] = ' + '{\n'
    text += f'    // {"nu_forward":>{4*max_participants}}      {"nu_backward":>{4*max_participants}}      {"nu":>{4*max_participants}}\n'
    for i, reaction in enumerate(reactions):
        comma = ',' if i < len(reactions)-1 else ' '
        text += '    { ' + print_array(nu_forward_small[i], 3) + ',    ' + print_array(nu_backward_small[i], 3) + ',    ' + print_array(nu_small[i], 3) + ' }' + comma + '    // ' + f'{i:>2}. {reaction}\n'
    text += '};\n\n'

    # Three-body reactions
    text += line_start + 'THIRD-BODY REACTIONS' + line_end
    text += f'static constexpr index_t num_third_bodies = {len(ThirdBodyIndexes)};\n'
    if len(ThirdBodyIndexes) != 0:
        text += f'static constexpr index_t third_body_indexes[num_third_bodies] =  {print_array(ThirdBodyIndexes, 6)};\n'
        text += f'static constexpr bool is_pressure_dependent[num_third_bodies] = {print_array([bool(i in PressureDependentIndexes) for i in ThirdBodyIndexes], 6)};\n\n'
        text += f'static constexpr double alfa[num_third_bodies][num_species] = '+ print_array(alfa, 8, [f'{x:>2}. {reactions[x]}' for x in ThirdBodyIndexes], species) + ';\n\n'
    else:
        text += f'static constexpr index_t *third_body_indexes = nullptr;\n'
        text += f'static constexpr bool *is_pressure_dependent = nullptr;\n'
        text += f'static constexpr double *alfa = nullptr;\n\n'

    # Irreversible reactions
    text += line_start + 'Irreversible reactions' + line_end
    text += f'static constexpr index_t num_irreversible = {len(IrreversibleIndexes)};\n'
    if len(IrreversibleIndexes) != 0:
        text += f'static constexpr index_t irreversible_indexes[num_irreversible] = {print_array(IrreversibleIndexes, 4)};\n\n'
    else:
        text += f'static constexpr index_t *irreversible_indexes = nullptr;\n\n'

    # Pressure-dependent reactions
    text += line_start + 'Pressure-dependent reactions' + line_end
    text += f'static constexpr index_t num_pressure_dependent = {len(PressureDependentIndexes)};\n'
    text += f'static constexpr index_t num_lindemann = {len(LindemannIndexes)};\n'
    text += f'static constexpr index_t num_troe = {len(TroeIndexes)};\n'
    text += f'static constexpr index_t num_sri = {len(SRIIndexes)};\n'
    if len(PressureDependentIndexes) != 0:
        text += f'static constexpr index_t pressure_dependent_indexes[num_pressure_dependent] = {print_array(PressureDependentIndexes, 4)};\n'
        text += f'static constexpr Parameters::reac_type pressure_dependent_reac_types[num_pressure_dependent] = '
        def reac_type(index):
            if index in LindemannIndexes:
                return 'lindemann_reac'
            elif index in TroeIndexes:
                return 'troe_reac'
            elif index in SRIIndexes:
                return 'sri_reac'
            else:
                print(colored(f'Error, reaction {index} is not in LindemannIndexes, TroeIndexes or SRIIndexes', 'red'))
        text += '{' + ''.join([f'Parameters::reac_type::{reac_type(index)}, ' for index in PressureDependentIndexes])[:-2] + '};\n'
        isThirdBodyIndexes = []
    else:
        text += f'static constexpr index_t *pressure_dependent_indexes = nullptr;\n'
        text += f'static constexpr Parameters::reac_type *pressure_dependent_reac_types = nullptr;\n'
    ThirdBodyIndexes = np.array(ThirdBodyIndexes)
    for index in PressureDependentIndexes:
        if index in ThirdBodyIndexes:
            isThirdBodyIndexes.append(np.where(ThirdBodyIndexes == index)[0][0])
        else:
            isThirdBodyIndexes.append(invalid_index)
    if len(PressureDependentIndexes) != 0:
        text += f'static constexpr index_t is_third_body_indexes[num_pressure_dependent] = {print_array(isThirdBodyIndexes, 6)};\n\n'
        text += f'static constexpr double reac_const[num_pressure_dependent][3] = '+ print_array(ReacConst, 18, [f'{x:>2}. {reactions[x]}' for x in PressureDependentIndexes], ['A_0', 'b_0', 'E_0']) + ';\n\n'
    else:
        text += f'static constexpr index_t *is_third_body_indexes = nullptr;\n'
        text += f'static constexpr double *reac_const = nullptr;\n\n'
    if len(Troe) != 0:
        text += f'static constexpr double troe[{len(Troe)}][4] = '+ print_array(Troe, 18, [f'{x:>2}. {reactions[x]}' for x in TroeIndexes], ['alfa',  'T***',  'T*',  'T**']) + ';\n\n'
    else:
        text += f'static constexpr double *troe = nullptr;\n\n'
    if len(SRI[0]) != 0:
        text += f'static constexpr double sri[{len(SRI)}][5] = '+ print_array(SRI, 18, [f'{x:>2}. {reactions[x]}' for x in SRIIndexes], ['a',  'b',  'c',  'd', 'e']) + ';\n\n'
    else:
        text += f'static constexpr double *sri = nullptr;\n\n'

    PlogSperators = [0]
    PlogFlattened = []
    for i in range(len(Plog)):
        PlogSperators.append(PlogSperators[-1] + len(Plog[i]))
        for PlogLine in Plog[i]:
            PlogFlattened.append(PlogLine)
    PlogFlattened = np.array(PlogFlattened, dtype=np.float64)
    text += f'static constexpr index_t num_plog = {len(PlogIndexes)};\n'
    text += f'static constexpr index_t num_plog_levels = {len(PlogFlattened) if len(PlogIndexes)!=0 else 0};\n'
    if len(PlogIndexes) != 0:
        text += f'static constexpr index_t plog_indexes[num_plog] = {print_array(PlogIndexes, 4)};\n'
    else:
        text += f'static constexpr index_t *plog_indexes = nullptr;\n'
    text += f'static constexpr index_t plog_seperators[{"num_plog_levels+1" if len(PlogIndexes)!=0 else "2"}] = {print_array(PlogSperators, 4)};\n\n'
    PlogComments = [f'{x:>2}. {reactions[x]}' for x in PlogIndexes]
    PlogComments = sum([[comment] + (len(Plog[i])-1) * [''] for i, comment in enumerate(PlogComments)], [])
    if len(PlogIndexes) != 0:
        text += f'static constexpr double plog[num_plog_levels][4] = ' + print_array(PlogFlattened, 18, PlogComments, ['P_1',  'A_1',  'b_1',  'E_1']) + ';\n\n'
    else:
        text += f'static constexpr double *plog = nullptr;\n\n'

    text = text.replace('\t', '    ')
    write_to_file = header
    for line in text.split('\n'):
        if not line.startswith('//'):
            write_to_file += '    ' + line + '\n'
        else:
            write_to_file += line + '\n'
    write_to_file += '}' + f';    // struct {name}_struct\n\n'
    write_to_file += f'#endif   // {define}\n'

    file = open(f'./mechanism/{name}.h', 'w', encoding='utf8')
    file.write(write_to_file)
    file.close()
    
    print(f'model: {name}')
    print(f'File \'{name}.h\' succesfully created')


def main():
    parser = argparse.ArgumentParser(description="Create mechanism .h file from .inp file")
    parser.add_argument('-f', '--file', type=str, help='path to .inp file')
    parser.add_argument('-n', '--name', type=str, help='name of the output file (without extension)')
    args = parser.parse_args()

    # User chooses a file
    if args.file is None:
        # Find .inp files
        inp_files = []
        for root, dirs, files in os.walk('.'):
            for file in files:
                if file.endswith('.inp'):
                    inp_files.append(os.path.join(root, file))
        if not inp_files:
            print('No .inp files found in the current directory.')
            return None
        print('Available .inp files:')
        for idx, file in enumerate(inp_files):
            print(f"{idx: >2}: {file}")

        if inp_files is None:
            return
        choice = int(input("Choose a file by number: "))
        if choice < 0 or choice >= len(inp_files):
            print("Invalid choice.")
            return
        args.file = inp_files[choice]
    # Extract data
    if args.name is not None:
        extract(args.file, name=args.name)
    else:
        extract(args.file)

if __name__ == '__main__':
    main()