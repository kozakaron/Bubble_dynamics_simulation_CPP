import cantera as ct
import math
import os

invalid_index = 65535
use_comments = True  # whether to use comments in JSON files
comment = '//'       # comment in JSON files
yaml_dir = './mechanism/yaml_files'
json_dir = './mechanism/json_files'

# colors for terminal output
red =       '\033[91m'
light_red = '\033[31m'
green =     '\033[92m'
yellow =    '\033[93m'
blue =      '\033[94m'
magenta =   '\033[95m'
cyan =      '\033[96m'
grey =      '\033[90m'
reset =     '\033[0m'
white =     '\033[97m'
bold =      '\033[1m'
italic =    '\033[3m'

# thermal conductivity [W/m/K]
lambdas_dict = dict(
    HE=0.151,
    NE=0.0491,
    AR=0.0177,
    KR=0.00943,
    H2=0.1805,
    O2=0.02658,
    H2O=0.016,
    H2O2=0.5863,
    O3=0.019854,
    CO=0.024,
    CO2=0.01663,
    CH3=0.0156,
    CH4=0.03281,
    C2H4=0.020,
    C2H5=0.019,
    C2H6=0.018,
    C3H4=0.015,
    C3H5=0.018,
    C3H6=0.017,
    C3H7=0.016,
    C3H8=0.01674,
    CH2O=0.0165,
    C2OH=0.015,
    C3OH=0.01407,
    CH3OH=0.018,
    N2=0.02583,
    NH3=0.00244,
    NO2=0.00988,
)

def log(message: str, level: str = 'error'):
    """
    Log messages to terminal with color coding based on severity level.
    Levels: 'error', 'warning', 'info', 'debug'
    """

    if level == 'error':
        print(f'{red}[ERROR]{reset} {message}')
    elif level == 'warning':
        print(f'{yellow}[WARNING]{reset} {message}')
    elif level == 'info':
        print(f'{cyan}[INFO]{reset} {message}')
    elif level == 'debug':
        print(f'{cyan}[DEBUG]{reset} {message}')
    else:
        log(f'Unknown log level: {level}', level='error')
        print(message)


def indent(text, num_spaces=4):
    """
    Indents each line of the input text by num_spaces spaces.
    """

    spaces = ' ' * num_spaces
    return '\n'.join([spaces + line for line in text.split('\n')])


def print_element(element, width=0):
    """
    Formats an element index for JSON output. Arguments:
     * element: element
     * width: number of characters per element (for alignment)
    """

    if element is None: element = 'null'
    if isinstance(element, str): element = f'"{element}"'
    if isinstance(element, bool): element = 'true' if element else 'false'
    element = str(element)
    return f'{element: >{width}}'


def print_1D_array(array, width=0, remark='', max_len=0, comma=False):
    """
    Prints 1D list in JSON style. Arguments:
     * array: input array (1D list)
     * width: number of characters per element (for alignment)
     * remark: comment string to append at the end of the line
     * max_len: maximum number of elements to print in one line (0 for all)
     * comma: whether to append a comma at the end of the line
    """

    max_len = max_len if max_len > 0 else max(len(array), 1)
    multi_line = True if max_len < len(array) else False
    array = [
        print_element(element, width) if not isinstance(element, list) else print_1D_array(element, width)
        for element in array
    ]
    lines = [', '.join(array[i:i+max_len]) for i in range(0, len(array), max_len)]
    text = '[\n    ' if multi_line else '['
    text += ',\n    '.join(lines)
    text += '\n]' if multi_line else ']'
    text += ',' if comma else (' ' if remark and use_comments else '')
    text += f'    {comment} {remark}' if remark and use_comments else ''
    
    return text


def print_2D_array(array, width=0, lines=None, columns1=None, columns2=None):
    """
    Prints 2D list in JSON style. An element can be a 1D list. Arguments:
     * array: input array (2D list)
     * width: width of each element in characters
     * lines: list of comments after each line
     * columns1: list of comments above each column for names
     * columns2: list of comments above each column for units
    """

    header = ''
    if columns1 and use_comments:
        header += f'{comment: >6}' + ''.join([f'{col:>{width+2}}' for col in columns1])[2:] + '\n'
    if columns2 and use_comments:
        header += f'{comment: >6}' + ''.join([f'{col:>{width+2}}' for col in columns2])[2:] + '\n'

    rows = []
    for i, row in enumerate(array):
        remark = '' if not lines or i >= len(lines) else lines[i]
        comma = i < len(array) - 1
        rows.append(print_1D_array(row, width=width, remark=remark, comma=comma))

    text = '[\n' + header + '    ' + '\n    '.join(rows) + '\n]'
    return text


def remark(text):
    return f'{comment} {text}' if use_comments else ''


def reaction_order(reaction):
    order1 = 1
    if 'm^3' in str(reaction.rate_coeff_units):
        order1 = 2
    elif 'm^6' in str(reaction.rate_coeff_units):
        order1 = 3
    elif 'm^9' in str(reaction.rate_coeff_units):
        order1 = 4

    order2 = round(sum(reaction.reactants.values()))
    order2 += 1 if reaction.third_body and not 'falloff' in reaction.reaction_type else 0

    if order1 != order2:
        log(f'Reaction order mismatch: {order1} (from units) vs {order2} (from stoichiometry)', level='warning')
    return order1


def transform_arrhenius_parameters(A, b, E, order, reaction):
    if A <= 0.0:
        log(f'Negative or zero A factor: {A} ({reaction})', level='error')
        A = abs(A)

    E_over_R = E / ct.gas_constant
    A *= (1e-3)**(order - 1)   # convert m^(3n-3) / kmol^(n-1) / s --> m^(3n-3) / mol^(n-1) / s
    ln_A = math.log(A)

    return [ln_A, b, E_over_R]


def get_arrhenius_parameters(reaction, order):
    """
    Returns [ln(A), b, E/R] for a given reaction, where:
        * A: pre-exponential factor in m^(3n-3) / mol^(n-1) / s
        * n: reaction order
        * b: temperature exponent
        * E/R: activation energy divided by universal gas constant in K
    """

    if 'pressure-dependent-Arrhenius' in reaction.reaction_type:    # PLOG
        A = b = E = 1.0
    elif 'falloff' in reaction.reaction_type:                       # FALL-OFF
        A = reaction.input_data['high-P-rate-constant']['A']
        b = reaction.input_data['high-P-rate-constant']['b']
        E = reaction.input_data['high-P-rate-constant']['Ea']
    elif 'Arrhenius' in reaction.reaction_type:                   # ARRHENIUS
        A = reaction.input_data['rate-constant']['A']
        b = reaction.input_data['rate-constant']['b']
        E = reaction.input_data['rate-constant']['Ea']
    else:
        log(f'Unknown reaction type: {reaction.reaction_type} ({reaction})', level='error')
        A = b = E = 1.0

    return transform_arrhenius_parameters(A, b, E, order, reaction)


if __name__ == '__main__':
    yaml_files = [file for file in os.listdir(yaml_dir) if file.lower().endswith('.yaml')]
    for yaml_file in yaml_files:
        model = os.path.splitext(yaml_file)[0]
        yaml_path = os.path.join(yaml_dir, yaml_file).replace('\\', '/').replace('//', '/')
        json_file = model + '.json'
        json_path = os.path.join(json_dir, json_file).replace('\\', '/').replace('//', '/')
        print('_' * len(f'Converting {yaml_file}'))
        print(f'{bold}Converting {yaml_file}{reset}')
        mechanism = ct.Solution(yaml_path)


# SPECIES

        lambdas = [lambdas_dict.get(species, 0.0) for species in mechanism.species_names]

        species_text = f'''"num_elements": {mechanism.n_elements},
"num_species": {mechanism.n_species},
"index_of_water": {mechanism.species_index('H2O') if 'H2O' in mechanism.species_names else invalid_index},
"invalid_index": {invalid_index},
"element_names": {print_1D_array(mechanism.element_names, width=5, max_len=10, comma=True)}
                       {remark(print_1D_array(list(range(mechanism.n_species)), width=10))}
"species_names":          {print_1D_array(mechanism.species_names, width=10, comma=True)}
"molar_weights":          {print_1D_array([round(1e-3*sp.molecular_weight, 8) for sp in mechanism.species()], width=10, remark='[kg/mol]', comma=True)}
"thermal_conductivities": {print_1D_array(lambdas, width=10, remark='[W/m/K]', comma=False)}
'''


# NASA7 POLYNOMIALS

        species_comments = [f'{idx: >3}. {species}' for idx, species in enumerate(mechanism.species_names)]
        thermo_list = [species.thermo for species in mechanism.species()]
        for species, th in zip(mechanism.species(), thermo_list):
            if th.n_coeffs != 15:
                log(f'Only the 7 coefficient NASA polynomials and 2 temperature levels are supported: {species} has {th.n_coeffs} coefficients', level='error')

        nasa7_text = f'''"temp_ranges": {print_2D_array(
    [[th.min_temp, th.coeffs[0], th.max_temp] for th in thermo_list],
    width=8, lines=species_comments, columns1=["T_low", "T_mid", "T_high"], columns2=["[K]", "[K]", "[K]"]
)},
"a_low": {print_2D_array(
    [th.coeffs[8:15] for th in thermo_list],
    width=16, lines=species_comments, columns1=[f"a_{i}" for i in range(7)]
)},
"a_high": {print_2D_array(
    [th.coeffs[1:8] for th in thermo_list],
    width=16, lines=species_comments, columns1=[f"a_{i}" for i in range(7)]
)}'''


# ARRHENIUS PARAMETERS

        reaction_orders = [reaction_order(reaction) for reaction in mechanism.reactions()]
        reaction_comments = [f'{idx: >3}. {str(reaction).replace(" ", "")}{f" (duplicate)" if reaction.duplicate else ""}' for idx, reaction in enumerate(mechanism.reactions())]
        arrhenius_params_with_orders = [
            [get_arrhenius_parameters(reaction, order), order] for reaction, order in zip(mechanism.reactions(), reaction_orders)
        ]
        reaction_participants = [set(reaction.products.keys()).union(set(reaction.reactants.keys())) for reaction in mechanism.reactions()]
        num_max_species_per_reaction = max(len(participants) for participants in reaction_participants)

        reaction_participant_indexes = [
            [mechanism.species_index(species) for species in participants] +        # indexes of species participating in the reaction
            [invalid_index] * (num_max_species_per_reaction - len(participants))     # padding with invalid_index
            for participants in reaction_participants
        ]
        for participants in reaction_participant_indexes: participants.sort()

        nu_forward = [[0]*num_max_species_per_reaction for _ in range(len(mechanism.reactions()))]
        nu_backward = [[0]*num_max_species_per_reaction for _ in range(len(mechanism.reactions()))]
        for idx, reaction in enumerate(mechanism.reactions()):
            for species, coeff in reaction.reactants.items():
                species_index = mechanism.species_index(species)
                compression_index = reaction_participant_indexes[idx].index(species_index)
                nu_forward[idx][compression_index] += int(coeff)
            for species, coeff in reaction.products.items():
                species_index = mechanism.species_index(species)
                compression_index = reaction_participant_indexes[idx].index(species_index)
                nu_backward[idx][compression_index] += int(coeff)
        nu = [[nu_backward[i][j] - nu_forward[i][j] for j in range(num_max_species_per_reaction)] for i in range(len(mechanism.reactions()))]
        nu_comment_padding = [''] * (num_max_species_per_reaction - 2)

        arrhenius_text = f'''"num_reactions": {mechanism.n_reactions},
"num_max_species_per_reaction": {num_max_species_per_reaction},
"arrhenius_params_with_orders": {print_2D_array(
    arrhenius_params_with_orders, width=22, lines=reaction_comments, columns1=["ln(A)", "b", "E/R", "order (n)"], columns2=["[m^(3n-3)/mol^(n-1)/s]", "[-]", "[K]", "[-]"]
)},
"stochiometric_coeffs": {print_2D_array(
    [[participants,nf, nb, n] for (participants, nf, nb, n) in zip(reaction_participant_indexes, nu_forward, nu_backward, nu)],
    width=5, lines=reaction_comments,
    columns1=nu_comment_padding + ["", "indexes  "] + nu_comment_padding + ["nu ", "forward "] + nu_comment_padding + ["nu ", "backward  "] + nu_comment_padding + ["", "nu"]
)}'''


# THIRD_BODY REACTIONS

        third_body_indexes = [idx for idx, reaction in enumerate(mechanism.reactions()) if reaction.third_body]
        falloff_indexes = [idx for idx, reaction in enumerate(mechanism.reactions()) if 'falloff' in reaction.reaction_type]
        plog_indexes = [idx for idx, reaction in enumerate(mechanism.reactions()) if 'pressure-dependent-Arrhenius' in reaction.reaction_type]
        third_body_efficiencies = []
        for i, reaction_idx in enumerate(third_body_indexes):
            reaction = mechanism.reaction(reaction_idx)
            third_body_efficiencies.append([reaction.third_body.default_efficiency]*mechanism.n_species)
            for species, efficiency in reaction.third_body.efficiencies.items():
                species_index = mechanism.species_index(species)
                third_body_efficiencies[i][species_index] = efficiency

        third_body_text = f'''"num_third_body_reactions": {len(third_body_indexes)},
"third_body_reaction_indexes": {print_1D_array(third_body_indexes, width=5, comma=True)}
"is_falloff_reaction":         {print_1D_array([idx in falloff_indexes for idx in third_body_indexes], width=5, comma=True)}
"third_body_efficiencies": {print_2D_array(
    third_body_efficiencies, width=8,
    lines=[reaction_comments[idx] for idx in third_body_indexes], columns1=mechanism.species_names
)}'''


# IRREVERSIBLE REACTIONS

        irreversible_indexes = [idx for idx, reaction in enumerate(mechanism.reactions()) if reaction.reversible is False]
        irreversible_text = f'''"num_irreversible_reactions": {len(irreversible_indexes)},
"irreversible_reaction_indexes": {print_1D_array(irreversible_indexes, width=5, comma=False)}
'''


# FALLOFF REACTIONS

        falloff_reactions = [mechanism.reaction(idx) for idx in falloff_indexes]
        lindemann_indexes = [idx for idx in falloff_indexes if mechanism.reaction(idx).reaction_type == 'falloff-Lindemann']
        troe_indexes = [idx for idx in falloff_indexes if mechanism.reaction(idx).reaction_type == 'falloff-Troe']
        troe_reactions = [mechanism.reaction(idx) for idx in troe_indexes]
        sri_indexes = [idx for idx in falloff_indexes if mechanism.reaction(idx).reaction_type == 'falloff-SRI']
        sri_reactions = [mechanism.reaction(idx) for idx in sri_indexes]
        falloff_reaction_types = [
            'lindemann' if idx in lindemann_indexes else
            'troe' if idx in troe_indexes else
            'sri' if idx in sri_indexes else
            log(f'Unknown falloff reaction type: {mechanism.reaction(idx).reaction_type} ({mechanism.reaction(idx)})', level='error')
            for idx in falloff_indexes
        ]
        is_third_body_indexes = [(third_body_indexes.index(idx) if idx in third_body_indexes else invalid_index) for idx in falloff_indexes]
        falloff_parameters = [
            transform_arrhenius_parameters(
                reac.input_data['low-P-rate-constant']['A'],
                reac.input_data['low-P-rate-constant']['b'],
                reac.input_data['low-P-rate-constant']['Ea'],
                reaction_orders[idx]+1, # +1 for (+M) third body, which doesn't affect high-P rate
                reac
            )
            for idx, reac in zip(falloff_indexes, falloff_reactions)
        ]
        troe_parameters = [
            [
                reaction.input_data['Troe']['A'],
                reaction.input_data['Troe'].get('T3', 1.0e-30),
                reaction.input_data['Troe'].get('T1', 1.0e30),
                reaction.input_data['Troe'].get('T2', 1.0e30)
            ]
            for idx, reaction in zip(troe_indexes, troe_reactions)
        ]
        sri_parameters = [
            [
                reaction.input_data['SRI']['A'],
                reaction.input_data['SRI']['B'],
                reaction.input_data['SRI']['C'],
                reaction.input_data['SRI'].get('D', 1.0),
                reaction.input_data['SRI'].get('E', 0.0)
            ]
            for idx, reaction in zip(sri_indexes, sri_reactions)
        ]

        falloff_text = f'''"num_falloff_reactions": {len(falloff_indexes)},
"num_lindemann_reactions": {len(lindemann_indexes)},
"num_troe_reactions": {len(troe_indexes)},
"num_sri_reactions": {len(sri_indexes)},
"falloff_reaction_indexes": {print_1D_array(falloff_indexes, width=10, comma=True)}
"falloff_reaction_types":   {print_1D_array(falloff_reaction_types, width=10, comma=True)}
"is_third_body_indexes":    {print_1D_array(is_third_body_indexes, width=10, comma=True)}
"falloff_parameters": {print_2D_array(
    falloff_parameters, width=22, lines=[reaction_comments[idx] for idx in falloff_indexes],
    columns1=["ln(A)", "b", "E/R"], columns2=["[m^(3n-3)/mol^(n-1)/s]", "[-]", "[K]"]
)},
"troe_parameters": {print_2D_array(
    troe_parameters, width=22, lines=[reaction_comments[idx] for idx in troe_indexes],
    columns1=["alpha", "T***", "T*", "T**"], columns2=["[-]", "[K]", "[K]", "[K]"]
)},
"sri_parameters": {print_2D_array(
    sri_parameters, width=22, lines=[reaction_comments[idx] for idx in sri_indexes],
    columns1=["a", "b", "c", "d", "e"], columns2=["[-]", "[K]", "[K]", "[-]", "[-]"]
)}
'''

# PLOG REACTIONS

        plog_reactions = [mechanism.reaction(idx) for idx in plog_indexes]
        plog_seperators = [0]
        plog_levels = []
        plog_comments = ['']
        for idx, reaction in enumerate(plog_reactions):
            reaction_index = plog_indexes[idx]
            pressures = reaction.input_data['rate-constants']
            pressures.sort(key=lambda x: x['P'])
            plog_seperators.append(plog_seperators[-1] + len(pressures))
            plog_comments[-1] = reaction_comments[reaction_index]
            for rate in pressures:
                P = rate['P']
                A = rate['A']
                b = rate['b']
                E = rate['Ea']
                plog_levels.append([round(P, 1)] + transform_arrhenius_parameters(A, b, E, reaction_orders[reaction_index], reaction))
                plog_comments.append('')

        plog_text = f'''"num_plog_reactions": {len(plog_indexes)},
"num_plog_levels": {len(plog_levels)},
"plog_reaction_indexes": {print_1D_array(plog_indexes, width=10, comma=True)}
"plog_seperators": {print_1D_array(plog_seperators, width=10, comma=True)}
"plog_parameters": {print_2D_array(
    plog_levels, width=22, lines=plog_comments,
    columns1=["P", "ln(A)", "b", "E/R"], columns2=["[Pa]", "[m^(3n-3)/mol^(n-1)/s]", "[-]", "[K]"]
)}'''


# COMBINE ALL TEXT

        json_text = f'''{{
"model": "{model}",
"species": {{
{indent(species_text, 4)}
}},
"thermodynamics": {{
{indent(nasa7_text, 4)}
}},
"arrhenius_parameters": {{
{indent(arrhenius_text, 4)}
}},
"third_body_reactions": {{
{indent(third_body_text, 4)}
}},
"irreversible_reactions": {{
{indent(irreversible_text, 4)}
}},
"falloff_reactions": {{
{indent(falloff_text, 4)}
}},
"plog_reactions": {{
{indent(plog_text, 4)}
}}
}}'''

        if not os.path.exists(json_dir):
            os.makedirs(json_dir)
        with open(json_path, 'w') as json_file:
            json_file.write(json_text)
