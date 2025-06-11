"""Use this script to automatically fill the mechanism folder.
Inputs: INP files from base_folder (see below)
Outputs: Mechanism files in mechanism
"""

import inp_data_extractor as inp
import os

base_folder = './Bubble_dynamics_simulation/INP file examples/'

file_list = [
    ('chemkin_AR_HE_FIXED_by_Cantera.inp',                              'chemkin_ar_he'),
    ('chem_Otomo2018_without_O_FIXED_by_Cantera.inp',                   'chemkin_otomo2018_without_o'),
    ('chem_Otomo2018_FIXED_by_Cantera.inp',                             'chemkin_otomo2018'),
    ('chem_KAUST2023_N2_carbonfree_FIXED_by_Cantera.inp',               'chemkin_kaust2023_n2'),
    ('chem_KAUST2023_N2_carbonfree_FIXED_by_Cantera_without_O.inp',     'chemkin_kaust2023_n2_without_o'),
]

for file, name in file_list:
    file = os.path.join(base_folder, file).replace('\\', '/')
    inp.extract(file, name=name)