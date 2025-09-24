"""Use this script to automatically fill the mechanism folder.
Inputs: INP files from base_folder (see below)
Outputs: Mechanism files in mechanism
"""

import inp_data_extractor as inp
import os

base_folder1 = './Bubble_dynamics_simulation/INP file examples/'
base_folder2 = './Bubble_dynamics_simulation/temp/'

file_list1 = [
    ('chemkin_AR_HE_FIXED_by_Cantera.inp',                              'chemkin_ar_he'),
    ('chem_Otomo2018_without_O_FIXED_by_Cantera.inp',                   'chemkin_otomo2018_without_o'),
    ('chem_Otomo2018_FIXED_by_Cantera.inp',                             'chemkin_otomo2018'),
    ('chem_KAUST2023_N2_carbonfree_FIXED_by_Cantera.inp',               'chemkin_kaust2023_n2'),
    ('chem_KAUST2023_N2_carbonfree_FIXED_by_Cantera_without_O.inp',     'chemkin_kaust2023_n2_without_o'),
]

file_list2 = [
    ('chem_ELTE2016_ethanol.inp',                                       'chemkin_elte2016_ethanol'),
    ('chem_ELTE2016_syngas.inp',                                        'chemkin_elte2016_syngas'),
    ('chem_ELTE2017_methanol.inp',                                      'chemkin_elte2017_methanol'),
]

for file, name in file_list1:
    #base_folder = os.path.abspath(base_folder)
    file = os.path.join(base_folder1, file).replace('\\', '/')
    inp.extract(file, name=name)
    print('-'*50)

for file, name in file_list2:
    #base_folder = os.path.abspath(base_folder)
    file = os.path.join(base_folder2, file).replace('\\', '/')
    inp.extract(file, name=name)
    print('-'*50)