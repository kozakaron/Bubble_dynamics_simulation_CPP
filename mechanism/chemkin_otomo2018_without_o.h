#ifndef CHEMKIN_OTOMO2018_WITHOUT_O_H
#define CHEMKIN_OTOMO2018_WITHOUT_O_H

#include "parameters.h"

struct chemkin_otomo2018_without_o_struct{

    static constexpr char model[] = "chemkin_otomo2018_without_o";
    static constexpr char input_file[] = "chem_Otomo2018_without_O_FIXED_by_Cantera.inp";
    
// SPECIES
    
    static constexpr index_t num_elements = 5;
    static constexpr index_t num_species = 12;
    static constexpr index_t index_of_water = 12;
    static constexpr index_t invalid_index = 65535;
    static constexpr std::pair<const char*, index_t> elements[5] = {{"O", 0}, {"H", 1}, {"N", 2}, {"HE", 3}, {"AR", 4}};
    static constexpr std::pair<const char*, index_t> species[12] = {{"NH3", 0}, {"H2", 1}, {"H", 2}, {"NH2", 3}, {"NH", 4}, {"N", 5}, {"NNH", 6}, {"N2H4", 7}, {"N2H3", 8}, {"N2H2", 9}, {"H2NN", 10}, {"N2", 11}};
    static constexpr const char *species_names[] = {"NH3", "H2", "H", "NH2", "NH", "N", "NNH", "N2H4", "N2H3", "N2H2", "H2NN", "N2"};
        //                                                    NH3           H2            H          NH2           NH            N          NNH         N2H4         N2H3         N2H2         H2NN           N2  
    static constexpr double W[num_species] =        {    17.03061,     2.01594,     1.00797,    16.02264,    15.01467,     14.0067,    29.02137,    32.04528,    31.03731,    30.02934,    30.02934,     28.0134};
    static constexpr double lambdas[num_species] =  {     0.00244,      0.1805,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,     0.02583};
    
    
// NASA POLYNOMIALS
    
    static constexpr index_t NASA_order = 5;
    static constexpr double temp_range[num_species][3] = {
        //   T_low   T_high    T_mid 
        {   200.0,  6000.0,  1000.0},    // NH3
        {   200.0,  6000.0,  1000.0},    // H2
        {   200.0,  6000.0,  1000.0},    // H
        {   200.0,  3000.0,  1000.0},    // NH2
        {   200.0,  6000.0,  1000.0},    // NH
        {   200.0,  6000.0,  1000.0},    // N
        {   200.0,  6000.0,  1000.0},    // NNH
        {   300.0,  5000.0,  1000.0},    // N2H4
        {   300.0,  5000.0,  1000.0},    // N2H3
        {   300.0,  5000.0,  1000.0},    // N2H2
        {   300.0,  5000.0,  1695.0},    // H2NN
        {   200.0,  6000.0,  1000.0}     // N2
    };
    
    static constexpr double a_low[num_species][NASA_order+2] = {
        //             a_1              a_2              a_3              a_4              a_5              a_6              a_7 
        {      4.25973162,   -0.0039914953,  1.61261377e-05, -1.41967916e-08,   4.3265074e-12,     -6689.50881,    -0.530747736},    // NH3
        {      2.68763434,   0.00508352924, -1.09134795e-05,  9.75972638e-09, -2.98961778e-12,     -948.720664,    -0.706735279},    // H2
        {             2.5, -1.21434192e-16,  3.14642268e-19,   -3.308722e-22,  1.21191563e-25,        25473.66,     -0.44668285},    // H
        {      4.16323757,  -0.00180348166,   5.9604961e-06, -4.37855755e-09,  1.18864591e-12,      21501.7162,     0.025876001},    // NH2
        {      3.52662396,  2.76384425e-05, -6.48983088e-07,  1.46181571e-09, -6.05450344e-13,      42102.9525,      1.71203374},    // NH
        {      2.51135408, -9.58123439e-05,  2.83259838e-07, -3.43876564e-10,  1.45074992e-13,      56103.6198,      4.14794568},    // N
        {      4.04626045,  -0.00169165832,  8.57759881e-06, -7.01204812e-09,  1.74633997e-12,      28851.0751,      4.12950127},    // NNH
        {    0.0644655703,    0.0274970029, -2.89937434e-05,  1.74515927e-08,  -4.4219888e-12,      10451.9161,      21.2777259},    // N2H4
        {      2.65784686,   0.00872895665,  2.48636948e-06, -6.98011085e-09,  2.63643786e-12,      17322.3998,      9.68717587},    // N2H3
        {      1.61793839,    0.0130635784, -1.71583135e-05,  1.60573529e-08, -6.09410984e-12,      24675.2659,      13.7949214},    // N2H2
        {      2.94840788,   0.00429905889,  1.48887782e-06, -2.14228978e-09,  5.38242847e-13,      33596.0682,      8.68034385},    // H2NN
        {      3.58256851, -0.000558781407,  7.83391479e-07,  8.73631719e-10, -7.49971666e-13,     -1051.60018,      2.75873431}     // N2
    };
    
    static constexpr double a_high[num_species][NASA_order+2] = {
        //             a_1              a_2              a_3              a_4              a_5              a_6              a_7 
        {      2.09566674,   0.00614750045, -2.00328925e-06,  3.01334626e-10, -1.71227204e-14,     -6307.60502,      9.59699937},    // NH3
        {      2.93286575,  0.000826608026, -1.46402364e-07,  1.54100414e-11,   -6.888048e-16,       -816.2239,     -1.02647801},    // H2
        {             2.5,             0.0,             0.0,             0.0,             0.0,        25473.66,     -0.44668285},    // H
        {      2.59263049,   0.00347683597, -1.08271624e-06,  1.49342558e-10, -5.75241187e-15,      21886.8065,      7.90583344},    // NH2
        {      2.78372644,   0.00132985888, -4.24785573e-07,  7.83494442e-11,  -5.5045131e-15,      42345.8847,      5.74063777},    // NH
        {       2.4159429,   0.00017489065,  -1.1902369e-07,   3.0226244e-11,  -2.0360983e-15,      56133.6705,       4.6495384},    // N
        {      3.42744423,   0.00323295234, -1.17296299e-06,  1.90508356e-10, -1.14491506e-14,       28808.692,      6.39339762},    // NNH
        {        4.977317,     0.009595519,   -3.547639e-06,    6.124299e-10,   -4.029795e-14,      9341.22436,      -2.9629862},    // N2H4
        {        4.441846,     0.007214271,   -2.495684e-06,    3.920565e-10,    -2.29895e-14,      16645.2716,     -0.42307347},    // N2H3
        {        3.371185,     0.006039968,   -2.303854e-06,    4.062789e-10,   -2.713144e-14,      24181.7109,      4.98058364},    // N2H2
        {      3.13531032,   0.00568632569, -1.93983467e-06,  3.01290501e-10, -1.74978144e-14,      33364.7095,       7.0448582},    // H2NN
        {      2.95257637,    0.0013969004, -4.92631603e-07,  7.86010195e-11, -4.60755204e-15,     -924.423063,      5.87156477}     // N2
    };
    
    
// REACTION CONSTANTS
    
    static constexpr index_t num_reactions = 35;
    static constexpr double A[num_reactions] = {
                     4.6e+19,             2.2e+16,            640000.0,           1000000.0,                 5.6,
                      9600.0,    70000000000000.0,   100000000000000.0,                0.57,    30000000000000.0,
                1000000000.0,   100000000000000.0,    50000000000000.0,    50000000000000.0,         170000000.0,
                     72000.0,  1500000000000000.0,   560000000000000.0,     7000000000000.0,     3900000000000.0,
                     3.6e+47,         240000000.0,            920000.0,    30000000000000.0,    20000000000000.0,
                     1.8e+40,             85000.0,               0.088,           2400000.0,             3.4e+26,
                     3.4e+26,         480000000.0,    70000000000000.0,           1800000.0,            1.89e+18
    };
    
    static constexpr double b[num_reactions] = {
                        -1.4,                 0.0,                2.39,                2.32,                3.53,
                        2.46,                 0.0,                 0.0,                3.88,                 0.0,
                         0.0,                 0.0,                 0.0,                 0.0,                1.02,
                        1.88,                -0.5,              -0.414,                 0.0,                 0.0,
                      -10.38,                 1.5,                1.94,                 0.0,                 0.0,
                       -8.41,                2.63,                4.05,                 2.0,               -4.83,
                       -4.83,                 1.5,                 0.0,                1.94,               -0.85
    };
    
    static constexpr double E[num_reactions] = {
                    104380.0,             93470.0,             10171.0,               799.0,               552.0,
                       107.0,                 0.0,                 0.0,               342.0,                 0.0,
                         0.0,                 0.0,                 0.0,                 0.0,             11783.0,
                      8802.0,                 0.0,                66.0,              2500.0,              1500.0,
                     69009.0,               -10.0,             -1152.0,                 0.0,                 0.0,
                     73320.0,               230.0,              1610.0,             -1192.0,             46228.0,
                     46228.0,              -894.0,                 0.0,             -1152.0,            224950.0
    };
    
    
// REACTION MATRIXES
    
    static constexpr index_t num_max_specie_per_reaction = 4;
    static constexpr index_t nu_indexes[num_reactions][num_max_specie_per_reaction] = {
        {     1,     2, 65535, 65535},    //  0. H2+M=H+H+M
        {     0,     2,     3, 65535},    //  1. NH3+M=NH2+H+M
        {     0,     1,     2,     3},    //  2. NH3+H=NH2+H2
        {     1,     2,     3,     4},    //  3. NH2+H=NH+H2
        {     0,     3,     4, 65535},    //  4. NH2+NH2=NH3+NH
        {     0,     3,     4,     5},    //  5. NH2+NH=NH3+N
        {     2,     3,     5,    11},    //  6. NH2+N=N2+H+H
        {     1,     2,     4,     5},    //  7. NH+H=N+H2
        {     3,     4,     5, 65535},    //  8. NH+NH=NH2+N
        {     2,     4,     5,    11},    //  9. NH+N=N2+H
        {     2,     6,    11, 65535},    // 10. NNH=N2+H
        {     1,     2,     6,    11},    // 11. NNH+H=N2+H2
        {     3,     4,     6,    11},    // 12. NNH+NH=N2+NH2
        {     0,     3,     6,    11},    // 13. NNH+NH2=N2+NH3
        {     1,     3,     9, 65535},    // 14. NH2+NH2=N2H2+H2
        {     1,     3,    10, 65535},    // 15. NH2+NH2=H2NN+H2
        {     2,     3,     4,     9},    // 16. NH2+NH=N2H2+H
        {     3,     7, 65535, 65535},    // 17. NH2+NH2(+M)=N2H4(+M)
        {     1,     2,     7,     8},    // 18. N2H4+H=N2H3+H2
        {     0,     3,     7,     8},    // 19. N2H4+NH2=N2H3+NH3
        {     2,     8,     9, 65535},    // 20. N2H3=N2H2+H
        {     1,     2,     8,     9},    // 21. N2H3+H=N2H2+H2
        {     0,     3,     8,     9},    // 22. N2H3+NH2=N2H2+NH3
        {     0,     3,     8,    10},    // 23. N2H3+NH2=H2NN+NH3
        {     3,     4,     8,     9},    // 24. N2H3+NH=N2H2+NH2
        {     2,     6,     9, 65535},    // 25. N2H2=NNH+H
        {     1,     2,     6,     9},    // 26. N2H2+H=NNH+H2
        {     0,     3,     6,     9},    // 27. N2H2+NH2=NNH+NH3
        {     3,     4,     6,     9},    // 28. N2H2+NH=NNH+NH2
        {     2,     6,    10, 65535},    // 29. H2NN=NNH+H
        {     2,     6,    10, 65535},    // 30. H2NN=NNH+H
        {     1,     2,     6,    10},    // 31. H2NN+H=NNH+H2
        {     2,     9,    10, 65535},    // 32. H2NN+H=N2H2+H
        {     0,     3,     6,    10},    // 33. H2NN+NH2=NNH+NH3
        {     5,    11, 65535, 65535}     // 34. N2+M=N+N+M
    };
    
    static constexpr stoich_t nu[num_reactions][3][num_max_specie_per_reaction] = {
        //       nu_forward           nu_backward                    nu
        { {  1,  0,  0,  0},    {  0,  2,  0,  0},    { -1,  2,  0,  0} },    //  0. H2+M=H+H+M
        { {  1,  0,  0,  0},    {  0,  1,  1,  0},    { -1,  1,  1,  0} },    //  1. NH3+M=NH2+H+M
        { {  1,  0,  1,  0},    {  0,  1,  0,  1},    { -1,  1, -1,  1} },    //  2. NH3+H=NH2+H2
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    //  3. NH2+H=NH+H2
        { {  0,  2,  0,  0},    {  1,  0,  1,  0},    {  1, -2,  1,  0} },    //  4. NH2+NH2=NH3+NH
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    //  5. NH2+NH=NH3+N
        { {  0,  1,  1,  0},    {  2,  0,  0,  1},    {  2, -1, -1,  1} },    //  6. NH2+N=N2+H+H
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    //  7. NH+H=N+H2
        { {  0,  2,  0,  0},    {  1,  0,  1,  0},    {  1, -2,  1,  0} },    //  8. NH+NH=NH2+N
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    //  9. NH+N=N2+H
        { {  0,  1,  0,  0},    {  1,  0,  1,  0},    {  1, -1,  1,  0} },    // 10. NNH=N2+H
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 11. NNH+H=N2+H2
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 12. NNH+NH=N2+NH2
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 13. NNH+NH2=N2+NH3
        { {  0,  2,  0,  0},    {  1,  0,  1,  0},    {  1, -2,  1,  0} },    // 14. NH2+NH2=N2H2+H2
        { {  0,  2,  0,  0},    {  1,  0,  1,  0},    {  1, -2,  1,  0} },    // 15. NH2+NH2=H2NN+H2
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 16. NH2+NH=N2H2+H
        { {  2,  0,  0,  0},    {  0,  1,  0,  0},    { -2,  1,  0,  0} },    // 17. NH2+NH2(+M)=N2H4(+M)
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 18. N2H4+H=N2H3+H2
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 19. N2H4+NH2=N2H3+NH3
        { {  0,  1,  0,  0},    {  1,  0,  1,  0},    {  1, -1,  1,  0} },    // 20. N2H3=N2H2+H
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 21. N2H3+H=N2H2+H2
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 22. N2H3+NH2=N2H2+NH3
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 23. N2H3+NH2=H2NN+NH3
        { {  0,  1,  1,  0},    {  1,  0,  0,  1},    {  1, -1, -1,  1} },    // 24. N2H3+NH=N2H2+NH2
        { {  0,  0,  1,  0},    {  1,  1,  0,  0},    {  1,  1, -1,  0} },    // 25. N2H2=NNH+H
        { {  0,  1,  0,  1},    {  1,  0,  1,  0},    {  1, -1,  1, -1} },    // 26. N2H2+H=NNH+H2
        { {  0,  1,  0,  1},    {  1,  0,  1,  0},    {  1, -1,  1, -1} },    // 27. N2H2+NH2=NNH+NH3
        { {  0,  1,  0,  1},    {  1,  0,  1,  0},    {  1, -1,  1, -1} },    // 28. N2H2+NH=NNH+NH2
        { {  0,  0,  1,  0},    {  1,  1,  0,  0},    {  1,  1, -1,  0} },    // 29. H2NN=NNH+H
        { {  0,  0,  1,  0},    {  1,  1,  0,  0},    {  1,  1, -1,  0} },    // 30. H2NN=NNH+H
        { {  0,  1,  0,  1},    {  1,  0,  1,  0},    {  1, -1,  1, -1} },    // 31. H2NN+H=NNH+H2
        { {  1,  0,  1,  0},    {  1,  1,  0,  0},    {  0,  1, -1,  0} },    // 32. H2NN+H=N2H2+H
        { {  0,  1,  0,  1},    {  1,  0,  1,  0},    {  1, -1,  1, -1} },    // 33. H2NN+NH2=NNH+NH3
        { {  0,  1,  0,  0},    {  2,  0,  0,  0},    {  2, -1,  0,  0} }     // 34. N2+M=N+N+M
    };
    
    
// THIRD-BODY REACTIONS
    
    static constexpr index_t num_third_bodies = 4;
    static constexpr index_t third_body_indexes[num_third_bodies] =  {     0,     1,    17,    34};
    static constexpr bool is_pressure_dependent[num_third_bodies] = { false, false,  true, false};
    
    static constexpr double alfa[num_third_bodies][num_species] = {
        //     NH3       H2        H      NH2       NH        N      NNH     N2H4     N2H3     N2H2     H2NN       N2 
        {     1.0,     2.5,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0},    //  0. H2+M=H+H+M
        {     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0},    //  1. NH3+M=NH2+H+M
        {     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0},    // 17. NH2+NH2(+M)=N2H4(+M)
        {     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0}     // 34. N2+M=N+N+M
    };
    
    
// Irreversible reactions
    
    static constexpr index_t num_irreversible = 0;
    static constexpr index_t *irreversible_indexes = nullptr;
    
    
// Pressure-dependent reactions
    
    static constexpr index_t num_pressure_dependent = 1;
    static constexpr index_t num_lindemann = 0;
    static constexpr index_t num_troe = 1;
    static constexpr index_t num_sri = 0;
    static constexpr index_t pressure_dependent_indexes[num_pressure_dependent] = {  17};
    static constexpr Parameters::reac_type pressure_dependent_reac_types[num_pressure_dependent] = {Parameters::reac_type::troe_reac};
    static constexpr index_t is_third_body_indexes[num_pressure_dependent] = {     2};
    
    static constexpr double reac_const[num_pressure_dependent][3] = {
        //               A_0                b_0                E_0 
        {           1.6e+34,             -5.49,            1987.0}     // 17. NH2+NH2(+M)=N2H4(+M)
    };
    
    static constexpr double troe[1][4] = {
        //              alfa               T***                 T*                T** 
        {              0.31,             1e-30,             1e+30,             1e+30}     // 17. NH2+NH2(+M)=N2H4(+M)
    };
    
    static constexpr double *sri = nullptr;
    
    static constexpr index_t num_plog = 4;
    static constexpr index_t num_plog_levels = 12;
    static constexpr index_t plog_indexes[num_plog] = {  20,  25,  29,  30};
    static constexpr index_t plog_seperators[num_plog_levels+1] = {   0,   3,   6,   9,  12};
    
    static constexpr double plog[num_plog_levels][4] = {
        //               P_1                A_1                b_1                E_1 
        {           10000.0,           2.3e+43,             -9.55,           64468.0},    // 20. N2H3=N2H2+H
        {          100000.0,           3.6e+47,            -10.38,           69009.0},    // 
        {         1000000.0,           1.8e+45,             -9.39,           70141.0},    // 
        {           10000.0,           5.6e+36,             -7.75,           70340.0},    // 25. N2H2=NNH+H
        {          100000.0,           1.8e+40,             -8.41,           73320.0},    // 
        {         1000000.0,           3.1e+41,             -8.42,           76102.0},    // 
        {           10000.0,           5.9e+32,             -6.99,           51791.0},    // 29. H2NN=NNH+H
        {          100000.0,           9.6e+35,             -7.57,           54841.0},    // 
        {         1000000.0,             5e+36,             -7.43,           57295.0},    // 
        {           10000.0,           7.2e+28,             -7.77,           50758.0},    // 30. H2NN=NNH+H
        {          100000.0,           3.2e+31,             -6.22,           52318.0},    // 
        {         1000000.0,           5.1e+33,             -6.52,           54215.0}     // 
    };
    
    
};    // struct chemkin_otomo2018_without_o_struct

#endif   // CHEMKIN_OTOMO2018_WITHOUT_O_H
