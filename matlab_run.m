clc;clear;close all;
addpath('interface')
cpar=matlab_interface.example_cpar();
%
ps=matlab_interface.example_parameter_study();
ps.mechanism='chemkin_kaust2023_n2';
ps.target_specie='NH3';
ps.enable_evaporation=false;
% 
ps.species={'h2','n2'};
ps.fractions={0.75,0.25};
% 
ps.R_E.start=10.0e-6;
ps.R_E.end=600.0e-6;
ps.R_E.num_steps=591;
% 
ps.P_amb.value=10.0e5; %P_infty [Pa]

ps.excitation_params{2, 1}.value=1.0e4; %f
ps.excitation_params{1, 1}.start=-9.0e5; %p_A
ps.excitation_params{1, 1}.end=-30.0e5; %p_A
ps.excitation_params{1, 1}.num_steps=250; %p_A


matlab_interface.run_parameter_study(ps);
data=matlab_interface.read_parameter_study('_parameter_studies/test1/');

help matlab_interface
asd=matlab_interface.line_to_dict(data(100,:));