function text = print_data(data, print_it)
    % PRINT_DATA Prints the data struct in an organized way.
    % Arguments:
    %   data: struct containing the simulation data.
    %   print_it: if true, prints the text; if false, returns the text as a string.

    if nargin < 2
        print_it = true; % Default to printing the text
    end

    % Extract control parameters and solution
    cpar = getfield_or_default(data, 'cpar', struct());
    sol = getfield_or_default(data, 'sol', struct());

    % Control Parameters
    text = sprintf("Control Parameters:\n");
    text = text + sprintf("  ID: %s\n", getfield_or_default(cpar, 'ID', 'N/A'));
    text = text + sprintf("  Mechanism: %s\n", getfield_or_default(cpar, 'mechanism', 'N/A'));
    text = text + sprintf("  R_E: %.2f [um]\n", 1e6 * getfield_or_default(cpar, 'R_E', NaN));
    text = text + sprintf("  Species: %s\n", strjoin(getfield_or_default(cpar, 'species', {}), ', '));
    text = text + sprintf("  Fractions: %s\n", mat2str(getfield_or_default(cpar, 'fractions', [])));
    text = text + sprintf("  P_amb: %.2f [Pa]\n", getfield_or_default(cpar, 'P_amb', NaN));
    text = text + sprintf("  T_inf: %.2f [K]\n", getfield_or_default(cpar, 'T_inf', NaN));
    text = text + sprintf("  alfa_M: %.2f [-]\n", getfield_or_default(cpar, 'alfa_M', NaN));
    text = text + sprintf("  P_v: %.2f [Pa]\n", getfield_or_default(cpar, 'P_v', NaN));
    text = text + sprintf("  mu_L: %.4f [Pa·s]\n", getfield_or_default(cpar, 'mu_L', NaN));
    text = text + sprintf("  rho_L: %.2f [kg/m³]\n", getfield_or_default(cpar, 'rho_L', NaN));
    text = text + sprintf("  c_L: %.2f [m/s]\n", getfield_or_default(cpar, 'c_L', NaN));
    text = text + sprintf("  Surfactant: %.2f\n", getfield_or_default(cpar, 'surfactant', NaN));
    text = text + sprintf("  Enable Heat Transfer: %s\n", bool_to_string(getfield_or_default(cpar, 'enable_heat_transfer', false)));
    text = text + sprintf("  Enable Evaporation: %s\n", bool_to_string(getfield_or_default(cpar, 'enable_evaporation', false)));
    text = text + sprintf("  Enable Reactions: %s\n", bool_to_string(getfield_or_default(cpar, 'enable_reactions', false)));
    text = text + sprintf("  Enable Dissipated Energy: %s\n", bool_to_string(getfield_or_default(cpar, 'enable_dissipated_energy', false)));
    text = text + sprintf("  Target Specie: %s\n", getfield_or_default(cpar, 'target_specie', 'N/A'));
    text = text + sprintf("  Excitation Params: %s\n", mat2str(getfield_or_default(cpar, 'excitation_params', [])));
    text = text + sprintf("  Excitation Type: %s\n", getfield_or_default(cpar, 'excitation_type', 'N/A'));

    % Simulation Info
    text = text + sprintf("\nSimulation Info:\n");
    text = text + sprintf("  Success: %s\n", bool_to_string(getfield_or_default(sol, 'success', false)));
    text = text + sprintf("  Error: %s\n", getfield_or_default(sol, 'error', 'N/A'));
    text = text + sprintf("  Runtime: %.2f [s]\n", getfield_or_default(sol, 'runtime', NaN));
    text = text + sprintf("  Num Steps: %d\n", getfield_or_default(sol, 'num_steps', NaN));
    text = text + sprintf("  Num Repeats: %d\n", getfield_or_default(sol, 'num_repeats', NaN));
    text = text + sprintf("  Num Function Evaluations: %d\n", getfield_or_default(sol, 'num_fun_evals', NaN));
    text = text + sprintf("  Num Jacobian Evaluations: %d\n", getfield_or_default(sol, 'num_jac_evals', NaN));

    % Handle t_last
    t_last = getfield_or_default(sol, 't', NaN);
    if isnumeric(t_last) && all(~isnan(t_last)) % Ensure t_last is numeric and all elements are not NaN
        t_last = t_last(end); % Take the last element
    else
        t_last = NaN; % Default to NaN if invalid
    end
    text = text + sprintf("  t_last = %.6e [s]\n", t_last);

    % Results
    text = text + sprintf("\nResults:\n");
    text = text + sprintf("  Dissipated Energy: %.6e [J]\n", getfield_or_default(data, 'dissipated_energy', NaN));
    text = text + sprintf("  n_target_specie: %.6e [mol]\n", getfield_or_default(data, 'n_target_specie', NaN));
    text = text + sprintf("  Energy Demand: %.2f [MJ/kg]\n", getfield_or_default(data, 'energy_demand', NaN));

    % Print or return the text
    if print_it
        fprintf('%s', text);
    end
end

function value = getfield_or_default(structure, field, default)
    % GETFIELD_OR_DEFAULT Returns the value of a field in a struct or a default value if the field does not exist.
    if isfield(structure, field)
        value = structure.(field);
    else
        value = default;
    end
end

function str = bool_to_string(value)
    % BOOL_TO_STRING Converts a boolean value to 'true' or 'false' string.
    if value
        str = 'true';
    else
        str = 'false';
    end
end