function data = run_simulation(cpar, json_path, executable_path, t_max, timeout, save_steps)
    % RUN_SIMULATION Runs a single simulation using the C++ executable.
    % Usage: data = run_simulation(cpar, json_path, executable_path, t_max, timeout, save_steps);

    if nargin < 2, json_path = 'ignore.json'; end
    if nargin < 3
        if ispc
            executable_path = './bin/main.exe'; % Windows
        else
            executable_path = './bin/main'; % Linux/Mac
        end
    end
    if nargin < 4, t_max = 1.0; end
    if nargin < 5, timeout = 60.0; end
    if nargin < 6, save_steps = true; end

    % Ensure species, fractions, and excitation_params are arrays
    if ~isfield(cpar, 'species') || isempty(cpar.species)
        cpar.species = {}; % Ensure it's an empty cell array of strings
    elseif ischar(cpar.species)
        cpar.species = {cpar.species}; % Convert single string to cell array
    end

    if ~isfield(cpar, 'fractions') || isempty(cpar.fractions)
        cpar.fractions = []; % Ensure it's an empty array of doubles
    elseif isscalar(cpar.fractions)
        cpar.fractions = {double(cpar.fractions)}; % Ensure it's a double array
    end

    if ~isfield(cpar, 'excitation_params') || isempty(cpar.excitation_params)
        cpar.excitation_params = []; % Ensure it's an empty array of doubles
    elseif isscalar(cpar.excitation_params)
        cpar.excitation_params = {double(cpar.excitation_params)}; % Ensure it's a double array
    end

    % Write the JSON file
    json_data = struct('cpar', cpar);
    json_str = jsonencode(json_data, PrettyPrint=true);
    fid = fopen(json_path, 'w');
    fwrite(fid, json_str, 'char');
    fclose(fid);

    % Check paths
    json_path = check_path(json_path);
    executable_path = check_path(executable_path);

    % Run the simulation
    command_list = {executable_path, '--run', json_path, '--tmax', num2str(t_max), '--timeout', num2str(timeout)};
    if save_steps
        command_list{end+1} = '--save';
    end
    
    % Create a process to run the command
    [status, cmdout] = system(strjoin(command_list, ' '));
    disp(cmdout);
    fprintf('\n%s returned with code %d.\n', executable_path, status);

    % Read the JSON-binary file
    data = read_json_and_binary(json_path);
end