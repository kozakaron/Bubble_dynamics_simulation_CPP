function data = run_simulation(cpar, options)
    % RUN_SIMULATION Runs a single simulation using the C++ executable.
    % Usage: data = run_simulation(cpar, t_max=1.0, timeout=60.0, save_steps=true, json_path="ignore.json");
    arguments
        cpar struct
        options.json_path (1,1) string = "ignore.json"
        options.executable_path (1,1) string = ""
        options.t_max (1,1) double = 1.0
        options.timeout (1,1) double = 60.0
        options.save_steps (1,1) logical = true
    end

    if options.executable_path == ""
        if ispc
            executable_path = './bin/main.exe'; % Windows
        else
            executable_path = './bin/main'; % Linux/Mac
        end
    else
        executable_path = options.executable_path;
    end

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
    fid = fopen(options.json_path, 'w');
    fwrite(fid, json_str, 'char');
    fclose(fid);

    % Check paths
    json_path = check_path(options.json_path);
    executable_path = check_path(executable_path);

    % Run the simulation
    command_list = {char(executable_path), '--run', char(json_path), '--tmax', num2str(options.t_max), '--timeout', num2str(options.timeout)};
    if options.save_steps
        command_list{end+1} = '--save';
    end
    
    % Create a process to run the command
    [status, cmdout] = system(strjoin(command_list, ' '));
    cmdout = regexprep(cmdout, '\x1B\[[0-9;]*[A-Za-z]', '');
    disp(cmdout);
    fprintf('\n%s returned with code %d.\n', executable_path, status);

    % Read the JSON-binary file
    data = read_json_and_binary(json_path);
end