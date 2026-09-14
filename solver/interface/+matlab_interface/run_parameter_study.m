function run_parameter_study(parameter_study, options)
    % RUN_PARAMETER_STUDY Runs a parameter study using the C++ executable.
    % Saves the provided parameter_study dictionary to a JSON file and uses it as input to run the study.
    % Captures stdout and stderr in real-time and prints them as they arrive.
    %
    % Arguments:
    %   parameter_study: Struct containing the parameter study configuration.
    %   json_path: Path to the JSON file (default: 'ignore.json').
    %   executable_path: Path to the executable (default: './bin/main' or './bin/main.exe').
    %   t_max: Maximum simulation time (default: 1.0).
    %   timeout: Timeout for the simulation (default: 60.0).
    %   save_directory: Directory to save the results (default: './_parameter_studies/test').
    %   cpu_count: Number of CPU cores to use (0 = auto, default: 0).
    arguments
        parameter_study struct
        options.json_path (1,1) string = "ignore.json"
        options.executable_path (1,1) string = ""
        options.t_max (1,1) double = 1.0
        options.timeout (1,1) double = 60.0
        options.save_directory (1,1) string = "./_parameter_studies/test"
        options.cpu_count (1,1) double {mustBeInteger, mustBeNonnegative} = 0
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
    if ~isfield(parameter_study, 'species') || isempty(parameter_study.species)
        parameter_study.species = {}; % Ensure it's an empty cell array of strings
    elseif ischar(parameter_study.species)
        parameter_study.species = {parameter_study.species}; % Convert single string to cell array
    end

    if ~isfield(parameter_study, 'fractions') || isempty(parameter_study.fractions)
        parameter_study.fractions = []; % Ensure it's an empty array of doubles
    elseif isscalar(parameter_study.fractions)
        parameter_study.fractions = {double(parameter_study.fractions)}; % Ensure it's a double array
    end

    if ~isfield(parameter_study, 'excitation_params') || isempty(parameter_study.excitation_params)
        parameter_study.excitation_params = {}; % Ensure it's an empty cell array
    elseif isstruct(parameter_study.excitation_params)
        % If it's a single struct, wrap it in a cell array
        parameter_study.excitation_params = {parameter_study.excitation_params};
    elseif iscell(parameter_study.excitation_params)
        % Validate that all elements are structs
        for i = 1:length(parameter_study.excitation_params)
            if ~isstruct(parameter_study.excitation_params{i})
                error('All elements in excitation_params must be structs.');
            end
        end
    else
        error('excitation_params must be a struct or a cell array of structs.');
    end

    % Write the JSON file
    json_data = struct('parameter_study', parameter_study);
    json_str = jsonencode(json_data, 'PrettyPrint', true);
    fid = fopen(options.json_path, 'w');
    if fid == -1
        error('Could not open file "%s" for writing.', options.json_path);
    end
    fwrite(fid, json_str, 'char');
    fclose(fid);

    % Check paths
    json_path = check_path(options.json_path);
    executable_path = check_path(executable_path);

    % Run simulation (call the executable)
    start_time = tic;
    command_list = {
        char(executable_path), '--parameter_study', char(json_path), ...
        '--tmax', num2str(options.t_max), ...
        '--timeout', num2str(options.timeout), ...
        '--directory', char(options.save_directory)
    };
    if options.cpu_count > 0
        command_list{end + 1} = '--cpu';
        command_list{end + 1} = num2str(options.cpu_count);
    end

    % Create a process to run the command and capture stdout in real-time
    process = System.Diagnostics.Process();
    process.StartInfo.FileName = executable_path;
    process.StartInfo.Arguments = strjoin(command_list(2:end), ' ');
    process.StartInfo.UseShellExecute = true;
    process.StartInfo.CreateNoWindow = false;

    % Start the process
    process.Start();

    % Wait for the process to exit
    process.WaitForExit();
    elapsed_time = toc(start_time);

    % Get the exit code
    status = process.ExitCode;
    fprintf('\n%s returned with code %d after %.4f seconds.\n', executable_path, status, elapsed_time);
end