function run_parameter_study(parameter_study, json_path, executable_path, t_max, timeout, save_directory)
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
    if nargin < 6, save_directory = './_parameter_studies/test'; end

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
    fid = fopen(json_path, 'w');
    if fid == -1
        error('Could not open file "%s" for writing.', json_path);
    end
    fwrite(fid, json_str, 'char');
    fclose(fid);

    % Check paths
    json_path = check_path(json_path);
    executable_path = check_path(executable_path);

    % Run simulation (call the executable)
    start_time = tic;
    command_list = {
        executable_path, '--parameter_study', json_path, ...
        '--tmax', num2str(t_max), ...
        '--timeout', num2str(timeout), ...
        '--directory', save_directory
    };

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