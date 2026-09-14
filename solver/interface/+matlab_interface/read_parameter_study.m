function all_data = read_parameter_study(directory)
    % READ_PARAMETER_STUDY Reads a parameter study directory and returns a table.
    % Usage:
    %   all_data = read_parameter_study('path/to/parameter/study');
    %   good_data = all_data(all_data.success == "true", :);
    arguments
        directory (1,1) string  % Path to the parameter study directory
    end

    % Ensure the directory exists
    if ~isfolder(directory)
        error('The directory "%s" does not exist.', directory);
    end

    % Initialize variables
    fprintf('Found files:\n');
    all_data = table();
    num_files = 0;

    % Get all files in the directory and subdirectories
    files = dir(fullfile(directory, '**', '*.csv'));

    % Iterate through all files
    for i = 1:length(files)
        file_path = fullfile(files(i).folder, files(i).name);

        % Skip unwanted files (e.g., Jupyter checkpoints)
        if contains(files(i).folder, 'ipynb_checkpoints')
            continue;
        end

        % Read the CSV file
        num_files = num_files + 1;
        opts = detectImportOptions(file_path, 'Delimiter', ','); % Force comma as delimiter
        opts.VariableNamingRule = 'preserve'; % Preserve original column headers
        opts.ExtraColumnsRule = 'ignore'; % Ignore extra columns
        opts.EmptyLineRule = 'read'; % Handle empty lines
        opts = setvaropts(opts, 'QuoteRule', 'remove'); % Handle quoted strings
        current_data = readtable(file_path, opts);

        % Print file details
        relative_path = strrep(file_path, directory, '');
        fprintf('\t%-64s (%4d rows)\n', relative_path, height(current_data));
        
        % Convert NaN to empty string for parameter columns (these should be empty, not NaN)
        if any(strcmp(current_data.Properties.VariableNames, 'excitation_params'))
            current_data.excitation_params = fillmissing(current_data.excitation_params, 'constant', '');
        end
        if any(strcmp(current_data.Properties.VariableNames, 'liquid_eos_params'))
            current_data.liquid_eos_params = fillmissing(current_data.liquid_eos_params, 'constant', '');
        end

        % Append the data to the main table
        all_data = [all_data; current_data]; %#ok<AGROW>
    end

    % Print summary statistics
    fprintf('_______________________________________\n');
    total_rows = height(all_data);
    fprintf('Total: %4d rows in %2d files\n', total_rows, num_files);

    % Sort the data by 'energy_demand' column if it exists
    if any(strcmp(all_data.Properties.VariableNames, 'energy_demand'))
        all_data = sortrows(all_data, 'energy_demand', 'ascend');
    end
end