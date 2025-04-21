function parameter_study = example_parameter_study()
    % EXAMPLE_PARAMETER_STUDY Returns an example dictionary for a parameter study.
    % Usage: parameter_study = example_parameter_study();

    file_name = fullfile('interface', 'example_parameter_study.json');

    if ~isfile(file_name)
        error('File "%s" does not exist.', file_name);
    end

    % Parse the JSON file
    json_data = jsondecode(fileread(file_name));

    % Extract the 'cpar' field
    if isfield(json_data, 'parameter_study')
        parameter_study = json_data.parameter_study;
    else
        parameter_study = struct(); % Return an empty struct if 'parameter_study' is not found
    end
end