function cpar = example_cpar()
    % EXAMPLE_CPAR Returns an example control parameter struct.
    % Usage: cpar = matlab_interface.example_cpar();

    % Define the file path
    file_name = fullfile('interface', 'example_cpar.json');

    % Read the JSON file
    if ~isfile(file_name)
        error('File "%s" does not exist.', file_name);
    end

    % Parse the JSON file
    json_data = jsondecode(fileread(file_name));

    % Extract the 'cpar' field
    if isfield(json_data, 'cpar')
        cpar = json_data.cpar;
    else
        cpar = struct(); % Return an empty struct if 'cpar' is not found
    end
end