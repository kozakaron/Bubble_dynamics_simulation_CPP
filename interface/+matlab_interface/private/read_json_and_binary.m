function data = read_json_and_binary(file_path)
    % READ_JSON_AND_BINARY Reads a JSON-binary mix file and extracts data.
    % The binary format contains: t array, x array (2D state matrix), 
    % p_excitation array, and p_internal array.
    % 
    % The output structure includes:
    %   - t: time array
    %   - x: full state array (for backward compatibility)
    %   - R, R_dot, T, E_diss: individual state variables
    %   - p_excitation, p_internal: pressure time series
    %   - c_<species_name>: individual species concentrations
    %
    % Usage: data = read_json_and_binary(file_path);

    if ~isfile(file_path)
        error('File "%s" does not exist.', file_path);
    end

    fid = fopen(file_path, 'rb');
    content = fread(fid, '*uint8')';
    fclose(fid);

    binary_marker = uint8('<BINARY>');
    marker_index = strfind(content, binary_marker);

    if isempty(marker_index)
        error('The file does not contain the "<BINARY>" marker.');
    end

    % Split the JSON and binary parts
    json_part = char(content(1:marker_index-1));
    binary_part = content(marker_index + length(binary_marker):end);

    % Parse the JSON part
    data = jsondecode(json_part);

    % Extract metadata from the JSON
    if ~isfield(data, 'sol') || ~isfield(data.sol, 'num_saved_steps') || ~isfield(data.sol, 'num_dim')
        error('The JSON part is missing required fields.');
    end

    num_saved_steps = data.sol.num_saved_steps;
    num_dim = data.sol.num_dim;

    % Read the binary part sequentially:
    offset = 0;
    t = typecast(binary_part(offset + 1:offset + num_saved_steps * 8), 'double');
    offset = offset + num_saved_steps * 8;
    
    x = typecast(binary_part(offset + 1:offset + num_saved_steps * num_dim * 8), 'double');
    x = reshape(x, [num_dim, num_saved_steps])';
    offset = offset + num_saved_steps * num_dim * 8;
    
    p_excitation = typecast(binary_part(offset + 1:offset + num_saved_steps * 8), 'double');
    offset = offset + num_saved_steps * 8;
    
    p_internal = typecast(binary_part(offset + 1:offset + num_saved_steps * 8), 'double');

    % Add the binary data to the struct
    data.sol.t = t;
    data.sol.x = x;  % Keep for backward compatibility
    data.sol.p_excitation = p_excitation;
    data.sol.p_internal = p_internal;
    
    data.sol.R = x(:, 1);           % radius [m]
    data.sol.R_dot = x(:, 2);       % radius velocity [m/s]
    data.sol.T = x(:, 3);           % temperature [K]
    data.sol.E_diss = x(:, end);    % dissipated energy [J]
    
    if isfield(data, 'mechanism') && isfield(data.mechanism, 'species_names')
        species_names = data.mechanism.species_names;
        concentrations = x(:, 4:end-1);
        for i = 1:length(species_names)
            species_name = species_names{i};
            field_name = sprintf('c_%s', species_name);
            data.sol.(field_name) = concentrations(:, i);
        end
    end
end