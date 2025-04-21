function data = read_json_and_binary(file_path)
    % READ_JSON_AND_BINARY Reads a JSON-binary mix file and extracts data.
    % Usage: data = _read_json_and_binary(file_path);

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

    % Read the binary part
    t = typecast(binary_part(1:num_saved_steps * 8), 'double');
    x = typecast(binary_part(num_saved_steps * 8 + 1:end), 'double');
    x = reshape(x, [num_dim, num_saved_steps])';

    % Add the binary data to the struct
    data.sol.t = t;
    data.sol.x = x;
end