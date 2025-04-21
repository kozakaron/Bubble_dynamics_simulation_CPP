function path = check_path(path)
    % CHECK_PATH Checks if the path exists and returns the absolute path.
    
    if ~isfolder(path) && ~isfile(path)
        error('The path "%s" does not exist.', path);
    end
    path = fullfile(path);
end