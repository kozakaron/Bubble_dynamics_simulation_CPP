function value = getfield_or_default(structure, field, default)
    % GETFIELD_OR_DEFAULT Returns the value of a field in a struct or a default value if the field does not exist.
    if isstruct(structure) && isfield(structure, field)
        value = structure.(field);
    else
        value = default;
    end
end
