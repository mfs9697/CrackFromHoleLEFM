function must_field(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('build_appended_hole_loop:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end