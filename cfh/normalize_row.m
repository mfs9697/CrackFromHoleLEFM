function v = normalize_row(v)
    v = v(:).';
    nv = norm(v);
    if nv <= 0
        error('build_domain_hole_pencil_polyline:ZeroVector', ...
            'Cannot normalize a zero vector.');
    end
    v = v / nv;
end
