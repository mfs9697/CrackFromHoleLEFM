function C = local_scale_stage1_mesh(C0, meshScale)

    C = C0;

    if isfield(C, 'mesh1')
        if isfield(C.mesh1, 'hmin') && ~isempty(C.mesh1.hmin)
            C.mesh1.hmin = meshScale * C0.mesh1.hmin;
        end
        if isfield(C.mesh1, 'hmax') && ~isempty(C.mesh1.hmax)
            C.mesh1.hmax = meshScale * C0.mesh1.hmax;
        end
        if isfield(C.mesh1, 'hhole') && ~isempty(C.mesh1.hhole)
            C.mesh1.hhole = meshScale * C0.mesh1.hhole;
        end
    end
end