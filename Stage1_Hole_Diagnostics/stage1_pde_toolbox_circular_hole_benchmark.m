function P = stage1_pde_toolbox_circular_hole_benchmark(C, varargin)
%STAGE1_PDE_TOOLBOX_CIRCULAR_HOLE_BENCHMARK
% PDE Toolbox benchmark inspired by MathWorks'
% StressConcentrationInPlateWithCircularHoleExample.
%
% This uses an exact circular hole in decsg, solves the same rectangular
% plate under remote y-tension, and samples tangential stress near the hole
% using StaticStructuralResults/interpolateStress.

    ip = inputParser;
    addParameter(ip, 'Hmax', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(ip, 'Hmin', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(ip, 'Hgrad', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(ip, 'Nphi', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 16));
    addParameter(ip, 'EpsShift', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
    parse(ip, varargin{:});

    must(C, 'A');
    must(C, 'B');
    must(C, 'E');
    must(C, 'nu');

    hole = get_single_circular_hole(C);

    A = C.A;
    B = C.B;
    c = hole.center(:).';
    R = hole.r;

    if getf(C, 'ps', 1) == 1
        analysisType = 'static-planestrain';
    else
        analysisType = 'static-planestress';
    end

    model = createpde('structural', analysisType);

    R1 = [3 4 0 A A 0 -B -B B B]';
    C1 = [1 c(1) c(2) R 0 0 0 0 0 0]';
    gd = [R1 C1];
    ns = char('R1','C1')';
    dl = decsg(gd, 'R1-C1', ns);
    geometryFromEdges(model, dl);

    structuralProperties(model, ...
        'YoungsModulus', C.E, ...
        'PoissonsRatio', C.nu);

    geom = model.Geometry;
    eTop = nearestEdge(geom, [0.5*A, B]);
    eBot = nearestEdge(geom, [0.5*A, -B]);

    sig0 = getf(getf(C, 'load', struct()), 'sig0', 1.0);
    structuralBoundaryLoad(model, 'Edge', eTop, 'SurfaceTraction', [0; sig0]);
    structuralBoundaryLoad(model, 'Edge', eBot, 'SurfaceTraction', [0; -sig0]);

    vLB = nearest_vertex(geom.Vertices, [0, -B]);
    vRB = nearest_vertex(geom.Vertices, [A, -B]);
    structuralBC(model, 'Vertex', vLB, 'XDisplacement', 0, 'YDisplacement', 0);
    structuralBC(model, 'Vertex', vRB, 'YDisplacement', 0);

    Hmax = ip.Results.Hmax;
    if isempty(Hmax)
        Hmax = getf(getf(C, 'mesh1', struct()), 'hmax', R/6);
    end
    Hmin = ip.Results.Hmin;
    if isempty(Hmin)
        Hmin = getf(getf(C, 'mesh1', struct()), 'hmin', []);
    end
    Hgrad = ip.Results.Hgrad;
    if isempty(Hgrad)
        Hgrad = getf(getf(C, 'mesh1', struct()), 'hgrad', []);
    end

    gmArgs = {'Hmax', Hmax};
    if ~isempty(Hmin)
        gmArgs = [gmArgs, {'Hmin', Hmin}]; %#ok<AGROW>
    end
    if ~isempty(Hgrad)
        gmArgs = [gmArgs, {'Hgrad', Hgrad}]; %#ok<AGROW>
    end
    generateMesh(model, gmArgs{:});

    result = solve(model);

    nphi = ip.Results.Nphi;
    if isempty(nphi)
        nphi = getf(getf(C, 'stage1', struct()), 'nphi', 720);
    end

    epsShift = ip.Results.EpsShift;
    if isempty(epsShift)
        mesh1 = getf(C, 'mesh1', struct());
        hhole = getf(mesh1, 'hhole', []);
        hmax = getf(mesh1, 'hmax', []);
        if ~isempty(hhole)
            epsShift = 0.25 * hhole;
        elseif ~isempty(hmax)
            epsShift = 0.10 * hmax;
        else
            epsShift = 1e-3 * R;
        end
        epsShift = min(epsShift, 0.10 * R);
        epsShift = max(epsShift, 1e-8 * max(R, 1));
    end

    phi = linspace(0, 2*pi, nphi+1).';
    phi(end) = [];

    n_mat = [cos(phi), sin(phi)];
    t_hat = [-sin(phi), cos(phi)];
    xb = c + R*n_mat;
    xq = xb + epsShift*n_mat;

    stressHole = interpolateStress(result, xq(:,1), xq(:,2));
    sig_xx = stressHole.sxx(:);
    sig_yy = stressHole.syy(:);
    sig_xy = stressHole.sxy(:);

    sig_tt = nan(nphi,1);
    for k = 1:nphi
        S = [sig_xx(k), sig_xy(k);
             sig_xy(k), sig_yy(k)];
        t = t_hat(k,:).';
        sig_tt(k) = t.' * S * t;
    end

    sig_pos = max(sig_tt, 0);
    [sigPeak, idxPeak] = max(sig_pos);

    P = struct();
    P.model = model;
    P.result = result;
    P.phi = phi;
    P.x = xb;
    P.xq = xq;
    P.sig_tt = sig_tt;
    P.sig_tt_peak = sigPeak;
    P.Kt = sigPeak / sig0;
    P.idx_peak = idxPeak;
    P.phi_peak = phi(idxPeak);
    P.phi_peak_deg = rad2deg(P.phi_peak);
    P.eps_shift = epsShift;
    P.mesh = model.Mesh;
    P.analysisType = analysisType;
    P.edgeIDs = struct('top', eTop, 'bottom', eBot);
    P.vertexIDs = struct('left_bottom', vLB, 'right_bottom', vRB);
end


function hole = get_single_circular_hole(C)
    if isfield(C, 'hole') && ~isempty(C.hole)
        hole = C.hole;
    elseif isfield(C, 'holes') && numel(C.holes) == 1
        hole = C.holes{1};
    else
        error('stage1_pde_toolbox_circular_hole_benchmark:HoleSpec', ...
            'Expected one circular hole in C.hole or C.holes.');
    end

    if ~strcmpi(strtrim(hole.type), 'circle')
        error('stage1_pde_toolbox_circular_hole_benchmark:UnsupportedHole', ...
            'Only circular holes are supported.');
    end
end


function idx = nearest_vertex(V, xq)
    if size(V,1) == 2
        P = V.';
    else
        P = V;
    end

    [~, idx] = min(sum((P - xq).^2, 2));
end


function must(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('stage1_pde_toolbox_circular_hole_benchmark:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end


function v = getf(S, field, default)
    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        v = S.(field);
    else
        v = default;
    end
end
