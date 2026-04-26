function G2 = geom_hole_shortcrack(C, I, theta, varargin)
%GEOM_HOLE_SHORTCRACK  Build geometry/mesh for plate + hole + short crack.
%
%   G2 = geom_hole_shortcrack(C, I, theta)
%   G2 = geom_hole_shortcrack(C, I, theta, 'a0', a0)
%
% Input:
%   C       configuration struct from cfg_hole_initiation
%   I       initiation struct from find_hole_initiation_point
%   theta   crack angle relative to inward normal, in radians
%
% Name-value pairs:
%   'a0'    short crack length (default = C.a0)
%
% Output:
%   G2      cracked geometry struct with fields:
%       .p             node coordinates [np x 2]
%       .t             element connectivity [nt x 3] or [nt x 6]
%       .geom          PDE model / geometry object if available
%       .meshobj       mesh object if available
%       .outerPoly     outer boundary polygon
%       .holeLoops     hole polygon loops
%       .hole          single circular hole spec
%       .crack         crack data struct
%       .tip           crack-tip metadata
%       .edgeSets      boundary node sets (if available)
%       .meta          auxiliary metadata
%
% Crack direction convention:
%   e_dir = cos(theta)*n_in_star + sin(theta)*t_hat_star
%
% so that:
%   theta = 0   -> inward normal direction
%   theta > 0   -> rotated toward positive tangent
%
% Notes:
%   This is a first implementation draft. It tries to use the CrackPath
%   "pencil-domain" backend if available. If the exact helper signature
%   differs in your local code, only the internal call block near
%   "TRY 1 / TRY 2" should need adjustment.

    % ------------------------------------------------------------
    % parse inputs
    % ------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'a0', getf(C, 'a0', []), @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    a0 = p.Results.a0;

    if isempty(a0)
        error('geom_hole_shortcrack:MissingA0', ...
            'Short crack length a0 must be provided in C.a0 or as a name-value pair.');
    end

    % ------------------------------------------------------------
    % checks
    % ------------------------------------------------------------
    must(C, 'A');
    must(C, 'B');
    must(I, 'x_star');
    must(I, 'n_mat_star');
    must(I, 't_hat_star');

    A = C.A;
    B = C.B;

    hole = get_single_circular_hole(C);

    % ------------------------------------------------------------
    % crack start, direction, tip
    % ------------------------------------------------------------
    x0   = I.x_star(:).';
    n_mat = normalize_row(I.n_mat_star(:).');
    tHat = normalize_row(I.t_hat_star(:).');

    e_dir = cos(theta) * n_mat + sin(theta) * tHat;
    e_dir = normalize_row(e_dir);

    xtip = x0 + a0 * e_dir;

    % basic sanity check: tip should lie inside the plate box
    if ~(xtip(1) > 0 && xtip(1) < A && xtip(2) > -B && xtip(2) < B)
        warning('geom_hole_shortcrack:TipOutsideBox', ...
            ['Short crack tip lies outside or very near the outer plate box. ', ...
             'theta=%.6f rad, xtip=(%.6g, %.6g).'], theta, xtip(1), xtip(2));
    end

    % crack polyline: just one straight segment in this first paper
    nseg = 2;
    s = linspace(0,1,nseg).';
    P = x0 + s.*(xtip - x0);

    % ------------------------------------------------------------
    % outer boundary and holes
    % ------------------------------------------------------------
    outerPoly = [ ...
        0, -B;
        A, -B;
        A,  B;
        0,  B ];

    holes = normalize_holes(C);
    holeLoops = holes_to_loops(holes);

    % ------------------------------------------------------------
    % mesh controls
    % ------------------------------------------------------------
    mesh2 = getf(C, 'mesh2', struct());

    hmax       = getf(mesh2, 'hmax', 0.010);
    hhole      = getf(mesh2, 'hhole', []);
    hcrack     = getf(mesh2, 'hcrack', []);
    hgrad      = getf(mesh2, 'hgrad', 1.25);
    chw        = getf(mesh2, 'chw', 0.0015);
    tip_radius = getf(mesh2, 'tip_radius', 0.0020);

    % plotting option
    plotStruct = getf(C, 'plot', struct());
    showMesh   = getf(plotStruct, 'show_mesh2', false);

    % ------------------------------------------------------------
    % initialize output fields
    % ------------------------------------------------------------
    G2 = struct();
    G2.p         = [];
    G2.t         = [];
    G2.geom      = [];
    G2.meshobj   = [];
    G2.outerPoly = outerPoly;
    G2.holeLoops = holeLoops;
    G2.hole      = hole;
    G2.edgeSets  = struct();

    % ------------------------------------------------------------
    % crack metadata
    % ------------------------------------------------------------
    crack = struct();
    crack.x0       = x0;
    crack.xtip     = xtip;
    crack.a0       = a0;
    crack.theta    = theta;
    crack.thetaDeg = theta * 180 / pi;
    crack.e_dir    = e_dir;
    crack.n_mat    = n_mat;

    crack.t_hat    = tHat;
    crack.polyline = P;

    G2.crack = crack;

    % tip frame:
    % tangent = crack direction
    % normal  = +90 degree rotation of tangent
    t_tip = e_dir;
    n_tip = [-t_tip(2), t_tip(1)];

    tip = struct();
    tip.x          = xtip;
    tip.tangent    = t_tip;
    tip.normal     = n_tip;
    tip.radiusJ    = tip_radius;

    G2.tip = tip;

    % ------------------------------------------------------------
    % TRY 1: use pencil-domain crack meshing backend
    % ------------------------------------------------------------
    ok_mesh = false;
    lastErr = [];

    try
        % This block assumes your CrackPath helper can build a domain with:
        %   outer polygon + hole loops + crack polyline + channel width + sizes
        %
        % Depending on the exact helper signature in your local project,
        % you may only need to edit this call.
        %
        % Suggested internal argument package:
        S = struct();
        S.outerPoly = outerPoly;
        S.holeLoops = holeLoops;
        S.P         = P;
        S.chw       = chw;
        S.hmax      = hmax;
        S.hhole     = hhole;
        S.hcrack    = hcrack;
        S.hgrad     = hgrad;

        [dom, meta_dom] = call_build_domain_pencil_polyline(S);
        [p, t, meshobj, geomobj, edgeSets] = call_mesh_pencil_domain(dom, S);

        G2.p       = p;
        G2.t       = t;
        G2.meshobj = meshobj;
        G2.geom    = geomobj;
        G2.edgeSets = edgeSets;

        G2.meta = struct();
        G2.meta.backend   = 'pencil_polyline';
        G2.meta.domain    = dom;
        G2.meta.domainAux = meta_dom;
        G2.meta.mesh2     = mesh2;

        ok_mesh = true;

    catch ME1
        lastErr = ME1;
    end

    % ------------------------------------------------------------
    % TRY 2: fallback diagnostic geometry only
    % ------------------------------------------------------------
    if ~ok_mesh
        try
            [mdl, dl, bt, gd, ns, sf] = build_pde_geometry_with_holes(outerPoly, holeLoops);

            % We store the geometry, but a plain PDE domain cannot represent
            % an internal open crack directly. So this fallback is diagnostic
            % only, not the final cracked mesh.
            G2.geom = mdl;

            G2.meta = struct();
            G2.meta.backend = 'fallback_geometry_only';
            G2.meta.decsg = struct('dl', dl, 'bt', bt, 'gd', gd, 'ns', ns, 'sf', sf);
            G2.meta.mesh2 = mesh2;
            G2.meta.fallbackNote = [ ...
                'Hole geometry built successfully, but no cracked mesh was created. ', ...
                'Please connect the local CrackPath pencil meshing helper signature.'];

        catch ME2
            if isempty(lastErr)
                rethrow(ME2);
            else
                error('geom_hole_shortcrack:BackendFailure', ...
                    ['Failed to build cracked geometry.\n\n', ...
                     'Primary pencil-domain backend error:\n%s\n\n', ...
                     'Fallback geometry error:\n%s'], ...
                     lastErr.message, ME2.message);
            end
        end
    end

    % ------------------------------------------------------------
    % edge sets fallback
    % ------------------------------------------------------------
    if isempty(fieldnames(G2.edgeSets)) && ~isempty(G2.p)
        G2.edgeSets = build_edge_sets_basic(G2.p, A, B, holes);
    end

    % ------------------------------------------------------------
    % optional plot
    % ------------------------------------------------------------
    if showMesh
        figure; clf; hold on; axis equal; box on

        if ~isempty(G2.t) && ~isempty(G2.p)
            if size(G2.t,2) >= 3
                triplot(G2.t(:,1:3), G2.p(:,1), G2.p(:,2), ...
                    'Color', [0.75 0.75 0.75]);
            end
        end

        % outer boundary
        plot([outerPoly(:,1); outerPoly(1,1)], ...
             [outerPoly(:,2); outerPoly(1,2)], ...
             'k-', 'LineWidth', 1.2);

        % hole boundaries
        for ih = 1:numel(holeLoops)
            H = holeLoops{ih};
            plot([H(:,1); H(1,1)], [H(:,2); H(1,2)], ...
                'r-', 'LineWidth', 1.2);
        end

        % crack
        plot(P(:,1), P(:,2), 'b-', 'LineWidth', 2.0);
        plot(x0(1),   x0(2),   'ko', 'MarkerSize', 6, 'LineWidth', 1.2);
        plot(xtip(1), xtip(2), 'bo', 'MarkerSize', 6, 'LineWidth', 1.2);

        % local frame at start point
        scl = max(a0, 1e-6);
        quiver(x0(1), x0(2), scl*n_mat(1), scl*n_mat(2), 0, ...
            'Color', [0 0.5 0], 'LineWidth', 1.4, 'MaxHeadSize', 0.8);
        quiver(x0(1), x0(2), scl*tHat(1), scl*tHat(2), 0, ...
            'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, 'MaxHeadSize', 0.8);

        title(sprintf('Stage II geometry: hole + short crack, \\theta = %.2f^\\circ', ...
            crack.thetaDeg));
        xlabel('x');
        ylabel('y');
        xlim([0, A]);
        ylim([-B, B]);
    end

end


% =========================================================================
% backend adapters
% =========================================================================

function [dom, meta] = call_build_domain_pencil_polyline(S)
%CALL_BUILD_DOMAIN_PENCIL_POLYLINE  Thin adapter around local CrackPath helper.
%
% Edit this function only if your local helper signature differs.

    if exist('build_domain_pencil_polyline', 'file') ~= 2
        error('geom_hole_shortcrack:MissingHelper', ...
            'build_domain_pencil_polyline.m is not on the MATLAB path.');
    end

    % --- attempt 1: pass one struct ---
    try
        out = build_domain_pencil_polyline(S);
        if isstruct(out)
            dom = out;
            meta = struct();
            return;
        end
    catch
    end

    % --- attempt 2: canonical expanded signature guess ---
    try
        [dom, meta] = build_domain_pencil_polyline( ...
            S.outerPoly, S.P, ...
            'holes', S.holeLoops, ...
            'chw', S.chw, ...
            'hmax', S.hmax, ...
            'hgrad', S.hgrad, ...
            'hhole', S.hhole, ...
            'hcrack', S.hcrack);
        return;
    catch ME
        error('geom_hole_shortcrack:BuildDomainCall', ...
            ['Could not call build_domain_pencil_polyline with the guessed ', ...
             'signatures. Please adapt call_build_domain_pencil_polyline().\n%s'], ...
             ME.message);
    end
end


function [p, t, meshobj, geomobj, edgeSets] = call_mesh_pencil_domain(dom, S)
%CALL_MESH_PENCIL_DOMAIN  Thin adapter around local CrackPath helper.
%
% Edit this function only if your local helper signature differs.

    if exist('mesh_pencil_domain', 'file') ~= 2
        error('geom_hole_shortcrack:MissingHelper', ...
            'mesh_pencil_domain.m is not on the MATLAB path.');
    end

    meshobj = [];
    geomobj = [];
    edgeSets = struct();

    % --- attempt 1: one input / one struct output ---
    try
        out = mesh_pencil_domain(dom);
        if isstruct(out) && isfield(out, 'p') && isfield(out, 't')
            p = out.p;
            t = out.t;
            if isfield(out, 'meshobj'),  meshobj = out.meshobj; end
            if isfield(out, 'geom'),     geomobj = out.geom;    end
            if isfield(out, 'edgeSets'), edgeSets = out.edgeSets; end
            return;
        end
    catch
    end

    % --- attempt 2: guessed expanded return signature ---
    try
        [p, t, meshobj, geomobj, edgeSets] = mesh_pencil_domain( ...
            dom, ...
            'hmax', S.hmax, ...
            'hgrad', S.hgrad, ...
            'hhole', S.hhole, ...
            'hcrack', S.hcrack);
        return;
    catch ME
        error('geom_hole_shortcrack:MeshDomainCall', ...
            ['Could not call mesh_pencil_domain with the guessed signatures. ', ...
             'Please adapt call_mesh_pencil_domain().\n%s'], ...
             ME.message);
    end
end


% =========================================================================
% helpers
% =========================================================================

function holes = normalize_holes(C)
%NORMALIZE_HOLES Normalize hole input to a cell array of structs.

    if isfield(C, 'holes') && ~isempty(C.holes)
        holes = C.holes;
    elseif isfield(C, 'hole') && ~isempty(C.hole)
        holes = {C.hole};
    else
        holes = {};
    end

    if isstruct(holes)
        holes = num2cell(holes);
    end

    if ~iscell(holes)
        error('geom_hole_shortcrack:BadHoleInput', ...
            'C.holes must be a cell array or struct array; or provide C.hole.');
    end
end


function hole = get_single_circular_hole(C)
%GET_SINGLE_CIRCULAR_HOLE Return the single circular hole spec.

    holes = normalize_holes(C);

    if numel(holes) ~= 1
        error('geom_hole_shortcrack:HoleSpec', ...
            'This first draft supports exactly one circular hole.');
    end

    hole = holes{1};

    if ~isstruct(hole) || ~isfield(hole, 'type')
        error('geom_hole_shortcrack:BadHoleSpec', ...
            'Hole specification must be a struct with field "type".');
    end

    if ~strcmpi(strtrim(hole.type), 'circle')
        error('geom_hole_shortcrack:UnsupportedHoleType', ...
            'Only circular holes are supported in this first draft.');
    end

    must(hole, 'center');
    must(hole, 'r');
end

function edgeSets = build_edge_sets_basic(p, A, B, holes)
%BUILD_EDGE_SETS_BASIC Basic geometric boundary-node sets.
% Fallback only.

    x = p(:,1);
    y = p(:,2);

    Lx = max(A, 1.0);
    Ly = max(2*B, 1.0);
    tol = 1e-8 * max(Lx, Ly);

    edgeSets = struct();
    edgeSets.left   = find(abs(x - 0) < tol);
    edgeSets.right  = find(abs(x - A) < tol);
    edgeSets.bottom = find(abs(y + B) < tol);
    edgeSets.top    = find(abs(y - B) < tol);

    hole_ids = false(size(p,1),1);

    for ih = 1:numel(holes)
        hk = holes{ih};
        if strcmpi(strtrim(hk.type), 'circle')
            c = hk.center(:).';
            r = hk.r;
            rr = hypot(x - c(1), y - c(2));
            hole_ids = hole_ids | (abs(rr - r) < 5*tol);
        end
    end

    edgeSets.hole = find(hole_ids);
end


function must(S, field)
%MUST Error if field does not exist or is empty.

    if ~isfield(S, field) || isempty(S.(field))
        error('geom_hole_shortcrack:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end


function v = getf(S, field, default)
%GETF Get struct field or default.

    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        v = S.(field);
    else
        v = default;
    end
end