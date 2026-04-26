function M = mesh_hole_pencil_domain(D, varargin)
%MESH_HOLE_PENCIL_DOMAIN  Mesh the appended-hole geometry.
%
%   M = mesh_hole_pencil_domain(D)
%   M = mesh_hole_pencil_domain(D, 'Hmax', 0.005, 'Hgrad', 1.2)
%
% Input:
%   D   geometry-description struct from build_domain_hole_pencil_polyline
%       Required fields:
%         .outerPoly    [No x 2] outer boundary polygon (CCW)
%         .holeLoops    cell array of inner polygons
%
% Name-value options:
%   'Hmax'      : global max element size (default = [])
%   'Hgrad'     : mesh growth factor      (default = [])
%   'Hmin'      : optional min element size (default = [])
%   'GeometricOrder' : 'linear' (default) or 'quadratic'
%   'PlotGeom'  : logical, plot PDE geometry labels (default false)
%   'PlotMesh'  : logical, plot mesh (default false)
%   'Hvertex'   : optional manual PDE-geometry vertex refinement
%   'Hedge'     : optional manual PDE-geometry edge refinement
%
% Output:
%   M   struct with fields:
%       .p           node coordinates [np x 2]
%       .t           element connectivity [nt x 3] or [nt x 6]
%       .geom        PDE model
%       .meshobj     PDE mesh object
%       .dl          decsg decomposed geometry
%       .bt          boolean table from decsg
%       .gd          geometry description matrix
%       .ns          names matrix
%       .sf          set formula
%       .edgeSets    boundary node sets
%       .region      region metadata
%
% Notes:
%   - This version is merged-only: there is no separate channel geometry.
%   - If manual Hvertex/Hedge is supplied, automatic tip refinement is skipped.
%   - If manual Hvertex/Hedge is not supplied, the function first tries
%     pure-geometry sharp-pencil identification. If that fails, it falls back
%     to a temporary background mesh to recover the geometry IDs, then remeshes
%     with local tip refinement.

    % ------------------------------------------------------------
    % parse options
    % ------------------------------------------------------------
    ip = inputParser;
    addParameter(ip, 'Hmax', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(ip, 'Hgrad', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(ip, 'Hmin', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(ip, 'GeometricOrder', 'linear', @(s) ischar(s) || isstring(s));
    addParameter(ip, 'PlotGeom', false, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'PlotMesh', false, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'Hvertex', [], @(x) isempty(x) || iscell(x));
    addParameter(ip, 'Hedge',   [], @(x) isempty(x) || iscell(x));
    parse(ip, varargin{:});

    Hmax = ip.Results.Hmax;
    Hgrad = ip.Results.Hgrad;
    Hmin = ip.Results.Hmin;
    geometricOrder = char(ip.Results.GeometricOrder);
    plotGeom = ip.Results.PlotGeom;
    plotMesh = ip.Results.PlotMesh;
    Hvertex = ip.Results.Hvertex;
    Hedge   = ip.Results.Hedge;

    % ------------------------------------------------------------
    % checks
    % ------------------------------------------------------------
    must(D, 'outerPoly');
    must(D, 'holeLoops');

    outerPoly = D.outerPoly;
    holeLoops = D.holeLoops;

    validate_polygon(outerPoly, 'outerPoly');
    for k = 1:numel(holeLoops)
        validate_polygon(holeLoops{k}, sprintf('holeLoops{%d}', k));
    end

    % enforce orientations
    if signed_polygon_area(outerPoly) < 0
        outerPoly = flipud(outerPoly);
    end

    for k = 1:numel(holeLoops)
        if signed_polygon_area(holeLoops{k}) > 0
            holeLoops{k} = flipud(holeLoops{k});   % inner loops clockwise
        end
    end

    % ------------------------------------------------------------
    % build decomposed geometry: rectangle minus hole loops
    % ------------------------------------------------------------
    polys = [{outerPoly}, holeLoops(:).'];
    names = cell(1, numel(polys));
    names{1} = 'R1';
    for k = 2:numel(polys)
        names{k} = sprintf('H%d', k-1);
    end

    sf = 'R1';
    for k = 2:numel(names)
        sf = [sf, '-', names{k}]; %#ok<AGROW>
    end

    [gd, ns] = polygons_to_gd_ns(polys, names);
    [dl, bt] = decsg(gd, sf, ns);

    mdl = createpde();
    geometryFromEdges(mdl, dl);

    if plotGeom
        figure('Name', 'mesh_hole_pencil_domain: geometry', 'Color', 'w'); clf
        pdegplot(mdl, 'EdgeLabels', 'on', 'FaceLabels', 'on', 'VertexLabels', 'on');
        axis equal
        title('PDE geometry: rectangle - appended-hole loops')
    end

    % ------------------------------------------------------------
    % base mesh arguments (without local tip refinement)
    % ------------------------------------------------------------
    gmArgsBase = {'GeometricOrder', geometricOrder};

    if ~isempty(Hmax)
        gmArgsBase = [gmArgsBase, {'Hmax', Hmax}];
    end
    if ~isempty(Hgrad)
        gmArgsBase = [gmArgsBase, {'Hgrad', Hgrad}];
    end
    if ~isempty(Hmin)
        gmArgsBase = [gmArgsBase, {'Hmin', Hmin}];
    end

    % ------------------------------------------------------------
    % geometry IDs for the sharp pencil
    % ------------------------------------------------------------
    geomIDs = [];
    useManualRefinement = ~isempty(Hvertex) || ~isempty(Hedge);

    if ~useManualRefinement
        % First try the pure-geometry route
        try
            geomIDs = identify_sharp_pencil_geom_ids(mdl, D, 'Verbose', true);
        catch MEgeom
            fprintf(['mesh_hole_pencil_domain: geometry-only sharp-pencil identification failed:\n', ...
            '  %s\n', ...
            'Falling back to a temporary background mesh to recover geometry IDs.\n'], ...
                MEgeom.message);

            % Fallback: temporary background mesh, then identify IDs from it
            try
                msh0 = generateMesh(mdl, gmArgsBase{:});
                geomIDs = identify_pencil_ids_from_temp_mesh_local(mdl, msh0, D, 'Verbose', true);
            catch MEtmp
                warning('mesh_hole_pencil_domain:GeomIDFallbackFailed', ...
                    ['Temporary-mesh geometry-ID recovery also failed:\n  %s\n', ...
                     'Proceeding without automatic local tip refinement.'], ...
                     MEtmp.message);
                geomIDs = [];
            end
        end
    end

    % ------------------------------------------------------------
    % final mesh arguments
    % ------------------------------------------------------------
    gmArgs = gmArgsBase;

    % Manual local refinement overrides automatic
    if ~isempty(Hedge)
        gmArgs = [gmArgs, {'Hedge', Hedge}];
    elseif ~isempty(geomIDs)
        [Hface, ~] = choose_auto_refinement_sizes(Hmin, Hmax);
        gmArgs = [gmArgs, {'Hedge', {geomIDs.e_tip, Hface}}];
    end

    if ~isempty(Hvertex)
        gmArgs = [gmArgs, {'Hvertex', Hvertex}];
    elseif ~isempty(geomIDs)
        [~, Htip] = choose_auto_refinement_sizes(Hmin, Hmax);
        gmArgs = [gmArgs, {'Hvertex', {[geomIDs.v_tip], Htip}}];
    end

    % ------------------------------------------------------------
    % generate final mesh
    % ------------------------------------------------------------
    msh = generateMesh(mdl, gmArgs{:});

    p = msh.Nodes.';      % [np x 2]
    t = msh.Elements.';   % [nt x 3] or [nt x 6]

    % If geomIDs were not available earlier, try once more on the final mesh
    if isempty(geomIDs)
        try
            geomIDs = identify_pencil_ids_from_temp_mesh_local(mdl, msh, D, 'Verbose', false);
        catch
            geomIDs = [];
        end
    end

    % ------------------------------------------------------------
    % basic boundary node sets
    % ------------------------------------------------------------
    edgeSets = build_edge_sets_from_geometry(p, outerPoly, holeLoops, D);

    % ------------------------------------------------------------
    % region metadata
    % ------------------------------------------------------------
    region = struct();
    region.outerPoly = outerPoly;
    region.holeLoops = holeLoops;
    region.names     = names;
    region.nHoles    = numel(holeLoops);
    region.mode      = 'merged_appended_hole';
    region.geomIDs   = geomIDs;

    % ------------------------------------------------------------
    % optional mesh plot
    % ------------------------------------------------------------
    if plotMesh
        figure('Name', 'mesh_hole_pencil_domain: mesh', 'Color', 'w'); clf
        triplot(t(:,1:3), p(:,1), p(:,2), 'Color', [0.75 0.75 0.75]); hold on
        axis equal; box on

        plot_closed_polygon(outerPoly, 'k-', 1.2);
        for k = 1:numel(holeLoops)
            plot_closed_polygon(holeLoops{k}, 'r-', 1.2);
        end

        title('Mesh of rectangle - appended-hole loops')
        xlabel('x'); ylabel('y')
    end

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    M = struct();
    M.p       = p;
    M.t       = t;
    M.geom    = mdl;
    M.meshobj = msh;

    M.dl = dl;
    M.bt = bt;
    M.gd = gd;
    M.ns = ns;
    M.sf = sf;

    M.edgeSets = edgeSets;
    M.region   = region;
end


% =========================================================================
% helpers
% =========================================================================

function [gd, ns] = polygons_to_gd_ns(polys, names)
    nObj = numel(polys);
    if numel(names) ~= nObj
        error('mesh_hole_pencil_domain:NameCount', ...
            'names and polys must have the same length.');
    end

    gcols = cell(1, nObj);
    maxLen = 0;

    for j = 1:nObj
        P = polys{j};
        validate_polygon(P, sprintf('polys{%d}', j));

        if size(P,1) >= 2 && norm(P(end,:) - P(1,:), inf) < 1e-14
            P = P(1:end-1,:);
        end

        Nv = size(P,1);
        col = [2; Nv; P(:,1); P(:,2)];
        gcols{j} = col;
        maxLen = max(maxLen, numel(col));
    end

    gd = zeros(maxLen, nObj);
    for j = 1:nObj
        col = gcols{j};
        gd(1:numel(col), j) = col;
    end

    ns = char(names);
    ns = ns.';
end


function [Hface, Htip] = choose_auto_refinement_sizes(Hmin, Hmax)
    if ~isempty(Hmin)
        Hface = Hmin;
        Htip  = 0.5 * Hmin;
    elseif ~isempty(Hmax)
        Hface = 0.25 * Hmax;
        Htip  = 0.125 * Hmax;
    else
        error('mesh_hole_pencil_domain:NeedMeshScale', ...
            'Automatic tip refinement requires Hmin or Hmax.');
    end
end


function edgeSets = build_edge_sets_from_geometry(p, outerPoly, holeLoops, D)
    edgeSets = struct();

    if isfield(D, 'A') && isfield(D, 'B')
        A = D.A;
        B = D.B;
    else
        bb = [min(outerPoly(:,1)), max(outerPoly(:,1)), ...
              min(outerPoly(:,2)), max(outerPoly(:,2))];
        A = bb(2);
        B = max(abs([bb(3), bb(4)]));
    end

    x = p(:,1);
    y = p(:,2);
    tol = 1e-8 * max([A, 2*B, 1]);

    edgeSets.left   = find(abs(x - 0) < tol);
    edgeSets.right  = find(abs(x - A) < tol);
    edgeSets.bottom = find(abs(y + B) < tol);
    edgeSets.top    = find(abs(y - B) < tol);

    edgeSets.corners = struct();
    edgeSets.corners.left_bottom  = nearest_node(p, [0, -B]);
    edgeSets.corners.right_bottom = nearest_node(p, [A, -B]);
    edgeSets.corners.left_top     = nearest_node(p, [0,  B]);
    edgeSets.corners.right_top    = nearest_node(p, [A,  B]);

    edgeSets.holes = cell(size(holeLoops));
    for k = 1:numel(holeLoops)
        edgeSets.holes{k} = nodes_near_polygon_edges(p, holeLoops{k}, 20*tol);
    end

    if numel(edgeSets.holes) == 1
        edgeSets.hole = edgeSets.holes{1};
    else
        edgeSets.hole = [];
    end
end


function ids = identify_pencil_ids_from_temp_mesh_local(mdl, msh, D, varargin)
% Recover TipVertexID and the two sharp-pencil EdgeIDs from a temporary mesh.
% This is a fallback when pure geometry-edge evaluation is unreliable.

    ip = inputParser;
    addParameter(ip, 'Tol', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(ip, 'Verbose', false, @(x)islogical(x) && isscalar(x));
    parse(ip, varargin{:});

    tolIn   = ip.Results.Tol;
    verbose = ip.Results.Verbose;

    must_field_local(D, 'channelGeom');
    must_field_local(D.channelGeom, 'append');

    Gapp = D.channelGeom.append;
    must_field_local(Gapp, 'xtip');
    must_field_local(Gapp, 'face_upper');
    must_field_local(Gapp, 'face_lower');

    xtip = Gapp.xtip(:).';
    segU = Gapp.face_upper;
    segL = Gapp.face_lower;

    mouthUpper = farther_endpoint_local(segU, xtip);
    mouthLower = farther_endpoint_local(segL, xtip);

    p = msh.Nodes.';
    geom = mdl.Geometry;

    V = get_geom_vertices_2d_local(geom);
    if isempty(V)
        error('mesh_hole_pencil_domain:NoVerticesInFallback', ...
            'Could not access geometry vertices from mdl.Geometry.');
    end

    [v_tip, d_tip] = nearest_vertex_to_point_local(V, xtip);
    [v_up,  d_up ] = nearest_vertex_to_point_local(V, mouthUpper);
    [v_lo,  d_lo ] = nearest_vertex_to_point_local(V, mouthLower);

    if isempty(tolIn)
        scale = max([1, norm(xtip), norm(mouthUpper-xtip), norm(mouthLower-xtip)]);
        tol = 1e-8 * scale;
    else
        tol = tolIn;
    end

    n_tip = unique(findNodes(msh, 'region', 'Vertex', v_tip));
    if isempty(n_tip)
        n_tip = nearest_node(p, xtip);
    else
        n_tip = n_tip(1);
    end

    nEdges = get_num_edges_local(geom);

    candE   = [];
    candFar = [];

    for eid = 1:nEdges
        nEdge = unique(findNodes(msh, 'region', 'Edge', eid));
        if isempty(nEdge)
            continue;
        end

        isIncident = any(nEdge == n_tip) || min(vecnorm(p(nEdge,:) - xtip, 2, 2)) <= 50*tol;
        if ~isIncident
            continue;
        end

        [~, imax] = max(vecnorm(p(nEdge,:) - xtip, 2, 2));
        Pfar = p(nEdge(imax), :);

        candE(end+1,1)   = eid;   %#ok<AGROW>
        candFar(end+1,:) = Pfar;  %#ok<AGROW>
    end

    if numel(candE) < 2
        error('mesh_hole_pencil_domain:TooFewCandidateEdges', ...
            'Could not find two tip-incident candidate edges in the temporary mesh.');
    end

    costU = vecnorm(candFar - mouthUpper, 2, 2);
    costL = vecnorm(candFar - mouthLower, 2, 2);

    bestCost = inf;
    bestPair = [NaN NaN];

    for i = 1:numel(candE)
        for j = 1:numel(candE)
            if i == j
                continue;
            end
            cij = costU(i) + costL(j);
            if cij < bestCost
                bestCost = cij;
                bestPair = [i j];
            end
        end
    end

    e_upper = candE(bestPair(1));
    e_lower = candE(bestPair(2));

    if verbose
        fprintf('identify_pencil_ids_from_temp_mesh: v_tip = %d, dist = %.3e\n', v_tip, d_tip);
        fprintf('identify_pencil_ids_from_temp_mesh: v_up  = %d, dist = %.3e\n', v_up,  d_up);
        fprintf('identify_pencil_ids_from_temp_mesh: v_lo  = %d, dist = %.3e\n', v_lo,  d_lo);
        fprintf('identify_pencil_ids_from_temp_mesh: e_upper = %d\n', e_upper);
        fprintf('identify_pencil_ids_from_temp_mesh: e_lower = %d\n', e_lower);
    end

    ids = struct();
    ids.v_tip  = v_tip;
    ids.v_up   = v_up;
    ids.v_lo   = v_lo;
    ids.e_upper = e_upper;
    ids.e_lower = e_lower;
    ids.e_tip   = [e_upper e_lower];
    ids.xtip    = xtip;
    ids.mouthUpper = mouthUpper;
    ids.mouthLower = mouthLower;
    ids.vertex_dist_tip = d_tip;
    ids.vertex_dist_up  = d_up;
    ids.vertex_dist_lo  = d_lo;
end


function V = get_geom_vertices_2d_local(geom)
    V = [];
    if isprop(geom, 'Vertices')
        VV = geom.Vertices;
        if isnumeric(VV)
            if size(VV,1) == 2
                V = VV.';
            elseif size(VV,2) == 2
                V = VV;
            end
        end
    end
end


function nEdges = get_num_edges_local(geom)
    if isprop(geom, 'NumEdges')
        nEdges = geom.NumEdges;
    else
        error('mesh_hole_pencil_domain:NoNumEdges', ...
            'Could not access geom.NumEdges.');
    end
end


function [ivid, dmin] = nearest_vertex_to_point_local(V, x)
    d = vecnorm(V - x, 2, 2);
    [dmin, ivid] = min(d);
end


function Pfar = farther_endpoint_local(seg, xref)
    d1 = norm(seg(1,:) - xref);
    d2 = norm(seg(2,:) - xref);
    if d1 >= d2
        Pfar = seg(1,:);
    else
        Pfar = seg(2,:);
    end
end


function idx = nearest_node(p, xq)
    d2 = sum((p - xq).^2, 2);
    [~, idx] = min(d2);
end


function ids = nodes_near_polygon_edges(p, V, tol)
    if isempty(V)
        ids = [];
        return;
    end

    if norm(V(end,:) - V(1,:), inf) > 0
        Vc = [V; V(1,:)];
    else
        Vc = V;
    end

    keep = false(size(p,1),1);

    for k = 1:size(Vc,1)-1
        A = Vc(k,:);
        B = Vc(k+1,:);
        keep = keep | point_near_segment(p, A, B, tol);
    end

    ids = find(keep);
end


function mask = point_near_segment(P, A, B, tol)
    AB = B - A;
    L2 = max(dot(AB,AB), 1e-30);

    t = ((P - A) * AB.') / L2;
    t = max(0, min(1, t));

    Q = A + t .* AB;
    d2 = sum((P - Q).^2, 2);

    mask = d2 <= tol^2;
end


function validate_polygon(P, name)
    if ~isnumeric(P) || size(P,2) ~= 2 || size(P,1) < 3
        error('mesh_hole_pencil_domain:BadPolygon', ...
            '%s must be an [N x 2] numeric polygon with N >= 3.', name);
    end
end


function A = signed_polygon_area(P)
    x = P(:,1);
    y = P(:,2);
    x2 = [x(2:end); x(1)];
    y2 = [y(2:end); y(1)];
    A = 0.5 * sum(x .* y2 - x2 .* y);
end


function plot_closed_polygon(P, style, lw)
    plot([P(:,1); P(1,1)], [P(:,2); P(1,2)], style, 'LineWidth', lw);
end


function must(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('mesh_hole_pencil_domain:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end


function must_field_local(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('mesh_hole_pencil_domain:MissingNestedField', ...
            'Required field "%s" is missing or empty.', field);
    end
end