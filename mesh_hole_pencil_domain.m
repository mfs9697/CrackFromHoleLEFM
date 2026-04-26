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
%   - If manual Hvertex/Hedge is supplied, automatic tip identification is skipped.
%   - If manual Hvertex/Hedge is not supplied, the function attempts to
%     identify the sharp-pencil tip and two appended edges automatically.

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

    % ------------------------------------------------------------
    % optional automatic sharp-tip identification
    % ------------------------------------------------------------
    tipIDs = [];
    useManualRefinement = ~isempty(Hvertex) || ~isempty(Hedge);

    if ~useManualRefinement && isfield(D, 'channelGeom') && isfield(D.channelGeom, 'append')
        tipIDs = identify_sharp_pencil_geom_ids(mdl, D, 'Verbose', true);
    end

    if plotGeom
        figure('Name', 'mesh_hole_pencil_domain: geometry', 'Color', 'w'); clf
        pdegplot(mdl, 'EdgeLabels', 'on', 'FaceLabels', 'on', 'VertexLabels', 'on');
        axis equal
        title('PDE geometry: rectangle - appended-hole loops')
    end

    % ------------------------------------------------------------
    % mesh generation
    % ------------------------------------------------------------
    gmArgs = {'GeometricOrder', geometricOrder};

    if ~isempty(Hmax)
        gmArgs = [gmArgs, {'Hmax', Hmax}];
    end
    if ~isempty(Hgrad)
        gmArgs = [gmArgs, {'Hgrad', Hgrad}];
    end
    if ~isempty(Hmin)
        gmArgs = [gmArgs, {'Hmin', Hmin}];
    end

    if ~isempty(Hedge)
        gmArgs = [gmArgs, {'Hedge', Hedge}];
    end
    if ~isempty(Hvertex)
        gmArgs = [gmArgs, {'Hvertex', Hvertex}];
    end

    if ~isempty(tipIDs)
        [Hface, Htip] = choose_auto_refinement_sizes(Hmin, Hmax);
        gmArgs = [gmArgs, {'Hedge',   {[tipIDs.e_upper tipIDs.e_lower], Hface}}];
        gmArgs = [gmArgs, {'Hvertex', {[tipIDs.v_tip], Htip}}];
    end

    msh = generateMesh(mdl, gmArgs{:});

    p = msh.Nodes.';      % [np x 2]
    t = msh.Elements.';   % [nt x 3] or [nt x 6]

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

function plot_closed_polygon(P, style, lw)
    plot([P(:,1); P(1,1)], [P(:,2); P(1,2)], style, 'LineWidth', lw);
end


function must(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('mesh_hole_pencil_domain:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end