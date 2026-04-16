function M = mesh_hole_pencil_domain(D, varargin)
%MESH_HOLE_PENCIL_DOMAIN  Mesh geometry: rectangle - holes - pencil channel.
%
%   M = mesh_hole_pencil_domain(D)
%   M = mesh_hole_pencil_domain(D, 'Hmax', 0.005, 'Hgrad', 1.2)
%
% Input:
%   D   geometry-description struct from build_domain_hole_pencil_polyline
%       Required fields:
%         .outerPoly    [No x 2] outer boundary polygon (CCW)
%         .holeLoops    cell array of hole polygons
%         .channelPoly  [Nc x 2] local channel polygon
%
% Name-value options:
%   'Hmax'      : global max element size (default = [])
%   'Hgrad'     : mesh growth factor      (default = [])
%   'Hmin'      : optional min element size (default = [])
%   'GeometricOrder' : 'linear' (default) or 'quadratic'
%   'PlotGeom'  : logical, plot PDE geometry labels (default false)
%   'PlotMesh'  : logical, plot mesh (default false)
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
%   1) The domain meshed here is
%          outer rectangle \ (holes U channel)
%      i.e. the channel interior is removed.
%   2) This creates a slit-like removed zone, but PDE Toolbox does NOT
%      preserve paired crack-face topology explicitly. So this is suitable
%      as a first meshing stage / geometry-debug stage.
%   3) If later you need duplicated nodes on opposite crack faces for
%      strict LEFM crack-face treatment, a custom topology mesher will still
%      be needed. This function is the correct first bridge.
%
% Geometry conventions:
%   - outerPoly should be CCW
%   - hole loops and channel polygon are treated as interior voids

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
    parse(ip, varargin{:});

    Hmax = ip.Results.Hmax;
    Hgrad = ip.Results.Hgrad;
    Hmin = ip.Results.Hmin;
    geometricOrder = char(ip.Results.GeometricOrder);
    plotGeom = ip.Results.PlotGeom;
    plotMesh = ip.Results.PlotMesh;

    % ------------------------------------------------------------
    % checks
    % ------------------------------------------------------------
    must(D, 'outerPoly');
    must(D, 'holeLoops');
    must(D, 'channelPoly');

    outerPoly  = D.outerPoly;
    holeLoops  = D.holeLoops;
    channelPoly = D.channelPoly;

    validate_polygon(outerPoly,  'outerPoly');
    validate_polygon(channelPoly,'channelPoly');
    for k = 1:numel(holeLoops)
        validate_polygon(holeLoops{k}, sprintf('holeLoops{%d}', k));
    end

    % enforce orientations:
    %   outer boundary -> CCW
    %   inner void loops -> CW preferred for bookkeeping, though decsg uses SF
    if signed_polygon_area(outerPoly) < 0
        outerPoly = flipud(outerPoly);
    end

    for k = 1:numel(holeLoops)
        if signed_polygon_area(holeLoops{k}) > 0
            holeLoops{k} = flipud(holeLoops{k});
        end
    end

    if signed_polygon_area(channelPoly) > 0
        channelPoly = flipud(channelPoly);
    end

    % ------------------------------------------------------------
    % build decomposed geometry:
    %   R1 - H1 - H2 - ... - C1
    % ------------------------------------------------------------
    polys = [{outerPoly}, holeLoops(:).', {channelPoly}];
    names = cell(1, numel(polys));

    names{1} = 'R1';
    for k = 2:numel(polys)-1
        names{k} = sprintf('H%d', k-1);
    end
    names{end} = 'C1';

    [gd, ns] = polygons_to_gd_ns(polys, names);

    sf = 'R1';
    for k = 2:numel(names)
        sf = [sf, '-', names{k}]; %#ok<AGROW>
    end

    [dl, bt] = decsg(gd, sf, ns);

    mdl = createpde();
    geometryFromEdges(mdl, dl);

    if plotGeom
        figure('Name', 'mesh\_hole\_pencil\_domain: geometry', 'Color', 'w'); clf
        pdegplot(mdl, 'EdgeLabels', 'on', 'FaceLabels', 'on');
        axis equal
        title('PDE geometry with edge/face labels')
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

    msh = generateMesh(mdl, gmArgs{:});

    p = msh.Nodes.';      % [np x 2]
    t = msh.Elements.';   % [nt x 3] or [nt x 6]

    % ------------------------------------------------------------
    % basic boundary node sets
    % ------------------------------------------------------------
    edgeSets = build_edge_sets_from_geometry(p, D);

    % ------------------------------------------------------------
    % region metadata
    % ------------------------------------------------------------
    region = struct();
    region.outerPoly   = outerPoly;
    region.holeLoops   = holeLoops;
    region.channelPoly = channelPoly;
    region.names       = names;
    region.nHoles      = numel(holeLoops);

    % ------------------------------------------------------------
    % optional mesh plot
    % ------------------------------------------------------------
    if plotMesh
        figure('Name', 'mesh\_hole\_pencil\_domain: mesh', 'Color', 'w'); clf
        triplot(t(:,1:3), p(:,1), p(:,2), 'Color', [0.75 0.75 0.75]); hold on
        axis equal; box on

        plot_closed_polygon(outerPoly, 'k-', 1.2);
        for k = 1:numel(holeLoops)
            plot_closed_polygon(holeLoops{k}, 'r-', 1.2);
        end
        plot_closed_polygon(channelPoly, 'b-', 1.4);

        title('Mesh of rectangle - holes - channel')
        xlabel('x'); ylabel('y')
    end

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    M = struct();
    M.p        = p;
    M.t        = t;
    M.geom     = mdl;
    M.meshobj  = msh;

    M.dl       = dl;
    M.bt       = bt;
    M.gd       = gd;
    M.ns       = ns;
    M.sf       = sf;

    M.edgeSets = edgeSets;
    M.region   = region;
end


% =========================================================================
% helpers
% =========================================================================

function [gd, ns] = polygons_to_gd_ns(polys, names)
%POLYGONS_TO_GD_NS Build decsg matrices for a list of polygons.
%
% Each polygon is encoded as a decsg polygon object:
%   gd(:,j) = [2; Nv; x1;...;xNv; y1;...;yNv; zeros-padding]
%
% All columns are padded to equal length.

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

        % remove duplicated closing point if present
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


function edgeSets = build_edge_sets_from_geometry(p, D)
%BUILD_EDGE_SETS_FROM_GEOMETRY Build practical node sets from geometry.
%
% Fields:
%   .left, .right, .top, .bottom
%   .holes{k}
%   .channel
%   .corners.left_bottom, right_bottom, left_top, right_top

    edgeSets = struct();

    % outer box
    if isfield(D, 'A') && isfield(D, 'B')
        A = D.A;
        B = D.B;
    else
        bb = [min(D.outerPoly(:,1)), max(D.outerPoly(:,1)), ...
              min(D.outerPoly(:,2)), max(D.outerPoly(:,2))];
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

    % holes
    edgeSets.holes = cell(size(D.holeLoops));
    for k = 1:numel(D.holeLoops)
        edgeSets.holes{k} = nodes_near_polygon_edges(p, D.holeLoops{k}, 20*tol);
    end

    % backward-compatible single-hole field if only one hole exists
    if numel(edgeSets.holes) == 1
        edgeSets.hole = edgeSets.holes{1};
    else
        edgeSets.hole = [];
    end

    % channel boundary nodes
    edgeSets.channel = nodes_near_polygon_edges(p, D.channelPoly, 20*tol);
end


function idx = nearest_node(p, xq)
    d2 = sum((p - xq).^2, 2);
    [~, idx] = min(d2);
end


function ids = nodes_near_polygon_edges(p, V, tol)
%NODES_NEAR_POLYGON_EDGES Return node ids near polygon edges.

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
%SIGNED_POLYGON_AREA Positive for CCW.

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