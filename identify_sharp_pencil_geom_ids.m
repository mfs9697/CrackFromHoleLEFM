function ids = identify_sharp_pencil_geom_ids(mdl, D, varargin)
%IDENTIFY_SHARP_PENCIL_GEOM_IDS
% Identify the PDE geometry vertex and the two PDE geometry edges
% corresponding to the sharp appended pencil, BEFORE generateMesh.
%
% Usage:
%   ids = identify_sharp_pencil_geom_ids(mdl, D)
%   ids = identify_sharp_pencil_geom_ids(mdl, D, 'Tol', 1e-8, 'Verbose', true)
%
% Required inputs:
%   mdl : PDE model after geometryFromEdges(mdl, dl)
%   D   : geometry-description struct from build_domain_hole_pencil_polyline
%         Must contain:
%           D.channelGeom.append.xtip
%           D.channelGeom.append.face_upper   = [Mup; xtip]
%           D.channelGeom.append.face_lower   = [xtip; Mlo]
%
% Output:
%   ids struct with fields:
%       .v_tip
%       .v_up
%       .v_lo
%       .e_upper
%       .e_lower
%       .e_tip
%       .xtip
%       .mouthUpper
%       .mouthLower
%       .vertex_dist_tip
%       .vertex_dist_up
%       .vertex_dist_lo
%
% Notes:
%   - This helper uses GEOMETRY labels, not mesh node IDs.
%   - It matches the two sharp pencil edges by endpoint-vertex pairs:
%         upper face <-> {v_up, v_tip}
%         lower face <-> {v_tip, v_lo}
%   - This is intended for the sharp appended-hole geometry.

ip = inputParser;
addParameter(ip, 'Tol', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(ip, 'Verbose', false, @(x)islogical(x) && isscalar(x));
parse(ip, varargin{:});

tolIn   = ip.Results.Tol;
verbose = ip.Results.Verbose;

% ------------------------------------------------------------
% checks
% ------------------------------------------------------------
if ~isfield(D, 'channelGeom') || ~isfield(D.channelGeom, 'append')
    error('identify_sharp_pencil_geom_ids:MissingAppendData', ...
        'D.channelGeom.append is required.');
end

Gapp = D.channelGeom.append;

must_field(Gapp, 'xtip');
must_field(Gapp, 'face_upper');
must_field(Gapp, 'face_lower');

xtip = Gapp.xtip(:).';
segU = Gapp.face_upper;
segL = Gapp.face_lower;

if ~isnumeric(segU) || size(segU,1) ~= 2 || size(segU,2) ~= 2
    error('identify_sharp_pencil_geom_ids:BadUpperFace', ...
        'face_upper must be [2 x 2].');
end
if ~isnumeric(segL) || size(segL,1) ~= 2 || size(segL,2) ~= 2
    error('identify_sharp_pencil_geom_ids:BadLowerFace', ...
        'face_lower must be [2 x 2].');
end

geom = mdl.Geometry;

if isempty(tolIn)
    scale = max([1, norm(xtip), segment_length(segU), segment_length(segL)]);
    tol = 1e-8 * scale;
else
    tol = tolIn;
end

% ------------------------------------------------------------
% target geometry points
% ------------------------------------------------------------
mouthUpper = farther_endpoint(segU, xtip);
mouthLower = farther_endpoint(segL, xtip);

% ------------------------------------------------------------
% geometry vertices
% ------------------------------------------------------------
V = get_geom_vertices_2d(geom);   % [Nv x 2]
if isempty(V)
    error('identify_sharp_pencil_geom_ids:NoVertices', ...
        'Could not access geometry vertices from mdl.Geometry.');
end

[v_tip, d_tip] = nearest_vertex_to_point(V, xtip);
[v_up,  d_up ] = nearest_vertex_to_point(V, mouthUpper);
[v_lo,  d_lo ] = nearest_vertex_to_point(V, mouthLower);

if verbose
    fprintf('identify_sharp_pencil_geom_ids: v_tip = %d, dist = %.3e\n', v_tip, d_tip);
    fprintf('identify_sharp_pencil_geom_ids: v_up  = %d, dist = %.3e\n', v_up,  d_up);
    fprintf('identify_sharp_pencil_geom_ids: v_lo  = %d, dist = %.3e\n', v_lo,  d_lo);
end

if v_tip == v_up || v_tip == v_lo || v_up == v_lo
    error('identify_sharp_pencil_geom_ids:VertexCollision', ...
        ['The identified geometry vertices are not distinct. ', ...
        'Check Tol or the appended-hole geometry.']);
end

% ------------------------------------------------------------
% geometry edges: identify by endpoint vertex pairs
% ------------------------------------------------------------
nEdges = get_num_edges(geom);

e_upper = NaN;
e_lower = NaN;

for eid = 1:nEdges
    [Pend, ok] = sample_geom_edge_endpoints(geom, eid);
    if ~ok || size(Pend,1) ~= 2
        continue;
    end

    [va, da] = nearest_vertex_to_point(V, Pend(1,:));
    [vb, db] = nearest_vertex_to_point(V, Pend(2,:));

    % Skip if one of the sampled endpoints is not close to any geometry vertex
    if da > 100*tol || db > 100*tol
        continue;
    end

    if same_unordered_pair([va vb], [v_up v_tip])
        e_upper = eid;
    elseif same_unordered_pair([va vb], [v_tip v_lo])
        e_lower = eid;
    end
end


% ------------------------------------------------------------
% fallback: identify missing edge(s) from tip-incident candidate edges
% ------------------------------------------------------------
if ~isfinite(e_upper) || ~isfinite(e_lower) || e_upper == e_lower
    candE   = [];
    candFar = [];

    for eid = 1:nEdges
        [Pend, ok] = sample_geom_edge_endpoints(geom, eid);
        if ~ok || size(Pend,1) ~= 2
            continue;
        end

        [va, da] = nearest_vertex_to_point(V, Pend(1,:));
        [vb, db] = nearest_vertex_to_point(V, Pend(2,:));

        if da > 100*tol || db > 100*tol
            continue;
        end

        % Keep only edges incident to the tip vertex (or extremely close to it)
        isIncident = (va == v_tip) || (vb == v_tip) || ...
                     (min(vecnorm(Pend - xtip, 2, 2)) <= 100*tol);

        if ~isIncident
            continue;
        end

        % Characteristic far endpoint of this edge away from the tip
        if va == v_tip && vb ~= v_tip
            Pfar = V(vb,:);
        elseif vb == v_tip && va ~= v_tip
            Pfar = V(va,:);
        else
            [~, imax] = max(vecnorm(Pend - xtip, 2, 2));
            Pfar = Pend(imax,:);
        end

        candE(end+1,1)   = eid;   %#ok<AGROW>
        candFar(end+1,:) = Pfar;  %#ok<AGROW>
    end

    if numel(candE) < 2
        error('identify_sharp_pencil_geom_ids:TooFewCandidateEdges', ...
            'Could not find two tip-incident candidate edges.');
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
end

if ~isfinite(e_upper)
    error('identify_sharp_pencil_geom_ids:UpperEdgeNotFound', ...
        'Could not identify the upper sharp-pencil edge.');
end
if ~isfinite(e_lower)
    error('identify_sharp_pencil_geom_ids:LowerEdgeNotFound', ...
        'Could not identify the lower sharp-pencil edge.');
end
if e_upper == e_lower
    error('identify_sharp_pencil_geom_ids:EdgeCollision', ...
        'Upper and lower sharp-pencil edges were identified as the same edge.');
end



if verbose
    fprintf('identify_sharp_pencil_geom_ids: e_upper = %d\n', e_upper);
    fprintf('identify_sharp_pencil_geom_ids: e_lower = %d\n', e_lower);
end

ids = struct();
ids.v_tip            = v_tip;
ids.v_up             = v_up;
ids.v_lo             = v_lo;
ids.e_upper          = e_upper;
ids.e_lower          = e_lower;
ids.e_tip            = [e_upper e_lower];
ids.xtip             = xtip;
ids.mouthUpper       = mouthUpper;
ids.mouthLower       = mouthLower;
ids.vertex_dist_tip  = d_tip;
ids.vertex_dist_up   = d_up;
ids.vertex_dist_lo   = d_lo;
end


% =========================================================================
% helpers
% =========================================================================

function V = get_geom_vertices_2d(geom)
% Extract geometry vertices as [Nv x 2].

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

if isempty(V)
    try
        VV = geom.Vertices;
        if size(VV,1) == 2
            V = VV.';
        else
            V = VV;
        end
    catch
    end
end
end


function nEdges = get_num_edges(geom)
% Number of geometry edges.

if isprop(geom, 'NumEdges')
    nEdges = geom.NumEdges;
    return;
end

try
    nEdges = geom.NumEdges;
catch
    error('identify_sharp_pencil_geom_ids:NoNumEdges', ...
        'Could not access geom.NumEdges.');
end
end


function [Pend, ok] = sample_geom_edge_endpoints(geom, edgeID)
% Evaluate only the two endpoints of a geometry edge.
% Returns Pend as [2 x 2].

s = [0 1];
Pend = [];
ok = false;

% Try evaluate(geom, edgeID, s)
try
    XY = evaluate(geom, edgeID, s);
    Pend = normalize_xy_output(XY);
    if size(Pend,1) == 2
        ok = true;
        return;
    end
catch
end

% Try evaluate(geom, s, edgeID)
try
    XY = evaluate(geom, s, edgeID);
    Pend = normalize_xy_output(XY);
    if size(Pend,1) == 2
        ok = true;
        return;
    end
catch
end

% Try evaluateEdge(geom, edgeID, s)
try
    XY = evaluateEdge(geom, edgeID, s); %#ok<NASGU>
    Pend = normalize_xy_output(XY);
    if size(Pend,1) == 2
        ok = true;
        return;
    end
catch
end
end


function P = normalize_xy_output(XY)
% Normalize edge-evaluation output to [N x 2].

P = [];

if ~isnumeric(XY) || isempty(XY)
    return;
end

if size(XY,1) == 2
    P = XY.';
elseif size(XY,2) == 2
    P = XY;
end
end


function [ivid, dmin] = nearest_vertex_to_point(V, x)
% Find the nearest geometry vertex to point x.

d = vecnorm(V - x, 2, 2);
[dmin, ivid] = min(d);
end


function tf = same_unordered_pair(a, b)
% True if 1x2 integer arrays a and b match as unordered pairs.

tf = isequal(sort(a(:).'), sort(b(:).'));
end


function Pfar = farther_endpoint(seg, xref)
% Return the segment endpoint farther from xref.

d1 = norm(seg(1,:) - xref);
d2 = norm(seg(2,:) - xref);

if d1 >= d2
    Pfar = seg(1,:);
else
    Pfar = seg(2,:);
end
end


function L = segment_length(seg)
L = norm(seg(2,:) - seg(1,:));
end

function [P, ok] = sample_geom_edge_points(geom, edgeID, nSample)
% Sample points on a geometry edge. Returns P as [N x 2].

    s = linspace(0,1,nSample);
    P = [];
    ok = false;

    try
        XY = evaluate(geom, edgeID, s);
        P = normalize_xy_output(XY);
        if size(P,1) >= 2
            ok = true;
            return;
        end
    catch
    end

    try
        XY = evaluate(geom, s, edgeID);
        P = normalize_xy_output(XY);
        if size(P,1) >= 2
            ok = true;
            return;
        end
    catch
    end

    try
        XY = evaluateEdge(geom, edgeID, s); %#ok<NASGU>
        P = normalize_xy_output(XY);
        if size(P,1) >= 2
            ok = true;
            return;
        end
    catch
    end
end


function d = point_to_segment_distance(P, A, B)
% Distance from each point row in P to the segment AB.

    AB = B - A;
    L2 = max(dot(AB, AB), 1e-30);

    t = ((P - A) * AB.') / L2;
    t = max(0, min(1, t));

    Q = A + t .* AB;
    d = vecnorm(P - Q, 2, 2);
end
