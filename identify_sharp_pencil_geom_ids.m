function ids = identify_sharp_pencil_geom_ids(mdl, D, varargin)
%IDENTIFY_SHARP_PENCIL_GEOM_IDS
% Identify the PDE geometry vertex and the two PDE geometry edges
% corresponding to the sharp appended pencil, BEFORE generateMesh.
%
% Usage:
%   ids = identify_sharp_pencil_geom_ids(mdl, D)
%   ids = identify_sharp_pencil_geom_ids(mdl, D, 'NSample', 25, 'Tol', 1e-8)
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
%       .e_upper
%       .e_lower
%       .xtip
%       .face_upper
%       .face_lower
%       .vertex_dist
%       .edge_score_upper
%       .edge_score_lower
%
% Notes:
%   - This helper uses geometry labels, not mesh node IDs.
%   - It tries to sample PDE edges through the geometry object.
%   - If your MATLAB release does not support the required geometry-edge
%     evaluation calls, use PlotGeom with VertexLabels/EdgeLabels once
%     and compare manually.

    ip = inputParser;
    addParameter(ip, 'NSample', 25, @(x)isnumeric(x)&&isscalar(x)&&x>=3);
    addParameter(ip, 'Tol', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(ip, 'Verbose', false, @(x)islogical(x)&&isscalar(x));
    parse(ip, varargin{:});

    nSample = ip.Results.NSample;
    tolIn   = ip.Results.Tol;
    verbose = ip.Results.Verbose;

    % -------------------- checks --------------------
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

    % -------------------- tolerance --------------------
    if isempty(tolIn)
        scale = max([1, norm(xtip), segment_length(segU), segment_length(segL)]);
        tol = 1e-6 * scale;
    else
        tol = tolIn;
    end

    % -------------------- vertex identification --------------------
    V = get_geom_vertices_2d(geom);  % [Nv x 2]
    if isempty(V)
        error('identify_sharp_pencil_geom_ids:NoVertices', ...
            'Could not access geometry vertices from mdl.Geometry.');
    end

    dV = vecnorm(V - xtip, 2, 2);
    [dminV, v_tip] = min(dV);

    if verbose
        fprintf('identify_sharp_pencil_geom_ids: nearest vertex to tip = %d, dist = %.3e\n', ...
            v_tip, dminV);
    end

    % -------------------- edge identification --------------------
    nEdges = get_num_edges(geom);

    scoreU = inf(nEdges,1);
    scoreL = inf(nEdges,1);

    for eid = 1:nEdges
        [Pe, ok] = sample_geom_edge_points(geom, eid, nSample);
        if ~ok || size(Pe,1) < 2
            continue;
        end

        scoreU(eid) = edge_to_segment_score(Pe, segU);
        scoreL(eid) = edge_to_segment_score(Pe, segL);
    end

    if ~any(isfinite(scoreU)) || ~any(isfinite(scoreL))
        error('identify_sharp_pencil_geom_ids:NoFiniteEdgeMatches', ...
        'Could not identify appended-pencil edges from the PDE geometry.');
    end

    [bestU, e_upper] = min(scoreU);
    [bestL, e_lower] = min(scoreL);

    if ~isfinite(bestU) || ~isfinite(bestL)
        error('identify_sharp_pencil_geom_ids:BadEdgeMatches', ...
            'Edge identification failed: non-finite score.');
    end

    % If both segments accidentally chose the same edge, resolve by second best
    if e_upper == e_lower
        scoreL2 = scoreL;
        scoreL2(e_upper) = inf;
        [bestL2, e_lower2] = min(scoreL2);

        scoreU2 = scoreU;
        scoreU2(e_lower) = inf;
        [bestU2, e_upper2] = min(scoreU2);

        % pick the better of the two conflict resolutions
        if bestU + bestL2 <= bestL + bestU2
            e_lower = e_lower2;
            bestL   = bestL2;
        else
            e_upper = e_upper2;
            bestU   = bestU2;
        end
    end

    if verbose
        fprintf('identify_sharp_pencil_geom_ids: upper edge = %d, score = %.3e\n', ...
            e_upper, bestU);
        fprintf('identify_sharp_pencil_geom_ids: lower edge = %d, score = %.3e\n', ...
            e_lower, bestL);
    end

    ids = struct();
    ids.v_tip            = v_tip;
    ids.e_upper          = e_upper;
    ids.e_lower          = e_lower;
    ids.xtip             = xtip;
    ids.face_upper       = segU;
    ids.face_lower       = segL;
    ids.vertex_dist      = dminV;
    ids.edge_score_upper = bestU;
    ids.edge_score_lower = bestL;
end


% =========================================================================
% helpers
% =========================================================================

function V = get_geom_vertices_2d(geom)
% Try to extract geometry vertices as [Nv x 2].

    V = [];

    % Common case
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

    % Some releases may store vertices differently
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

    % fallback
    try
        nEdges = geom.NumEdges;
    catch
        error('identify_sharp_pencil_geom_ids:NoNumEdges', ...
            'Could not access geom.NumEdges.');
    end
end


function [P, ok] = sample_geom_edge_points(geom, edgeID, nSample)
% Try several calling conventions to evaluate points on a geometry edge.
% Returns sampled points P as [nSample x 2].

    s = linspace(0,1,nSample);
    P = [];
    ok = false;

    % Try evaluate(geom, edgeID, s)
    try
        XY = evaluate(geom, edgeID, s);
        P = normalize_xy_output(XY);
        if size(P,1) >= 2
            ok = true;
            return;
        end
    catch
    end

    % Try evaluate(geom, s, edgeID)
    try
        XY = evaluate(geom, s, edgeID);
        P = normalize_xy_output(XY);
        if size(P,1) >= 2
            ok = true;
            return;
        end
    catch
    end

    % Try parameterized edge evaluation via geometry methods if available
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


function P = normalize_xy_output(XY)
% Normalize sampled edge output to [N x 2].

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


function score = edge_to_segment_score(Pe, seg)
% Score how well sampled geometry edge Pe matches target segment seg=[A;B].

    A = seg(1,:);
    B = seg(2,:);

    % sample-to-segment distance
    dCurve = point_to_segment_distance(Pe, A, B);
    meanCurve = mean(dCurve);

    % endpoint match score, allowing reversed orientation
    E1 = Pe(1,:);
    E2 = Pe(end,:);

    sdir = norm(E1 - A) + norm(E2 - B);
    srev = norm(E1 - B) + norm(E2 - A);
    endScore = min(sdir, srev);

    % midpoint proximity is a useful stabilizer
    Mseg = 0.5*(A + B);
    Medg = 0.5*(E1 + E2);
    midScore = norm(Medg - Mseg);

    score = endScore + 0.5*midScore + 0.25*meanCurve;
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


function L = segment_length(seg)
    L = norm(seg(2,:) - seg(1,:));
end


function must_field(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('identify_sharp_pencil_geom_ids:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end