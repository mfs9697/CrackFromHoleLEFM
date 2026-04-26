function D = build_domain_hole_pencil_polyline(Pmid, A, B, holes, w, varargin)
%BUILD_DOMAIN_HOLE_PENCIL_POLYLINE
% Build geometric description for a rectangular plate with one or more holes,
% where the touched circular hole is replaced by a merged sharp appended-hole loop.
%
%   D = build_domain_hole_pencil_polyline(Pmid, A, B, holes, w)
%
% Inputs:
%   Pmid   : (N x 2) crack midline vertices, ordered mouth -> ... -> tip
%   A, B   : plate dimensions, domain x in [0,A], y in [-B,B]
%   holes  : cell array of hole structs, e.g.
%              { struct('type','circle','center',[xc yc],'r',R,'npoly',N), ... }
%            or empty {}
%   w      : appended-hole mouth half-shift (formerly legacy channel half-width)
%
% Name-value options:
%   'corner_tol'  : positive scalar (default 1e-10)
%   'mouth_eps'   : mouth half-shift used for appended-hole construction
%                   default = w
%   'epsMode'     : 'arclength' (default) or 'angle'
%   'nArc'        : retained-hole arc resolution (default 160)
%   'orientation' : 'cw' (default) or 'ccw' for appended inner loop
%
% Compatibility-only options (accepted but ignored):
%   'join'
%   'miter_limit'
%   'tip'
%   'mode'
%   'nFace'
%   'nTip'
%   'tipRadius'
%
% Output:
%   D      : geometry-description struct with fields
%       .outerPoly       rectangle polygon [4 x 2], CCW
%       .holeLoops       cell array of inner polygons
%       .channelPoly     []   (kept only for downstream compatibility)
%       .Pmid            crack midline
%       .A, .B, .w       copied inputs
%       .holes           normalized hole list
%       .channelGeom     debug struct
%       .topology        auxiliary geometric/topological info
%
% Notes:
%   - This version is merged-only: there is no separate channel geometry.
%   - The crack must start on a supported circular hole.

    % ------------------------------------------------------------
    % checks
    % ------------------------------------------------------------
    validateattributes(Pmid, {'numeric'}, {'2d','ncols',2,'nonempty','finite'}, ...
        mfilename, 'Pmid', 1);
    validateattributes(A, {'numeric'}, {'scalar','real','positive','finite'}, ...
        mfilename, 'A', 2);
    validateattributes(B, {'numeric'}, {'scalar','real','positive','finite'}, ...
        mfilename, 'B', 3);
    validateattributes(w, {'numeric'}, {'scalar','real','positive','finite'}, ...
        mfilename, 'w', 5);

    if size(Pmid,1) < 2
        error('build_domain_hole_pencil_polyline:BadPmid', ...
            'Pmid must contain at least two points.');
    end

    holes = normalize_holes(holes);

    % ------------------------------------------------------------
    % parse options
    % ------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'corner_tol',  1e-10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'mouth_eps',   [],     @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'epsMode',     'arclength', @(s) ischar(s) || isstring(s));
    addParameter(p, 'nArc',        160,    @(x) isnumeric(x) && isscalar(x) && x >= 8);
    addParameter(p, 'orientation', 'cw',   @(s) ischar(s) || isstring(s));

    % compatibility-only options: accepted but ignored
    addParameter(p, 'join',        [], @(x) true);
    addParameter(p, 'miter_limit', [], @(x) true);
    addParameter(p, 'tip',         [], @(x) true);
    addParameter(p, 'mode',        [], @(x) true);
    addParameter(p, 'nFace',       [], @(x) true);
    addParameter(p, 'nTip',        [], @(x) true);
    addParameter(p, 'tipRadius',   [], @(x) true);

    parse(p, varargin{:});

    corner_tol = p.Results.corner_tol;
    epsMode    = tolower_safe(p.Results.epsMode);
    orient     = tolower_safe(p.Results.orientation);
    nArc       = p.Results.nArc;

    mouth_eps = p.Results.mouth_eps;
    if isempty(mouth_eps)
        mouth_eps = w;
    end

    % ------------------------------------------------------------
    % outer rectangle
    % ------------------------------------------------------------
    outerPoly = [ ...
        0, -B;
        A, -B;
        A,  B;
        0,  B ];

    % ------------------------------------------------------------
    % initial hole loops
    % ------------------------------------------------------------
    holeLoops = holes_to_loops(holes);

    % ------------------------------------------------------------
    % topology detection
    % ------------------------------------------------------------
    topo = struct();

    topo.Pstart = Pmid(1,:);
    topo.Ptip   = Pmid(end,:);

    topo.channelTouchesHole = false;
    topo.touchingHoleIndex  = [];
    topo.touchDistance      = NaN;

    topo.channelArea        = NaN;   % compatibility only
    topo.appendedHoleArea   = NaN;

    if isempty(holes)
        error('build_domain_hole_pencil_polyline:NoHoles', ...
            'At least one circular hole is required for the appended-hole route.');
    end

    [itouch, dtouch] = detect_hole_touch(Pmid(1,:), holes);
    topo.channelTouchesHole = ~isempty(itouch);
    topo.touchingHoleIndex  = itouch;
    topo.touchDistance      = dtouch;

    if isempty(itouch)
        error('build_domain_hole_pencil_polyline:NoTouchedHole', ...
            'Pmid(1,:) was not detected on a supported hole.');
    end

    hk = holes{itouch};
    if ~strcmpi(strtrim(hk.type), 'circle')
        error('build_domain_hole_pencil_polyline:UnsupportedMergedHole', ...
            'Only circular holes are currently supported.');
    end

    if exist('build_appended_hole_loop', 'file') ~= 2
        error('build_domain_hole_pencil_polyline:MissingHelper', ...
            'build_appended_hole_loop.m is required on the MATLAB path.');
    end

    % ------------------------------------------------------------
    % infer local frame from the crack polyline
    % ------------------------------------------------------------
    [A0, n_in, t_hat, theta, a0, Gframe] = ...
        infer_local_frame_from_pmid(hk, Pmid, corner_tol);

    % ------------------------------------------------------------
    % build merged appended-hole loop
    % ------------------------------------------------------------
    [Vapp, Gapp] = build_appended_hole_loop(hk, A0, n_in, t_hat, theta, a0, mouth_eps, ...
        'epsMode', epsMode, ...
        'nArc', nArc, ...
        'orientation', orient);

    holeLoops{itouch} = Vapp;

    channelGeom = struct();
    channelGeom.mode   = 'merged_appended_hole';
    channelGeom.frame  = Gframe;
    channelGeom.append = Gapp;

    topo.mode               = 'merged_appended_hole';
    topo.appendedHoleIndex  = itouch;
    topo.appendedHoleArea   = abs(signed_polygon_area(Vapp));

    % ------------------------------------------------------------
    % diagnostics / flags
    % ------------------------------------------------------------
    topo.bbox = struct();
    topo.bbox.outer = polygon_bbox(outerPoly);
    topo.bbox.holes = cell(size(holeLoops));
    for k = 1:numel(holeLoops)
        topo.bbox.holes{k} = polygon_bbox(holeLoops{k});
    end
    topo.bbox.channel = [];   % compatibility only

    topo.flags = struct();
    topo.flags.tipInsidePlate = point_in_box(Pmid(end,:), [0, A, -B, B], 1e-12);
    topo.flags.midlineInsidePlate = all(Pmid(:,1) >= -1e-12 & Pmid(:,1) <= A + 1e-12 & ...
                                        Pmid(:,2) >= -B - 1e-12 & Pmid(:,2) <= B + 1e-12);
    topo.flags.startOnHole = topo.channelTouchesHole;
    topo.flags.channelInsidePlateBox = [];   % compatibility only

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    D = struct();

    D.outerPoly   = outerPoly;
    D.holeLoops   = holeLoops;
    D.channelPoly = [];       % compatibility only

    D.Pmid        = Pmid;
    D.A           = A;
    D.B           = B;
    D.w           = w;
    D.holes       = holes;

    D.channelGeom = channelGeom;
    D.topology    = topo;

    D.options = struct();
    D.options.corner_tol  = corner_tol;
    D.options.mouth_eps   = mouth_eps;
    D.options.epsMode     = epsMode;
    D.options.nArc        = nArc;
    D.options.orientation = orient;
end


% =========================================================================
% merged-hole frame inference
% =========================================================================

function [A0, n_in, t_hat, theta, a0, G] = infer_local_frame_from_pmid(hole, Pmid, tol)
% Infer the local hole frame and relative crack angle from the first crack segment.

    c = hole.center(:).';
    r = hole.r;

    A = Pmid(1,:);
    rc = A - c;
    nr = norm(rc);
    if nr <= tol
        error('build_domain_hole_pencil_polyline:BadStartPoint', ...
            'Pmid(1,:) is too close to the hole center.');
    end

    % Project mouth point onto the circle for robustness
    A0 = c + r * rc / nr;

    % Inward normal into the material
    n_in = normalize_row(A0 - c);

    % Circle tangent, sign chosen to align with the first crack segment
    t_circ = [-n_in(2), n_in(1)];

    d0 = Pmid(2,:) - Pmid(1,:);
    if norm(d0) <= tol
        error('build_domain_hole_pencil_polyline:DegenerateMouthSegment', ...
            'The first segment of Pmid is degenerate.');
    end
    e_dir = normalize_row(d0);

    if dot(t_circ, e_dir) >= 0
        t_hat = t_circ;
    else
        t_hat = -t_circ;
    end

    % Relative angle in the local frame
    theta = atan2(dot(e_dir, t_hat), dot(e_dir, n_in));

    % Straight effective appendix length
    a0 = norm(Pmid(end,:) - Pmid(1,:));

    G = struct();
    G.A_input = A;
    G.A = A0;
    G.center = c;
    G.r = r;
    G.e_dir = e_dir;
    G.n_in = n_in;
    G.t_hat = t_hat;
    G.theta = theta;
    G.a0 = a0;
end


% =========================================================================
% helpers
% =========================================================================

function holes = normalize_holes(holes)
    if isempty(holes)
        holes = {};
        return;
    end

    if isstruct(holes)
        holes = num2cell(holes);
    end

    if ~iscell(holes)
        error('build_domain_hole_pencil_polyline:BadHoleInput', ...
            'holes must be a cell array, struct array, or empty.');
    end

    for k = 1:numel(holes)
        hk = holes{k};
        if ~isstruct(hk) || ~isfield(hk, 'type')
            error('build_domain_hole_pencil_polyline:BadHoleSpec', ...
                'Each hole must be a struct with field "type".');
        end
    end
end


function [itouch, dmin] = detect_hole_touch(x0, holes)
    itouch = [];
    dmin = NaN;

    best = inf;

    for k = 1:numel(holes)
        hk = holes{k};

        switch lower(strtrim(hk.type))
            case 'circle'
                must_field(hk, 'center');
                must_field(hk, 'r');

                c = hk.center(:).';
                r = hk.r;

                d = abs(norm(x0 - c) - r);

            case 'polygon'
                must_field(hk, 'vertices');
                V = hk.vertices;
                d = point_polygon_distance(x0, V);

            otherwise
                error('build_domain_hole_pencil_polyline:UnsupportedHoleType', ...
                    'Unsupported hole type "%s".', hk.type);
        end

        if d < best
            best = d;
            itouch = k;
        end
    end

    tol = 1e-8 * max(1, norm(x0));
    if best <= tol
        dmin = best;
    else
        itouch = [];
        dmin = best;
    end
end


function bbox = polygon_bbox(P)
    bbox = [min(P(:,1)), max(P(:,1)), min(P(:,2)), max(P(:,2))];
end


function tf = point_in_box(x, box, tol)
    tf = (x(1) >= box(1) - tol) && (x(1) <= box(2) + tol) && ...
         (x(2) >= box(3) - tol) && (x(2) <= box(4) + tol);
end


function A = signed_polygon_area(P)
    x = P(:,1);
    y = P(:,2);
    x2 = [x(2:end); x(1)];
    y2 = [y(2:end); y(1)];
    A = 0.5 * sum(x .* y2 - x2 .* y);
end


function d = point_polygon_distance(x0, V)
    if isempty(V)
        d = inf;
        return;
    end

    if norm(V(end,:) - V(1,:), inf) > 0
        V = [V; V(1,:)];
    end

    d = inf;
    for k = 1:size(V,1)-1
        A = V(k,:);
        B = V(k+1,:);
        d = min(d, point_segment_distance(x0, A, B));
    end
end


function d = point_segment_distance(P, A, B)
    AB = B - A;
    L2 = dot(AB, AB);

    if L2 <= 0
        d = norm(P - A);
        return;
    end

    t = dot(P - A, AB) / L2;
    t = max(0, min(1, t));

    Q = A + t * AB;
    d = norm(P - Q);
end


function v = normalize_row(v)
    v = v(:).';
    nv = norm(v);
    if nv <= 0
        error('build_domain_hole_pencil_polyline:ZeroVector', ...
            'Cannot normalize a zero vector.');
    end
    v = v / nv;
end


function must_field(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('build_domain_hole_pencil_polyline:MissingHoleField', ...
            'Required field "%s" is missing or empty.', field);
    end
end


function s = tolower_safe(x)
    if isstring(x)
        x = char(x);
    end
    if ~ischar(x)
        error('build_domain_hole_pencil_polyline:BadStringInput', ...
            'Expected char or string input.');
    end

    s = x;
    mask = (s >= 'A') & (s <= 'Z');
    s(mask) = char(double(s(mask)) + ('a' - 'A'));
end