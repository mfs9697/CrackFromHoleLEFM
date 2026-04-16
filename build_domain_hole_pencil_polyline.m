function D = build_domain_hole_pencil_polyline(Pmid, A, B, holes, w, varargin)
%BUILD_DOMAIN_HOLE_PENCIL_POLYLINE
% Build geometric description for a rectangular plate with:
%   - one or more holes,
%   - a local pencil-like channel around a crack polyline Pmid.
%
%   D = build_domain_hole_pencil_polyline(Pmid, A, B, holes, w)
%   D = build_domain_hole_pencil_polyline(..., 'tip', 'point', ...)
%
% Inputs:
%   Pmid   : (N x 2) crack midline vertices, ordered mouth -> ... -> tip
%   A, B   : plate dimensions, domain x in [0,A], y in [-B,B]
%   holes  : cell array of hole structs, e.g.
%              { struct('type','circle','center',[xc yc],'r',R,'npoly',N), ... }
%            or empty {}
%   w      : channel half-width
%
% Name-value options:
%   'join'        : 'miter' (default) or 'bevel'
%   'miter_limit' : positive scalar (default 6)  [currently accepted, lightly used]
%   'corner_tol'  : positive scalar (default 1e-10)
%   'tip'         : 'point' (default) or 'flat'
%
% Output:
%   D      : geometry-description struct with fields
%       .outerPoly       rectangle polygon [4 x 2], CCW
%       .holeLoops       cell array of hole polygons
%       .channelPoly     local channel polygon, CCW
%       .Pmid            crack midline
%       .A, .B, .w       copied inputs
%       .holes           normalized hole list
%       .channelGeom     debug struct for the local channel construction
%       .topology        auxiliary geometric/topological info
%
% Notes:
%   - This function builds geometry only; it does not mesh.
%   - Unlike the old edge-crack builder, this version does NOT connect the
%     channel to the outer rectangle boundary.
%   - The future mesher should interpret the domain as
%
%         plate minus holes minus channel interior
%
%     with the channel mouth connected to the hole boundary.
%
%   - The start point Pmid(1,:) is expected to lie on the hole boundary.

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
    addParameter(p, 'join',        'miter', @(s) ischar(s) || isstring(s));
    addParameter(p, 'miter_limit', 6,       @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'corner_tol',  1e-10,   @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'tip',         'point', @(s) ischar(s) || isstring(s));
    parse(p, varargin{:});

    join       = tolower_safe(p.Results.join);
    tipMode    = tolower_safe(p.Results.tip);

    miter_limit = p.Results.miter_limit; %#ok<NASGU>
    corner_tol  = p.Results.corner_tol;

    % ------------------------------------------------------------
    % outer rectangle
    % ------------------------------------------------------------
    outerPoly = [ ...
        0, -B;
        A, -B;
        A,  B;
        0,  B ];

    % ------------------------------------------------------------
    % hole loops
    % ------------------------------------------------------------
    holeLoops = holes_to_loops(holes);

    % ------------------------------------------------------------
    % local channel polygon
    % ------------------------------------------------------------
    [channelPoly, Gc] = build_local_channel_polygon(Pmid, w, ...
        'join', join, ...
        'corner_tol', corner_tol, ...
        'tip', tipMode);

    % ensure CCW orientation
    if signed_polygon_area(channelPoly) < 0
        channelPoly = flipud(channelPoly);
    end

    % ------------------------------------------------------------
    % topology and diagnostics
    % ------------------------------------------------------------
    topo = struct();

    topo.Pstart = Pmid(1,:);
    topo.Ptip   = Pmid(end,:);

    topo.channelTouchesHole = false;
    topo.touchingHoleIndex  = [];
    topo.touchDistance      = NaN;

    if ~isempty(holes)
        [itouch, dtouch] = detect_hole_touch(Pmid(1,:), holes);
        topo.channelTouchesHole = ~isempty(itouch);
        topo.touchingHoleIndex  = itouch;
        topo.touchDistance      = dtouch;
    end

    topo.channelArea = polyarea(channelPoly(:,1), channelPoly(:,2));

    topo.bbox.outer   = polygon_bbox(outerPoly);
    topo.bbox.channel = polygon_bbox(channelPoly);
    topo.bbox.holes   = cell(size(holeLoops));
    for k = 1:numel(holeLoops)
        topo.bbox.holes{k} = polygon_bbox(holeLoops{k});
    end

    topo.flags = struct();
    topo.flags.tipInsidePlate = point_in_box(Pmid(end,:), [0, A, -B, B], 1e-12);
    topo.flags.midlineInsidePlate = all(Pmid(:,1) >= -1e-12 & Pmid(:,1) <= A + 1e-12 & ...
                                        Pmid(:,2) >= -B - 1e-12 & Pmid(:,2) <= B + 1e-12);
    topo.flags.startOnHole = topo.channelTouchesHole;

    % crude overlap checks
    topo.flags.channelInsidePlateBox = all(channelPoly(:,1) >= -1e-12 & channelPoly(:,1) <= A + 1e-12 & ...
                                           channelPoly(:,2) >= -B - 1e-12 & channelPoly(:,2) <= B + 1e-12);

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    D = struct();

    D.outerPoly   = outerPoly;
    D.holeLoops   = holeLoops;
    D.channelPoly = channelPoly;

    D.Pmid        = Pmid;
    D.A           = A;
    D.B           = B;
    D.w           = w;
    D.holes       = holes;

    D.channelGeom = Gc;
    D.topology    = topo;

    D.options = struct();
    D.options.join        = join;
    D.options.miter_limit = p.Results.miter_limit;
    D.options.corner_tol  = corner_tol;
    D.options.tip         = tipMode;
end


% =========================================================================
% local channel construction
% =========================================================================

function [Poly, G] = build_local_channel_polygon(Pmid, w, varargin)
%BUILD_LOCAL_CHANNEL_POLYGON Build a local pencil-like polygon around Pmid.
%
% Supports a general polyline, but is mainly intended for the current
% short-crack case. The polygon is built from offset chains on both sides
% of the polyline and closed at the mouth and tip.

    p = inputParser;
    addParameter(p, 'join',       'miter', @(s) ischar(s) || isstring(s));
    addParameter(p, 'corner_tol', 1e-10,   @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'tip',        'point', @(s) ischar(s) || isstring(s));
    parse(p, varargin{:});

    join        = tolower_safe(p.Results.join);
    tipMode     = tolower_safe(p.Results.tip);
    corner_tol = p.Results.corner_tol;

    N = size(Pmid,1);
    if N < 2
        error('build_domain_hole_pencil_polyline:BadPmid', ...
            'Pmid must contain at least two points.');
    end

    % segment tangents and left normals
    tseg = zeros(N-1,2);
    nseg = zeros(N-1,2);

    for i = 1:N-1
        d = Pmid(i+1,:) - Pmid(i,:);
        L = norm(d);
        if L <= corner_tol
            error('build_domain_hole_pencil_polyline:DegenerateSegment', ...
                'Degenerate segment detected in Pmid.');
        end
        t = d / L;
        n = [-t(2), t(1)];  % left normal

        tseg(i,:) = t;
        nseg(i,:) = n;
    end

    upper = zeros(N,2);
    lower = zeros(N,2);

    % mouth point: exact attachment to the hole boundary
    upper(1,:) = Pmid(1,:);
    lower(1,:) = Pmid(1,:);
    
    % interior vertices
    for i = 2:N-1
        s = (i-1)/(N-1);   % 0 at mouth, 1 at tip
        wi = w * s;

        switch join
            case 'miter'
                [xu, xl] = miter_join(Pmid(i,:), ...
                    tseg(i-1,:), nseg(i-1,:), ...
                    tseg(i,:),   nseg(i,:), wi, corner_tol);

                upper(i,:) = xu;
                lower(i,:) = xl;

            case 'bevel'
                navg = nseg(i-1,:) + nseg(i,:);
                if norm(navg) <= corner_tol
                    navg = nseg(i,:);
                end
                navg = navg / norm(navg);

                upper(i,:) = Pmid(i,:) + wi*navg;
                lower(i,:) = Pmid(i,:) - wi*navg;

            otherwise
                error('build_domain_hole_pencil_polyline:UnknownJoin', ...
                    'Unknown join option "%s".', join);
        end
    end

    % tip end: full width
    ntip = nseg(end,:);
    ttip = tseg(end,:);

    upper(end,:) = Pmid(end,:) + w*ntip;
    lower(end,:) = Pmid(end,:) - w*ntip;

    switch tipMode
        case 'flat'
            Poly = [upper; flipud(lower)];

        case 'point'
            tipPt = Pmid(end,:) + w*ttip;
            Poly = [ ...
                upper(1:end-1,:);
                upper(end,:);
                tipPt;
                lower(end,:);
                flipud(lower(1:end-1,:))];

        otherwise
            error('build_domain_hole_pencil_polyline:UnknownTipMode', ...
                'Unknown tip option "%s".', tipMode);
    end
    % remove accidental consecutive duplicates
    Poly = remove_consecutive_duplicates(Poly, corner_tol);

    G = struct();
    G.upper = upper;
    G.lower = lower;
    G.tseg  = tseg;
    G.nseg  = nseg;
    G.tipMode = tipMode;
end


function [xu, xl] = miter_join(xc, t1, n1, t2, n2, w, tol)
%MITER_JOIN Compute left/right offset points at a polyline vertex by
% intersecting adjacent offset lines.

    % upper side = +w*n
    [okU, xu] = line_intersection( ...
        xc + w*n1, t1, ...
        xc + w*n2, t2, tol);

    % lower side = -w*n
    [okL, xl] = line_intersection( ...
        xc - w*n1, t1, ...
        xc - w*n2, t2, tol);

    % fallback if nearly parallel
    if ~okU
        nu = n1 + n2;
        if norm(nu) <= tol
            nu = n2;
        end
        nu = nu / norm(nu);
        xu = xc + w*nu;
    end

    if ~okL
        nl = n1 + n2;
        if norm(nl) <= tol
            nl = n2;
        end
        nl = nl / norm(nl);
        xl = xc - w*nl;
    end
end


function [ok, x] = line_intersection(p1, d1, p2, d2, tol)
%LINE_INTERSECTION Intersect two parametric lines:
%   p1 + a*d1 = p2 + b*d2

    A = [d1(:), -d2(:)];
    b = (p2(:) - p1(:));

    if abs(det(A)) <= tol
        ok = false;
        x = [NaN, NaN];
        return;
    end

    ab = A \ b;
    a = ab(1);

    x = (p1(:) + a*d1(:)).';
    ok = true;
end


function P2 = remove_consecutive_duplicates(P, tol)
%REMOVE_CONSECUTIVE_DUPLICATES Remove repeated consecutive vertices.

    if isempty(P)
        P2 = P;
        return;
    end

    keep = true(size(P,1),1);
    for i = 2:size(P,1)
        if norm(P(i,:) - P(i-1,:), inf) <= tol
            keep(i) = false;
        end
    end
    P2 = P(keep,:);
end


% =========================================================================
% helpers
% =========================================================================

function holes = normalize_holes(holes)
%NORMALIZE_HOLES Normalize hole input to a cell array of structs.

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
%DETECT_HOLE_TOUCH Detect whether point x0 lies on/near a hole boundary.

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
%POLYGON_BBOX Bounding box of polygon vertices.
% Returns [xmin xmax ymin ymax].

    bbox = [min(P(:,1)), max(P(:,1)), min(P(:,2)), max(P(:,2))];
end


function tf = point_in_box(x, box, tol)
%POINT_IN_BOX Test if point x lies inside box = [xmin xmax ymin ymax].

    tf = (x(1) >= box(1) - tol) && (x(1) <= box(2) + tol) && ...
         (x(2) >= box(3) - tol) && (x(2) <= box(4) + tol);
end


function A = signed_polygon_area(P)
%SIGNED_POLYGON_AREA Signed area of polygon P.
% Positive for CCW orientation.

    x = P(:,1);
    y = P(:,2);

    x2 = [x(2:end); x(1)];
    y2 = [y(2:end); y(1)];

    A = 0.5 * sum(x .* y2 - x2 .* y);
end


function d = point_polygon_distance(x0, V)
%POINT_POLYGON_DISTANCE Minimum distance from point x0 to polygon edges.

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
%POINT_SEGMENT_DISTANCE Distance from point P to segment AB.

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


function must_field(S, field)
%MUST_FIELD Check struct field existence.

    if ~isfield(S, field) || isempty(S.(field))
        error('build_domain_hole_pencil_polyline:MissingHoleField', ...
            'Required field "%s" is missing or empty.', field);
    end
end

function s = tolower_safe(x)
%TOLOWER_SAFE Convert char/string to lowercase without relying on lower().
%
% This avoids rare name-resolution issues with lower in some MATLAB setups.

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