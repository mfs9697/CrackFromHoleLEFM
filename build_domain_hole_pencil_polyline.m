function D = build_domain_hole_pencil_polyline(Pmid, A, B, holes, w, varargin)
%BUILD_DOMAIN_HOLE_PENCIL_POLYLINE
% Build geometric description for a rectangular plate with:
%   - one or more holes,
%   - and, if the crack starts on a circular hole, a merged appended-hole loop
%     instead of a separate pencil channel.
%
%   D = build_domain_hole_pencil_polyline(Pmid, A, B, holes, w)
%
% Inputs:
%   Pmid   : (N x 2) crack midline vertices, ordered mouth -> ... -> tip
%   A, B   : plate dimensions, domain x in [0,A], y in [-B,B]
%   holes  : cell array of hole structs, e.g.
%              { struct('type','circle','center',[xc yc],'r',R,'npoly',N), ... }
%            or empty {}
%   w      : legacy channel half-width; in merged mode, used to set the
%            mouth half-shift eps unless overridden
%
% Name-value options:
%   'join'        : 'miter' (default) or 'bevel'   [legacy fallback only]
%   'miter_limit' : positive scalar (default 6)    [legacy fallback only]
%   'corner_tol'  : positive scalar (default 1e-10)
%   'tip'         : 'point' (default) or 'flat'    [legacy fallback only]
%
%   'mode'        : 'auto' (default), 'merged', or 'legacy'
%                   auto   -> merged if crack starts on a supported hole,
%                             otherwise legacy separate-channel geometry
%                   merged -> force appended-hole construction
%                   legacy -> force old separate channel construction
%
%   'mouth_eps'   : mouth half-shift used for appended hole construction
%                   default = w
%   'epsMode'     : 'arclength' (default) or 'angle'
%   'nArc'        : retained-hole arc resolution (default 160)
%   'nFace'       : appendix face resolution     (default 8)
%   'nTip'        : appendix tip-cap resolution  (default 21)
%   'tipRadius'   : tip-cap radius for appended hole (default = auto)
%   'orientation' : 'cw' (default) or 'ccw' for appended inner loops
%
% Output:
%   D      : geometry-description struct with fields
%       .outerPoly       rectangle polygon [4 x 2], CCW
%       .holeLoops       cell array of inner polygons
%       .channelPoly     [] in merged mode; legacy local channel polygon otherwise
%       .Pmid            crack midline
%       .A, .B, .w       copied inputs
%       .holes           normalized hole list
%       .channelGeom     debug struct
%       .topology        auxiliary geometric/topological info
%
% Notes:
%   - In merged mode, the touched circular hole is replaced by a single
%     appended-hole contour, and channelPoly is empty.
%   - In legacy mode, the function reproduces the old behavior:
%         outer rectangle - holes - separate channel
%
% Companion note:
%   mesh_hole_pencil_domain.m should be adjusted to allow D.channelPoly = []
%   and mesh only outerPoly minus holeLoops in merged mode.

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

    addParameter(p, 'mode',        'auto',  @(s) ischar(s) || isstring(s));
    addParameter(p, 'mouth_eps',   [],      @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'epsMode',     'arclength', @(s) ischar(s) || isstring(s));
    addParameter(p, 'nArc',        160,     @(x) isnumeric(x) && isscalar(x) && x >= 8);
    addParameter(p, 'nFace',       8,       @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'nTip',        21,      @(x) isnumeric(x) && isscalar(x) && x >= 5);
    addParameter(p, 'tipRadius',   [],      @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'orientation', 'cw',    @(s) ischar(s) || isstring(s));

    parse(p, varargin{:});

    join         = tolower_safe(p.Results.join);
    tipMode      = tolower_safe(p.Results.tip);
    corner_tol   = p.Results.corner_tol;
    modeReq      = tolower_safe(p.Results.mode);
    epsMode      = tolower_safe(p.Results.epsMode);
    orient       = tolower_safe(p.Results.orientation);

    mouth_eps    = p.Results.mouth_eps;
    nArc         = p.Results.nArc;
    nFace        = p.Results.nFace;
    nTip         = p.Results.nTip;
    tipRadius    = p.Results.tipRadius;

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
    % hole loops (initial, before possible replacement)
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

    if ~isempty(holes)
        [itouch, dtouch] = detect_hole_touch(Pmid(1,:), holes);
        topo.channelTouchesHole = ~isempty(itouch);
        topo.touchingHoleIndex  = itouch;
        topo.touchDistance      = dtouch;
    else
        itouch = [];
    end

    % ------------------------------------------------------------
    % choose geometry mode
    % ------------------------------------------------------------
    switch modeReq
        case 'legacy'
            useMerged = false;

        case 'merged'
            useMerged = true;

        case 'auto'
            useMerged = topo.channelTouchesHole && ...
                        ~isempty(itouch) && ...
                        strcmpi(strtrim(holes{itouch}.type), 'circle');

        otherwise
            error('build_domain_hole_pencil_polyline:BadMode', ...
                'Unknown mode "%s". Use auto, merged, or legacy.', modeReq);
    end

    % ------------------------------------------------------------
    % build merged appended-hole OR legacy separate channel
    % ------------------------------------------------------------
    channelPoly = [];
    channelGeom = struct();

    if useMerged
        if isempty(itouch)
            error('build_domain_hole_pencil_polyline:NoTouchedHole', ...
                'Merged mode requested, but Pmid(1,:) was not detected on a hole.');
        end

        hk = holes{itouch};
        if ~strcmpi(strtrim(hk.type), 'circle')
            error('build_domain_hole_pencil_polyline:UnsupportedMergedHole', ...
                'Merged mode currently supports only circular holes.');
        end

        if exist('build_appended_hole_loop', 'file') ~= 2
            error('build_domain_hole_pencil_polyline:MissingHelper', ...
                'build_appended_hole_loop.m is required on the MATLAB path.');
        end

        % Infer local frame and crack angle from Pmid and the touched circular hole
        [A0, n_in, t_hat, theta, a0, Gframe] = infer_local_frame_from_pmid(hk, Pmid, corner_tol);

        [Vapp, Gapp] = build_appended_hole_loop(hk, A0, n_in, t_hat, theta, a0, mouth_eps, ...
            'epsMode', epsMode, ...
            'nArc', nArc, ...
            'nFace', nFace, ...
            'nTip', nTip, ...
            'tipRadius', tipRadius, ...
            'orientation', orient);

        holeLoops{itouch} = Vapp;

        channelGeom.mode   = 'merged_appended_hole';
        channelGeom.frame  = Gframe;
        channelGeom.append = Gapp;

        topo.mode = 'merged_appended_hole';
        topo.appendedHoleIndex = itouch;
        topo.appendedHoleArea  = abs(signed_polygon_area(Vapp));

        % no separate channel in merged mode
        channelPoly = [];

    else
        % ---------- legacy separate-channel behavior ----------
        [channelPoly, Gc] = build_local_channel_polygon(Pmid, w, ...
            'join', join, ...
            'corner_tol', corner_tol, ...
            'tip', tipMode);

        if signed_polygon_area(channelPoly) < 0
            channelPoly = flipud(channelPoly);
        end

        channelGeom.mode = 'legacy_separate_channel';
        channelGeom.legacy = Gc;

        topo.mode = 'legacy_separate_channel';
        topo.channelArea = polyarea(channelPoly(:,1), channelPoly(:,2));
    end

    % ------------------------------------------------------------
    % diagnostics / flags
    % ------------------------------------------------------------
    topo.bbox.outer = polygon_bbox(outerPoly);
    topo.bbox.holes = cell(size(holeLoops));
    for k = 1:numel(holeLoops)
        topo.bbox.holes{k} = polygon_bbox(holeLoops{k});
    end

    if ~isempty(channelPoly)
        topo.bbox.channel = polygon_bbox(channelPoly);
    else
        topo.bbox.channel = [];
    end

    topo.flags = struct();
    topo.flags.tipInsidePlate = point_in_box(Pmid(end,:), [0, A, -B, B], 1e-12);
    topo.flags.midlineInsidePlate = all(Pmid(:,1) >= -1e-12 & Pmid(:,1) <= A + 1e-12 & ...
                                        Pmid(:,2) >= -B - 1e-12 & Pmid(:,2) <= B + 1e-12);

    topo.flags.startOnHole = topo.channelTouchesHole;

    if ~isempty(channelPoly)
        topo.flags.channelInsidePlateBox = all(channelPoly(:,1) >= -1e-12 & channelPoly(:,1) <= A + 1e-12 & ...
                                               channelPoly(:,2) >= -B - 1e-12 & channelPoly(:,2) <= B + 1e-12);
    else
        topo.flags.channelInsidePlateBox = [];
    end

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    D = struct();

    D.outerPoly   = outerPoly;
    D.holeLoops   = holeLoops;
    D.channelPoly = channelPoly;   % [] in merged mode

    D.Pmid        = Pmid;
    D.A           = A;
    D.B           = B;
    D.w           = w;
    D.holes       = holes;

    D.channelGeom = channelGeom;
    D.topology    = topo;

    D.options = struct();
    D.options.join        = join;
    D.options.miter_limit = p.Results.miter_limit;
    D.options.corner_tol  = corner_tol;
    D.options.tip         = tipMode;
    D.options.mode        = modeReq;
    D.options.mouth_eps   = mouth_eps;
    D.options.epsMode     = epsMode;
    D.options.nArc        = nArc;
    D.options.nFace       = nFace;
    D.options.nTip        = nTip;
    D.options.tipRadius   = tipRadius;
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
% legacy local channel construction
% =========================================================================

function [Poly, G] = build_local_channel_polygon(Pmid, w, varargin)
% Build the old separate pencil-like polygon around Pmid.

    p = inputParser;
    addParameter(p, 'join',       'miter', @(s) ischar(s) || isstring(s));
    addParameter(p, 'corner_tol', 1e-10,   @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'tip',        'point', @(s) ischar(s) || isstring(s));
    parse(p, varargin{:});

    join       = tolower_safe(p.Results.join);
    tipMode    = tolower_safe(p.Results.tip);
    corner_tol = p.Results.corner_tol;

    N = size(Pmid,1);
    if N < 2
        error('build_domain_hole_pencil_polyline:BadPmid', ...
            'Pmid must contain at least two points.');
    end

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
        n = [-t(2), t(1)];

        tseg(i,:) = t;
        nseg(i,:) = n;
    end

    upper = zeros(N,2);
    lower = zeros(N,2);

    upper(1,:) = Pmid(1,:);
    lower(1,:) = Pmid(1,:);

    for i = 2:N-1
        s = (i-1)/(N-1);
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
                flipud(lower(1:end-1,:)) ];

        otherwise
            error('build_domain_hole_pencil_polyline:UnknownTipMode', ...
                'Unknown tip option "%s".', tipMode);
    end

    Poly = remove_consecutive_duplicates(Poly, corner_tol);

    G = struct();
    G.upper = upper;
    G.lower = lower;
    G.tseg  = tseg;
    G.nseg  = nseg;
    G.tipMode = tipMode;
end


function [xu, xl] = miter_join(xc, t1, n1, t2, n2, w, tol)
    [okU, xu] = line_intersection( ...
        xc + w*n1, t1, ...
        xc + w*n2, t2, tol);

    [okL, xl] = line_intersection( ...
        xc - w*n1, t1, ...
        xc - w*n2, t2, tol);

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


function P2 = remove_consecutive_duplicates(P, tol)
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