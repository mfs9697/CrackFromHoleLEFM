function [Vapp, G] = build_appended_hole_loop(hole, A, n_in, t_hat, theta, a0, eps, varargin)
%BUILD_APPENDED_HOLE_LOOP
% Build one closed loop for a circular hole with an appended short crack.
%
% Geometry idea:
%   - remove a small hole arc around the initiation point A
%   - replace it by an appendix consisting of:
%       lower face + rounded tip cap + upper face
%
% Inputs:
%   hole   : struct with fields
%              .type   = 'circle'
%              .center = [xc yc]
%              .r      = radius
%   A      : initiation point on the hole boundary [1x2]
%   n_in   : inward unit normal at A
%   t_hat  : tangent unit vector at A
%   theta  : crack angle relative to inward normal, radians
%   a0     : appendix length measured along crack direction
%   eps    : small mouth half-shift
%            default interpretation is arc length along the circle
%
% Name-value options:
%   'epsMode'     : 'arclength' (default) or 'angle'
%   'nArc'        : number of points on retained hole arc     (default 160)
%   'nFace'       : number of points on each appendix face    (default 8)
%   'nTip'        : number of points on tip cap               (default 21)
%   'tipRadius'   : tip-cap radius (default = auto)
%   'orientation' : 'cw' (default) or 'ccw'
%
% Outputs:
%   Vapp   : [N x 2] closed polygon loop of the appended hole
%            (last point is NOT repeated)
%   G      : debug struct
%
% Notes:
%   - first draft: circular hole only
%   - first draft: straight appendix faces
%   - suitable for replacing the touched hole loop in Stage II geometry

    % -------------------- parse options --------------------
    ip = inputParser;
    addParameter(ip, 'epsMode', 'arclength', @(s)ischar(s)||isstring(s));
    addParameter(ip, 'nArc', 160, @(x)isnumeric(x)&&isscalar(x)&&x>=8);
    addParameter(ip, 'nFace', 8, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
    addParameter(ip, 'nTip', 21, @(x)isnumeric(x)&&isscalar(x)&&x>=5);
    addParameter(ip, 'tipRadius', [], @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(ip, 'orientation', 'cw', @(s)ischar(s)||isstring(s));
    parse(ip, varargin{:});

    epsMode     = lower_safe(ip.Results.epsMode);
    nArc        = ip.Results.nArc;
    nFace       = ip.Results.nFace;
    nTip        = ip.Results.nTip;
    tipRadiusIn = ip.Results.tipRadius;
    orient      = lower_safe(ip.Results.orientation);

    % -------------------- checks --------------------
    must_field(hole, 'type');
    must_field(hole, 'center');
    must_field(hole, 'r');

    if ~strcmpi(strtrim(hole.type), 'circle')
        error('build_appended_hole_loop:UnsupportedHoleType', ...
            'This first draft supports only circular holes.');
    end

    c = hole.center(:).';
    r = hole.r;

    A     = A(:).';
    n_in  = normalize_row(n_in);
    t_hat = normalize_row(t_hat);

    if ~(isnumeric(theta) && isscalar(theta) && isfinite(theta))
        error('build_appended_hole_loop:BadTheta', 'theta must be a finite scalar.');
    end
    if ~(isnumeric(a0) && isscalar(a0) && a0 > 0)
        error('build_appended_hole_loop:BadA0', 'a0 must be a positive scalar.');
    end
    if ~(isnumeric(eps) && isscalar(eps) && eps > 0)
        error('build_appended_hole_loop:BadEps', 'eps must be a positive scalar.');
    end

    % -------------------- project A to the circle --------------------
    rc = A - c;
    nr = norm(rc);
    if nr <= 0
        error('build_appended_hole_loop:BadA', ...
            'Initiation point A cannot coincide with the hole center.');
    end
    A0 = c + r * rc / nr;

    phiA = atan2(A0(2)-c(2), A0(1)-c(1));

    % Circle tangent for increasing polar angle
    tcirc = [-sin(phiA), cos(phiA)];

    % We want +psi direction to follow +t_hat along the hole
    if dot(tcirc, t_hat) >= 0
        sgnPsi = +1;
    else
        sgnPsi = -1;
    end

    % -------------------- mouth half-shift --------------------
    switch epsMode
        case 'arclength'
            dphi = eps / r;
        case 'angle'
            dphi = eps;
        otherwise
            error('build_appended_hole_loop:BadEpsMode', ...
                'epsMode must be ''arclength'' or ''angle''.');
    end

    if dphi <= 0 || dphi >= pi/2
        error('build_appended_hole_loop:BadDphi', ...
            'Computed half-angle dphi must satisfy 0 < dphi < pi/2.');
    end

    % Raw shifted points: one on each side of A along the circle
    phi1 = phiA + sgnPsi*dphi;
    phi2 = phiA - sgnPsi*dphi;

    M1 = c + r*[cos(phi1), sin(phi1)];
    M2 = c + r*[cos(phi2), sin(phi2)];

    % -------------------- crack direction and face normal --------------------
    % Same convention as geom_hole_shortcrack.m:
    % e_dir = cos(theta)*n_in + sin(theta)*t_hat
    e_dir = cos(theta)*n_in + sin(theta)*t_hat;
    e_dir = normalize_row(e_dir);

    n_face = [-e_dir(2), e_dir(1)];

    % Classify which mouth point is "upper" / "lower" relative to n_face
    if dot(M1 - A0, n_face) >= dot(M2 - A0, n_face)
        Mup = M1;
        Mlo = M2;
        phi_up = phi1;
        phi_lo = phi2;
    else
        Mup = M2;
        Mlo = M1;
        phi_up = phi2;
        phi_lo = phi1;
    end

    % -------------------- tip cap geometry --------------------
    mouthWidth = norm(Mup - Mlo);

    if isempty(tipRadiusIn)
        rtip = min(0.25*a0, 0.50*mouthWidth);
    else
        rtip = tipRadiusIn;
    end

    if ~(rtip > 0 && rtip < a0)
        error('build_appended_hole_loop:BadTipRadius', ...
            'tipRadius must satisfy 0 < tipRadius < a0.');
    end

    xtip_front = A0 + a0*e_dir;
    xcap       = xtip_front - rtip*e_dir;

    Pup_end = xcap + rtip*n_face;
    Plo_end = xcap - rtip*n_face;

    % -------------------- appendix faces --------------------
    su = linspace(0,1,nFace).';
    upper = Mup + su.*(Pup_end - Mup);

    sl = linspace(0,1,nFace).';
    lower = Mlo + sl.*(Plo_end - Mlo);

    % -------------------- retained hole arc --------------------
    % Build the long arc from Mup to Mlo that excludes the removed mouth arc around A.
    phiArc = long_circle_arc_excluding_point(phi_up, phi_lo, phiA, nArc);

    arcMain = c + r*[cos(phiArc), sin(phiArc)];

    % -------------------- tip cap --------------------
    % Arc from lower-face end to upper-face end around the front side
    alpha_lo = atan2(Plo_end(2)-xcap(2), Plo_end(1)-xcap(1));
    alpha_up = atan2(Pup_end(2)-xcap(2), Pup_end(1)-xcap(1));

    alphaCap = front_cap_angles(alpha_lo, alpha_up, e_dir, nTip);

    tipcap = xcap + rtip*[cos(alphaCap), sin(alphaCap)];

    % -------------------- assemble loop --------------------
    % Start at Mup, go around the retained circle to Mlo,
    % then along lower face to tip, around cap, and back along upper face.
    Vapp = [
        arcMain;
        lower(2:end,:);
        tipcap(2:end-1,:);
        flipud(upper(1:end-1,:))
    ];

    % Remove accidental duplicates
    Vapp = remove_consecutive_duplicates(Vapp, 1e-12);
    Vapp = remove_last_if_equal_first(Vapp, 1e-12);

    % Enforce requested orientation
    Apoly = signed_polygon_area(Vapp);
    switch orient
        case 'cw'
            if Apoly > 0
                Vapp = flipud(Vapp);
            end
        case 'ccw'
            if Apoly < 0
                Vapp = flipud(Vapp);
            end
        otherwise
            error('build_appended_hole_loop:BadOrientation', ...
                'orientation must be ''cw'' or ''ccw''.');
    end

    % -------------------- debug output --------------------
    G = struct();
    G.center       = c;
    G.r            = r;
    G.A_input      = A;
    G.A            = A0;
    G.phiA         = phiA;
    G.dphi         = dphi;
    G.M1           = M1;
    G.M2           = M2;
    G.Mup          = Mup;
    G.Mlo          = Mlo;
    G.phi_up       = phi_up;
    G.phi_lo       = phi_lo;
    G.n_in         = n_in;
    G.t_hat        = t_hat;
    G.theta        = theta;
    G.e_dir        = e_dir;
    G.n_face       = n_face;
    G.a0           = a0;
    G.xtip_front   = xtip_front;
    G.xcap         = xcap;
    G.rtip         = rtip;
    G.arcMain      = arcMain;
    G.upper        = upper;
    G.lower        = lower;
    G.tipcap       = tipcap;
    G.orientation  = orient;
end


% =========================================================================
% helpers
% =========================================================================

function phiArc = long_circle_arc_excluding_point(phi_start, phi_end, phi_excl, nArc)
% Return the arc from phi_start to phi_end that excludes phi_excl.
% Two possibilities exist; choose the one that does NOT contain phi_excl.

    phi_start = wrap_2pi(phi_start);
    phi_end   = wrap_2pi(phi_end);
    phi_excl  = wrap_2pi(phi_excl);

    % Candidate 1: increasing arc from start to end
    phi_end_inc = phi_end;
    while phi_end_inc < phi_start
        phi_end_inc = phi_end_inc + 2*pi;
    end
    phi_excl_inc = phi_excl;
    while phi_excl_inc < phi_start
        phi_excl_inc = phi_excl_inc + 2*pi;
    end

    excl_on_inc = (phi_excl_inc >= phi_start) && (phi_excl_inc <= phi_end_inc);

    if excl_on_inc
        % Use decreasing arc instead
        phi_end_dec = phi_end;
        while phi_end_dec > phi_start
            phi_end_dec = phi_end_dec - 2*pi;
        end
        phiArc = linspace(phi_start, phi_end_dec, nArc).';
    else
        % Increasing arc excludes phi_excl
        phiArc = linspace(phi_start, phi_end_inc, nArc).';
    end
end


function alphaCap = front_cap_angles(alpha_lo, alpha_up, e_dir, nTip)
% Choose the cap arc from lower to upper around the front side of the slot.

    % Candidate 1: increasing arc
    a1 = alpha_up;
    while a1 < alpha_lo
        a1 = a1 + 2*pi;
    end
    cand1 = linspace(alpha_lo, a1, nTip).';

    % Candidate 2: decreasing arc
    a2 = alpha_up;
    while a2 > alpha_lo
        a2 = a2 - 2*pi;
    end
    cand2 = linspace(alpha_lo, a2, nTip).';

    % Pick the one whose midpoint lies more in the forward e_dir direction
    mid1 = [cos(cand1(round(end/2))), sin(cand1(round(end/2)))];
    mid2 = [cos(cand2(round(end/2))), sin(cand2(round(end/2)))];

    if dot(mid1, e_dir) >= dot(mid2, e_dir)
        alphaCap = cand1;
    else
        alphaCap = cand2;
    end
end


function v = normalize_row(v)
    v = v(:).';
    nv = norm(v);
    if nv <= 0
        error('build_appended_hole_loop:ZeroVector', ...
            'Cannot normalize a zero vector.');
    end
    v = v / nv;
end


function s = lower_safe(x)
    if isstring(x)
        x = char(x);
    end
    if ~ischar(x)
        error('build_appended_hole_loop:BadStringInput', ...
            'Expected char or string input.');
    end
    s = x;
    mask = (s >= 'A') & (s <= 'Z');
    s(mask) = char(double(s(mask)) + ('a' - 'A'));
end


function x = wrap_2pi(x)
    x = mod(x, 2*pi);
    if x < 0
        x = x + 2*pi;
    end
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


function P2 = remove_last_if_equal_first(P, tol)
    if size(P,1) >= 2 && norm(P(end,:) - P(1,:), inf) <= tol
        P2 = P(1:end-1,:);
    else
        P2 = P;
    end
end


function A = signed_polygon_area(P)
    x = P(:,1);
    y = P(:,2);
    x2 = [x(2:end); x(1)];
    y2 = [y(2:end); y(1)];
    A = 0.5 * sum(x .* y2 - x2 .* y);
end


function must_field(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('build_appended_hole_loop:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end