function [Vapp, G] = build_appended_hole_loop(hole, A, n_in, t_hat, theta, a0, eps, varargin)
%BUILD_APPENDED_HOLE_LOOP
% Build one closed loop for a circular hole with a SHARP appended pencil.
%
% Geometry idea:
%   - remove a small hole arc around the initiation point A
%   - replace it by two straight appendix faces meeting at one sharp tip
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
%   'orientation' : 'cw' (default) or 'ccw'
%
% Outputs:
%   Vapp   : [N x 2] closed polygon loop of the appended hole
%            (last point is NOT repeated)
%   G      : debug struct
%
% Notes:
%   - circular hole only
%   - sharp pencil: no rounded tip
%   - intended for T3 meshing with local Hedge/Hvertex refinement near tip

    % -------------------- parse options --------------------
    ip = inputParser;
    addParameter(ip, 'epsMode', 'arclength', @(s)ischar(s)||isstring(s));
    addParameter(ip, 'nArc', 160, @(x)isnumeric(x)&&isscalar(x)&&x>=8);
    addParameter(ip, 'orientation', 'cw', @(s)ischar(s)||isstring(s));
    parse(ip, varargin{:});

    epsMode = lower_safe(ip.Results.epsMode);
    nArc    = ip.Results.nArc;
    orient  = lower_safe(ip.Results.orientation);

    % -------------------- checks --------------------
    must_field(hole, 'type');
    must_field(hole, 'center');
    must_field(hole, 'r');

    if ~strcmpi(strtrim(hole.type), 'circle')
        error('build_appended_hole_loop:UnsupportedHoleType', ...
            'This version supports only circular holes.');
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

    % Two mouth points on the circle around A
    phi1 = phiA + sgnPsi*dphi;
    phi2 = phiA - sgnPsi*dphi;

    M1 = c + r*[cos(phi1), sin(phi1)];
    M2 = c + r*[cos(phi2), sin(phi2)];

    % -------------------- crack direction --------------------
    % Same convention as geom_hole_shortcrack.m:
    % e_dir = cos(theta)*n_in + sin(theta)*t_hat
    e_dir = cos(theta)*n_in + sin(theta)*t_hat;
    e_dir = normalize_row(e_dir);

    n_face = [-e_dir(2), e_dir(1)];

    % Classify upper/lower mouth points relative to the crack-face normal
    if dot(M1 - A0, n_face) >= dot(M2 - A0, n_face)
        Mup    = M1;
        Mlo    = M2;
        phi_up = phi1;
        phi_lo = phi2;
    else
        Mup    = M2;
        Mlo    = M1;
        phi_up = phi2;
        phi_lo = phi1;
    end

    % -------------------- sharp tip --------------------
    xtip = A0 + a0*e_dir;

    % -------------------- retained hole arc --------------------
    % Build the long arc from Mup to Mlo that excludes the removed mouth arc around A
    phiArc = long_circle_arc_excluding_point(phi_up, phi_lo, phiA, nArc);
    arcMain = c + r*[cos(phiArc), sin(phiArc)];

    % -------------------- assemble sharp loop --------------------
    % Path:
    %   long hole arc: Mup -> ... -> Mlo
    %   then straight edge Mlo -> xtip
    %   and polygon closure gives xtip -> Mup
    Vapp = [
        arcMain;
        xtip
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
    G.center      = c;
    G.r           = r;
    G.A_input     = A;
    G.A           = A0;
    G.phiA        = phiA;
    G.dphi        = dphi;
    G.M1          = M1;
    G.M2          = M2;
    G.Mup         = Mup;
    G.Mlo         = Mlo;
    G.phi_up      = phi_up;
    G.phi_lo      = phi_lo;
    G.n_in        = n_in;
    G.t_hat       = t_hat;
    G.theta       = theta;
    G.e_dir       = e_dir;
    G.n_face      = n_face;
    G.a0          = a0;
    G.xtip        = xtip;
    G.arcMain     = arcMain;
    G.face_upper  = [Mup; xtip];
    G.face_lower  = [xtip; Mlo];
    G.orientation = orient;
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