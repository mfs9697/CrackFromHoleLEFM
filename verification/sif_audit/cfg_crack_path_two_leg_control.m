function C = cfg_crack_path_two_leg_control(theta2_deg, varargin)
%CFG_CRACK_PATH_TWO_LEG_CONTROL
% Reproducible two-leg elastic crack control case for the SIF audit.
%
% The baseline geometry/material/mesh values are taken from the current
% non-perforated Crack-Path configuration (cfg_nonperforated.m):
%   A = B = 0.10 m
%   a = 0.02 m
%   theta1 = -25 deg
%   delta = 0.08*a
%   E = 4e3 MPa, nu = 0.30, plane strain
%   ncoh = 40, hgrad = 1.15, B/hmax = 40
%
% IMPORTANT:
%   This file defines a CONTROL case for code-to-code verification.  The
%   second leg is treated as a traction-free elastic crack segment.  It is
%   not yet claimed to reproduce the older published two-segment benchmark.
%
% Usage:
%   C = cfg_crack_path_two_leg_control();
%   C = cfg_crack_path_two_leg_control(2.0);
%   C = cfg_crack_path_two_leg_control(2.0, 'ncoh', 80);
%
% theta2_deg is the GLOBAL angle of the second crack leg.

    if nargin < 1 || isempty(theta2_deg)
        theta2_deg = 2.0;
    end

    ip = inputParser;
    addParameter(ip, 'A', 0.10, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'B', 0.10, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'crackLength', 0.02, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'theta1_deg', -25.0, @(x)isnumeric(x) && isscalar(x));
    addParameter(ip, 'delta_over_a', 0.08, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'ncoh', 40, @(x)isnumeric(x) && isscalar(x) && x>=2);
    addParameter(ip, 'hgrad', 1.15, @(x)isnumeric(x) && isscalar(x) && x>1);
    addParameter(ip, 'hmax_ratio', 40, @(x)isnumeric(x) && isscalar(x) && x>1);
    addParameter(ip, 'chw_factor', 1/16, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'E', 4e3, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'nu', 0.30, @(x)isnumeric(x) && isscalar(x) && x>0 && x<0.5);
    addParameter(ip, 'sigma0', 1.0, @(x)isnumeric(x) && isscalar(x));
    addParameter(ip, 'plotGeom', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'plotMesh', false, @(x)islogical(x) || isnumeric(x));
    parse(ip, varargin{:});
    S = ip.Results;

    theta1 = deg2rad(S.theta1_deg);
    theta2 = deg2rad(theta2_deg);
    a      = S.crackLength;
    delta  = S.delta_over_a * a;

    P0 = [0.0, 0.0];
    P1 = P0 + a      * [cos(theta1), sin(theta1)];
    P2 = P1 + delta  * [cos(theta2), sin(theta2)];

    C = struct();

    C.A = S.A;
    C.B = S.B;

    C.P0   = P0;
    C.Pmid = [P0; P1; P2];
    C.L    = [a; delta];
    C.theta = [theta1; theta2];

    C.a = a;
    C.delta = delta;
    C.theta1 = theta1;
    C.theta2 = theta2;

    C.ncoh       = round(S.ncoh);
    C.hgrad      = S.hgrad;
    C.hmax_ratio = S.hmax_ratio;

    hlast = delta / C.ncoh;
    C.chw = hlast * S.chw_factor;

    C.join        = 'miter';
    C.miter_limit = 6;
    C.corner_tol  = 1e-10;
    C.tip         = 'point';

    C.E2  = S.E;
    C.E   = S.E;
    C.nu  = S.nu;
    C.ps  = 1;       % Crack-Path LEFM branch is plane strain.
    C.G12 = C.E2/(2*(1+C.nu));

    coef = C.E2 / ((1+C.nu)*(1-2*C.nu));
    C.Dmat = coef * [ ...
        1-C.nu, C.nu, 0; ...
        C.nu, 1-C.nu, 0; ...
        0, 0, (1-2*C.nu)/2 ];

    C.eps1 = 1e-9;
    C.sigma0 = S.sigma0;

    C.plotGeom = logical(S.plotGeom);
    C.plotMesh = logical(S.plotMesh);

    C.audit = struct();
    C.audit.source_geometry = 'mfs9697/Crack-Path:cfg_nonperforated.m';
    C.audit.source_mesh     = 'mfs9697/Crack-Path:geom_pencil.m';
    C.audit.role = 'code-to-code SIF control; not yet the published historical benchmark';
end
