function C = cfg_twoleg_historical()
%CFG_TWOLEG_HISTORICAL  Frozen reconstruction of the Crack-Path two-leg setup.
%
% Source provenance:
%   mfs9697/Crack-Path/cfg.m (historical working configuration)
%   IJF 2026 paper, Sections 2, 7, and 8.
%
% IMPORTANT:
%   The published CZM problem treats Leg 1 as the physical crack and Leg 2
%   as a short cohesive continuation.  For the present SIF verification
%   study we use exactly the same two-leg geometry but make BOTH legs
%   traction-free so that the SIFs are defined at V2.  This is a verification
%   benchmark derived from the historical geometry; it is not claimed to be
%   a published second-tip SIF result.

    cm = 0.01;

    % Plate: x in [0,A], y in [-B,B].
    C.A = 10*cm;
    C.B = 10*cm;

    % Historical two-leg geometry.
    C.a = 2*cm;
    C.delta = 0.06*C.a;

    C.theta1_deg = -20.0;
    C.theta2_deg =  20.0;
    C.theta1 = deg2rad(C.theta1_deg);
    C.theta2 = deg2rad(C.theta2_deg);

    C.V0 = [0,0];
    C.V1 = C.V0 + C.a*[cos(C.theta1),sin(C.theta1)];
    C.V2 = C.V1 + C.delta*[cos(C.theta2),sin(C.theta2)];
    C.Pmid = [C.V0;C.V1;C.V2];

    % Local frame at the second-leg tip.
    C.e1 = [cos(C.theta2);sin(C.theta2)];
    C.e2 = [-sin(C.theta2);cos(C.theta2)];
    C.R_gl = [C.e1,C.e2];
    C.R_loc = C.R_gl.';

    % Plane strain material; E is in MPa, geometry in metres.
    C.E = 4e3;
    C.nu = 0.3;
    C.ps = 1;
    coef = C.E/((1+C.nu)*(1-2*C.nu));
    C.Dmat = coef*[ ...
        1-C.nu, C.nu, 0; ...
        C.nu, 1-C.nu, 0; ...
        0,0,(1-2*C.nu)/2];

    % Historical mesh controls from Crack-Path/cfg.m.
    C.ncoh = 40;
    C.hgrad = 1.1;
    C.hmax = C.B/40;
    C.chw = (1/8)*(C.delta/C.ncoh);

    % Unit reference stress for a linear-elastic baseline. Results scale
    % linearly with this value.
    C.sigma0 = 1.0; % MPa

    % Historical point constraints.
    C.fix_points = [C.A,0; 0,C.B; 0,-C.B];

    % Verification sweeps.
    C.rI_over_delta = [0.20 0.30 0.40 0.50 0.60];
    C.edi_outer_over_delta = [0.30 0.40 0.50 0.60];
    C.edi_inner_fraction = [0.00 0.20 0.35];

    C.nthet = 240;
    C.eps_th = 1e-3;

    % Mesh-family controls. "control" uses the unperturbed PDE mesh;
    % asymmetric variants perturb only interior nodes in a local tip annulus.
    C.asymmetry_levels = [0 0.06 0.12 0.20];
    C.perturb_rmin = 0.10*C.delta;
    C.perturb_rmax = 0.85*C.delta;

    C.plotMesh = false;
end
