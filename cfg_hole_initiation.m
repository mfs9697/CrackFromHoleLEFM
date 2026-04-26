function C = cfg_hole_initiation()
%CFG_HOLE_INITIATION  Configuration for the two-stage hole-initiation study.
%
% Two-stage model:
%   Stage I  : elastic plate with hole, determine initiation point on hole boundary
%   Stage II : insert a short crack at that point, build sharp appended-hole geometry
%
% Output:
%   C   struct with geometry, material, loading, meshing, and solver settings

%% ================== GEOMETRY ==================
% Rectangular plate:
%   x in [0, A], y in [-B, B]
C.A = 0.30;
C.B = 0.10;

% Circular hole
C.hole.type   = 'circle';
C.hole.center = [0.150, -0.060];
C.hole.r      = 0.030;
C.hole.npoly  = 240;

% Current project uses one hole
C.holes = {C.hole};

%% ================== MATERIAL ==================
C.E  = 210e3;
C.nu = 0.30;
C.ps = 1;              % 1 = plane strain, 0 = plane stress

% Fracture / initiation parameters
C.sig_c = 300.0;
C.a0    = 0.004;

%% ================== LOADING ==================
C.load.type = 'remote_tension_y';
C.load.sig0 = 1.0;

%% ================== BOUNDARY CONDITIONS ==================
C.bc.anchor_mode = 'minimal';

%% ================== MESHING: HOLE-ONLY STAGE ==================
% Base mesh scale tied to the polygonal resolution of the circular hole
C.mesh1.hmin      = 2*pi*C.hole.r / C.hole.npoly;
C.mesh1.hmax      = 15*C.mesh1.hmin;
C.mesh1.hhole     = C.mesh1.hmin;
C.mesh1.hgrad     = 1.3;
C.mesh1.refineBox = [];

%% ================== MESHING: SHORT-CRACK / APPENDED-HOLE STAGE ==================
% Background Stage-II mesh uses the same scale as Stage I.
% Local sharp-tip refinement is imposed later in mesh_hole_pencil_domain
% through Hvertex/Hedge.
C.mesh2.hmax        = C.mesh1.hmax;
C.mesh2.hhole       = C.mesh1.hmin;
C.mesh2.hcrack      = C.mesh1.hmin;
C.mesh2.hgrad       = C.mesh1.hgrad;

% Current tested mouth half-shift for the appended hole
C.mesh2.chw         = 0.0008;

% Tip-circle radius for later J-integral / SIF work
C.mesh2.tip_radius  = 0.0020;
C.mesh2.tip_ncircle = 80;

%% ================== STAGE I: BOUNDARY STRESS SAMPLING ==================
C.stage1.nphi     = 1440;
C.stage1.avg_mode = 'none';
C.stage1.avg_rho  = 0.0;

%% ================== STAGE II: DIRECTION SEARCH ==================
C.stage2.angle_mode      = 'relative_to_inward_normal';
C.stage2.thmin_deg       = -75;
C.stage2.thmax_deg       =  75;
C.stage2.nth_coarse      = 61;
C.stage2.nth_fine        = 41;
C.stage2.fine_halfwin_deg = 5.0;
C.stage2.criterion       = 'local_symmetry';

%% ================== SIF / POSTPROCESSING ==================
C.sif.method     = 'circle2';
C.sif.decomp     = 'interaction_integral';
C.sif.tip_radius = C.mesh2.tip_radius;

%% ================== SOLVER ==================
C.solver.verbose        = 1;
C.solver.check_symmetry = false;
C.solver.linear_solver  = 'backslash';

%% ================== PLOTTING ==================
% Centralized here so the driver no longer overrides them
C.plot.show_mesh1        = true;
C.plot.show_mesh2        = true;
C.plot.show_stage1_stress = true;
C.plot.show_stage2_scan   = true;

end