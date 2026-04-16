function C = cfg_hole_initiation()
%CFG_HOLE_INITIATION  Configuration for the two-stage hole-initiation study.
%
% Two-stage model:
%   Stage I  : elastic plate with hole, determine initiation point on hole boundary
%   Stage II : insert a short crack at that point, determine initial direction
%
% Output:
%   C   struct with geometry, material, loading, meshing, and solver settings

%% ================== GEOMETRY ==================
% Rectangular plate:
%   x in [0, A], y in [-B, B]
C.A = 0.30;                  % plate length
C.B = 0.10;                  % half-height

% Circular hole
C.hole.type   = 'circle';
C.hole.center = [0.150, -0.060];
C.hole.r      = 0.030;
C.hole.npoly  = 240;         % polygonal resolution for geometry construction

% Optional: future extension to multiple holes
C.holes = {C.hole};

%% ================== MATERIAL ==================
C.E  = 210e3;                % Young's modulus
C.nu = 0.30;                 % Poisson ratio
C.ps = 1;                    % 1 = plane strain, 0 = plane stress

% Fracture / initiation parameters
C.sig_c = 300.0;             % critical tensile stress for initiation
C.a0    = 0.004;             % short inserted crack length

%% ================== LOADING ==================
% Use proportional loading:
%   t = lambda * tbar
%
% For now: remote uniform vertical traction on top/bottom edges
C.load.type = 'remote_tension_y';

% Unit reference load. Actual initiation load is lambda_ini * sig0
C.load.sig0 = 1.0;

%% ================== BOUNDARY CONDITIONS ==================
% Recommended default:
%   - left-bottom corner: ux = 0, uy = 0
%   - right-bottom corner: uy = 0
% This removes rigid motions without polluting the stress field too much.
C.bc.anchor_mode = 'minimal';

%% ================== MESHING: HOLE-ONLY STAGE ==================
C.mesh1.hmin      = 0.001; 
C.mesh1.hmax      = 0.01;   % global max element size
C.mesh1.hhole     = 0.001;  % target size near hole
C.mesh1.hgrad     = 1.3;    % mesh growth factor
C.mesh1.refineBox = [];      % optional [xmin xmax ymin ymax]

%% ================== MESHING: SHORT-CRACK STAGE ==================
C.mesh2.hmax         = 0.010;
C.mesh2.hhole        = 0.0025;
C.mesh2.hcrack       = 0.0012;   % target size near short crack/tip
C.mesh2.hgrad        = 1.25;
C.mesh2.chw          = 0.0015;   % half-width of pencil channel around short crack
C.mesh2.tip_radius   = 0.0020;   % radius of J-integral / SIF extraction contour
C.mesh2.tip_ncircle  = 80;       % points for circle sampling if needed

%% ================== STAGE I: BOUNDARY STRESS SAMPLING ==================
C.stage1.nphi       = 1440;   % number of angular sampling points on hole boundary
C.stage1.avg_mode   = 'none'; % 'none' or 'arc_average'
C.stage1.avg_rho    = 0.0;    % arc half-length for averaging, if used

%% ================== STAGE II: DIRECTION SEARCH ==================
% Angle theta is measured from the outward hole normal at the initiation point.
% Physically admissible crack directions usually go into the body,
% so default search is centered around the inward normal = pi from outward normal.
C.stage2.angle_mode   = 'relative_to_inward_normal';
C.stage2.thmin_deg    = -75;
C.stage2.thmax_deg    =  75;
C.stage2.nth_coarse   = 61;
C.stage2.nth_fine     = 41;
C.stage2.fine_halfwin_deg = 5.0;

% Direction criterion
C.stage2.criterion = 'local_symmetry';   % 'local_symmetry' or 'MTS'

%% ================== SIF / POSTPROCESSING ==================
C.sif.method = 'circle2';    % use SIF_LEFM_circle2
C.sif.decomp = 'interaction_integral';  % documentation only
C.sif.tip_radius = C.mesh2.tip_radius;

%% ================== SOLVER ==================
C.solver.verbose = 1;
C.solver.check_symmetry = false;
C.solver.linear_solver = 'backslash';

%% ================== PLOTTING ==================
C.plot.show_mesh1 = false;
C.plot.show_mesh2 = true;
C.plot.show_stage1_stress = true;
C.plot.show_stage2_scan = true;

end