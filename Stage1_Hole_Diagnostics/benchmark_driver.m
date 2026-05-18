addpath(genpath(pwd))

C = cfg_hole_initiation();

C.hole.npoly = 360;
C.mesh1.hmin  = 2*pi*C.hole.r / C.hole.npoly;
C.mesh1.hhole = C.mesh1.hmin;
C.mesh1.hmax  = 10*C.mesh1.hmin;

P = stage1_pde_toolbox_circular_hole_benchmark(C);

fprintf('PDE Toolbox benchmark Kt = %.8f\n', P.Kt);
fprintf('Peak phi = %.5f deg\n', P.phi_peak_deg);