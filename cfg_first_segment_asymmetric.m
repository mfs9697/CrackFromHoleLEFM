function C=cfg_first_segment_asymmetric()
%CFG_FIRST_SEGMENT_ASYMMETRIC Frozen starting configuration for crack-path work.
%
% This configuration is the accepted asymmetric full-domain Stage-I benchmark
% used in audit Steps 19--38, promoted here as the starting point for the
% first-segment local-symmetry study.
%
% It is intentionally explicit about the parameters that define the physical
% starting state.  The generic cfg_hole_initiation.m is currently also used
% by centered-hole verification cases and must not silently redefine this
% asymmetric crack-path problem.

C=cfg_asymmetric_full_domain();

% -------------------- physical geometry --------------------
C.A=0.30;
C.B=0.10;
C.hole.type='circle';
C.hole.center=[0.17,-0.02];
C.hole.r=0.030;
C.hole.npoly=480;
C.holes={C.hole};

% -------------------- material and initiation --------------------
C.E=210e3;
C.nu=0.30;
C.ps=1;
C.sig_c=300.0;
C.a0=0.004; % 4-mm first-segment regularization for the next research stage

% -------------------- loading / anchoring --------------------
C.load.type='remote_tension_y';
C.load.sig0=1.0;
C.bc.anchor_mode='minimal';

% -------------------- Stage-I mesh --------------------
hArc=2*pi*C.hole.r/C.hole.npoly;
C.mesh1.hmin=hArc;
C.mesh1.hhole=hArc;
C.mesh1.hmax=20*hArc;
C.mesh1.hgrad=1.2;
C.mesh1.refineBox=[];

% Stage-II background values are frozen now only as metadata.  They are NOT
% yet accepted as the final first-segment mesh; the EDI-compatible paired
% tip/core construction will be qualified separately before the angle sweep.
C.mesh2.hmax=C.mesh1.hmax;
C.mesh2.hhole=C.mesh1.hhole;
C.mesh2.hcrack=C.mesh1.hmin;
C.mesh2.hgrad=C.mesh1.hgrad;
C.mesh2.chw=0.0001;

% -------------------- Stage-I estimator --------------------
C.stage1.method='boundary_extrapolated_t6';
C.stage1.nphi=1440;
C.stage1.periodic=true;
C.stage1.shift_fractions=[0.05 0.10 0.25];
C.stage1.radial_fit_order=1;
C.stage1.angular_fit_enable=true;
C.stage1.angular_fit_points=5;
C.stage1.angular_fit_halfwidth_factor=3.0;
C.stage1.max_tie_rel_tol=1e-8;

% -------------------- future Stage-II direction search --------------------
C.stage2.angle_mode='relative_to_inward_normal';
C.stage2.thmin_deg=-75;
C.stage2.thmax_deg=75;
C.stage2.nth_coarse=61;
C.stage2.nth_fine=41;
C.stage2.fine_halfwin_deg=5.0;
C.stage2.criterion='local_symmetry';

% The old circular-J selector is deliberately NOT promoted here as the
% first-segment observable.  The next stage will define a separately
% qualified interaction-EDI configuration for the 4-mm trial crack.
C.sif.method='pending_qualified_interaction_EDI';

% Quiet Stage-I solve for reproducible baseline generation.
C.solver.verbose=0;
C.plot.show_mesh1=false;
C.plot.show_mesh2=false;
C.plot.show_stage1_stress=false;
C.plot.show_stage2_scan=false;
end
