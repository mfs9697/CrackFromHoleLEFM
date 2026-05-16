function Results = main_stage1_centered_hole_stress_check()
%MAIN_STAGE1_CENTERED_HOLE_STRESS_CHECK
% Stage-I diagnostic driver for a centered circular hole.
%
% Purpose:
%   1) Put the circular hole at the plate center.
%   2) Check the boundary tangential-stress distribution.
%   3) Check the stress concentration factor.
%   4) Verify whether the numerical maximum occurs at the expected
%      symmetry point.
%   5) Optionally force the initiation point for Stage-II validation.

    clc;
    close all;

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(repoRoot));

    fprintf('\n============================================================\n');
    fprintf('MAIN_STAGE1_CENTERED_HOLE_STRESS_CHECK\n');
    fprintf('Centered-hole Stage-I stress-concentration diagnostic\n');
    fprintf('============================================================\n');

    %% Configuration
    C0 = cfg_hole_initiation();

    % Centered-hole benchmark
    C0.hole.center = [0.5*C0.A, 0.0];
    C0.holes = {C0.hole};

    % Symmetry-friendly polygon and sampling
    C0.hole.npoly = 360;
    C0.mesh1.hmin  = 2*pi*C0.hole.r / C0.hole.npoly;
    C0.mesh1.hhole = C0.mesh1.hmin;
    C0.mesh1.hmax  = 10*C0.mesh1.hmin;
    C0.stage1.nphi = 1440;

    % Mesh-density sensitivity
    meshScaleList = [1.0 0.75 0.5 0.25];

    % Expected symmetry point for vertical tension in our convention.
    % Adjust if we later decide the expected point is pi instead.
    phiExpected = 0.0;

    rows = [];

    CaseStore = struct();
    CaseStore.items = {};
    CaseStore.count = 0;

    for im = 1:numel(meshScaleList)

        meshScale = meshScaleList(im);
        C = local_scale_stage1_mesh(C0, meshScale);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Stage-I mesh scale = %.3f\n', meshScale);
        fprintf('------------------------------------------------------------\n');

        G  = geom_hole_only(C);
        S1 = solve_hole_only(C, G, 'lambda', 1.0);
        B  = sample_hole_boundary_stress(C, G, S1);
        Iauto = find_hole_initiation_point(C, B);
        
        C.mesh1.refineHoleAngles = Iauto.phi_star;
        C.mesh1.refineAngleHalfWidth = deg2rad(8);
        C.mesh1.hrefine = 0.5*C.mesh1.hhole;

        Iauto = find_hole_initiation_point(C, B);

        Iforced = stage1_force_initiation_phi(C, B, phiExpected);

        sig0 = C.load.sig0;

        Kt_auto   = Iauto.sig_tt_pos_unit / sig0;
        Kt_forced = Iforced.sig_tt_pos_unit / sig0;

        dphi_auto = local_angle_diff(Iauto.phi_star, phiExpected);

        rows = [rows; ...
            meshScale, ...
            size(G.p,1), size(G.t,1), ...
            C.hole.npoly, C.stage1.nphi, ...
            Iauto.phi_star, rad2deg(Iauto.phi_star), rad2deg(dphi_auto), ...
            Iauto.sig_tt_pos_unit, Kt_auto, ...
            Iforced.phi_star, rad2deg(Iforced.phi_star), ...
            Iforced.sig_tt_pos_unit, Kt_forced, ...
            B.eps_shift]; %#ok<AGROW>

        fprintf('  auto phi     = %.8f rad = %.5f deg\n', ...
            Iauto.phi_star, rad2deg(Iauto.phi_star));
        fprintf('  forced phi   = %.8f rad = %.5f deg\n', ...
            Iforced.phi_star, rad2deg(Iforced.phi_star));
        fprintf('  Kt auto      = %.8f\n', Kt_auto);
        fprintf('  Kt forced    = %.8f\n', Kt_forced);

        CaseStore.count = CaseStore.count + 1;
        CaseStore.items{CaseStore.count} = struct( ...
            'meshScale', meshScale, ...
            'C', C, 'G', G, 'S1', S1, 'B', B, ...
            'Iauto', Iauto, 'Iforced', Iforced);
    end

    T = array2table(rows, ...
        'VariableNames', { ...
            'meshScale', 'nNodeT3', 'nElemT3', ...
            'npoly', 'nphi', ...
            'phi_auto_rad', 'phi_auto_deg', 'phi_auto_minus_expected_deg', ...
            'sig_tt_auto', 'Kt_auto', ...
            'phi_forced_rad', 'phi_forced_deg', ...
            'sig_tt_forced', 'Kt_forced', ...
            'eps_shift'});

    fprintf('\n=== Stage-I centered-hole stress-concentration table ===\n');
    disp(T);

    Results = struct();
    Results.Table = T;
    Results.CaseStore = CaseStore;
    Results.phiExpected = phiExpected;

    % Plot the finest case
    if CaseStore.count > 0
        last = CaseStore.items{CaseStore.count};
        stage1_plot_boundary_stress(last.B, last.Iauto, last.Iforced);
    end
end