function Results = main_crack_path_LEFM_from_hole()
%MAIN_CRACK_PATH_LEFM_FROM_HOLE
% First LEFM crack-path driver for a crack initiated from a circular hole.
%
% Purpose:
%   1) Solve the hole-only problem and determine the initiation point.
%   2) Build a short appended crack normal to the hole boundary.
%   3) Solve the cracked problem and compute KI, KII.
%   4) Vary the appendix angle theta and find theta such that KII(theta)=0.
%
% Convention:
%   theta = 0 means crack extension along the material-side normal to the
%   hole boundary.
%
% Current scope:
%   This is the first one-step LEFM study. The MTS routine
%   kink_angle_LEFM_MTS.m is not yet used as the main propagation rule here.
%   It will be used later in the incremental crack-growth study, after the
%   first crack segment has been established.

    clc;
    close all;

    addpath(genpath(pwd));

    fprintf('\n============================================================\n');
    fprintf('MAIN_CRACK_PATH_LEFM_FROM_HOLE\n');
    fprintf('One-step LEFM crack initiation/extension from a hole\n');
    fprintf('============================================================\n');

    %% ============================================================
    % 0. Configuration
    %% ============================================================
    C = cfg_hole_initiation();

    % Trial appendix-angle sweep.
    % theta is measured from the inward material normal.
    thetaDegList = -12:1:8;
    thetaList    = deg2rad(thetaDegList);

    % Radius for circular J-integral / SIF contour.
    % If empty, SIF_LEFM_circle2 chooses its internal default.
    rI_user = [];

    % Plot controls
    doPlotStage1        = true;
    doPlotNormalMesh    = true;
    doPlotSweepCurves   = true;
    doPlotRootMesh      = true;

    %% ============================================================
    % 1. Stage I: hole-only problem
    %% ============================================================
    fprintf('\n=== Stage I: hole-only solution and initiation point ===\n');

    G  = geom_hole_only(C);
    S1 = solve_hole_only(C, G, 'lambda', 1.0);
    B  = sample_hole_boundary_stress(C, G, S1);
    I  = find_hole_initiation_point(C, B);

    fprintf('\nInitiation point:\n');
    fprintf('  phi_*          = %.8f rad = %.4f deg\n', ...
        I.phi_star, rad2deg(I.phi_star));
    fprintf('  x_*            = [%.10e, %.10e]\n', ...
        I.x_star(1), I.x_star(2));
    fprintf('  n_mat          = [%.10e, %.10e]\n', ...
        I.n_mat_star(1), I.n_mat_star(2));
    fprintf('  t_hat          = [%.10e, %.10e]\n', ...
        I.t_hat_star(1), I.t_hat_star(2));
    fprintf('  sigma_tt,max   = %.10e at unit load\n', ...
        I.sig_tt_pos_unit);
    fprintf('  lambda_ini     = %.10e\n', I.lambda_ini);
    fprintf('  applied_ini    = %.10e\n', I.sig_applied_ini);

    if doPlotStage1
        local_plot_stage1_boundary_stress(B, I);
    end

    %% ============================================================
    % 2. Normal appendix test: theta = 0
    %% ============================================================
    fprintf('\n=== Stage II-A: normal appendix test, theta = 0 deg ===\n');

    theta0 = 0.0;

    try
        [G2_0, D0, M0, Mc0] = build_stage2_cracked_mesh_for_theta( ...
            C, I, theta0, ...
            'PlotGeom', doPlotNormalMesh, ...
            'PlotMesh', doPlotNormalMesh, ...
            'PlotCollapsed', doPlotNormalMesh);

        S2_0 = solve_cracked_LEFM(C, Mc0);

        [KI0, KII0, Aux0] = compute_SIF_for_stage2( ...
            C, G2_0, Mc0, S2_0, ...
            'rI', rI_user);

        fprintf('\nNormal appendix SIFs:\n');
        fprintf('  theta          = %.6f deg\n', rad2deg(theta0));
        fprintf('  KI             = %.10e\n', KI0);
        fprintf('  KII            = %.10e\n', KII0);
        fprintf('  KII/KI         = %.10e\n', KII0 / KI0);

    catch ME
        fprintf('\nERROR during normal appendix SIF test:\n');
        fprintf('  %s\n', ME.message);
        rethrow(ME);
    end

    %% ============================================================
    % 3. Appendix-angle sweep
    %% ============================================================
    fprintf('\n=== Stage II-B: appendix-angle sweep ===\n');

    nTheta = numel(thetaList);

    Sweep = struct();
    Sweep.theta      = thetaList(:);
    Sweep.thetaDeg   = thetaDegList(:);
    Sweep.KI         = nan(nTheta,1);
    Sweep.KII        = nan(nTheta,1);
    Sweep.KII_over_KI = nan(nTheta,1);
    Sweep.valid      = false(nTheta,1);
    Sweep.message    = strings(nTheta,1);

    for it = 1:nTheta
        theta = thetaList(it);

        fprintf('\n--- theta %3d/%3d: %+8.3f deg ---\n', ...
            it, nTheta, rad2deg(theta));

        try
            [KI, KII, Aux] = local_eval_SIF_for_theta(C, I, theta, rI_user);

            Sweep.KI(it)          = KI;
            Sweep.KII(it)         = KII;
            Sweep.KII_over_KI(it) = KII / KI;
            Sweep.valid(it)       = true;
            Sweep.message(it)     = "ok";

            fprintf('  KI      = %.10e\n', KI);
            fprintf('  KII     = %.10e\n', KII);
            fprintf('  KII/KI  = %.10e\n', KII/KI);

        catch ME
            Sweep.valid(it)   = false;
            Sweep.message(it) = string(ME.message);

            fprintf('  FAILED: %s\n', ME.message);
        end
    end

    if doPlotSweepCurves
        local_plot_sweep(Sweep);
    end

    %% ============================================================
    % 4. Find theta such that KII(theta) = 0
    %% ============================================================
    fprintf('\n=== Stage II-C: root search for KII(theta) = 0 ===\n');

    idxValid = find(Sweep.valid & isfinite(Sweep.KII));

    if numel(idxValid) < 2
        warning('Not enough valid sweep points to search for KII=0.');
        Root = struct();
        Root.success = false;
        Root.message = 'Not enough valid sweep points.';
    else
        [bracket, hasBracket] = local_find_sign_change_bracket( ...
            Sweep.theta(idxValid), Sweep.KII(idxValid));

        if hasBracket
            fprintf('Sign-change bracket found:\n');
            fprintf('  theta_a = %.6f deg\n', rad2deg(bracket(1)));
            fprintf('  theta_b = %.6f deg\n', rad2deg(bracket(2)));

            fKII = @(th) local_eval_KII_for_theta(C, I, th, rI_user);

            try
                thetaStar = fzero(fKII, bracket);

                [G2s, Ds, Ms, Mcs] = build_stage2_cracked_mesh_for_theta( ...
                    C, I, thetaStar, ...
                    'PlotGeom', doPlotRootMesh, ...
                    'PlotMesh', doPlotRootMesh, ...
                    'PlotCollapsed', doPlotRootMesh);

                S2s = solve_cracked_LEFM(C, Mcs);

                [KIstar, KIIstar, AuxStar] = compute_SIF_for_stage2( ...
                    C, G2s, Mcs, S2s, ...
                    'rI', rI_user);

                Root = struct();
                Root.success      = true;
                Root.theta        = thetaStar;
                Root.thetaDeg     = rad2deg(thetaStar);
                Root.KI           = KIstar;
                Root.KII          = KIIstar;
                Root.KII_over_KI  = KIIstar / KIstar;
                Root.G2           = G2s;
                Root.D            = Ds;
                Root.M            = Ms;
                Root.Mc           = Mcs;
                Root.S2           = S2s;
                Root.AuxSIF       = AuxStar;
                Root.message      = 'ok';

                fprintf('\nRoot-refined appendix angle:\n');
                fprintf('  theta_*        = %.10f rad = %.6f deg\n', ...
                    Root.theta, Root.thetaDeg);
                fprintf('  KI(theta_*)    = %.10e\n', Root.KI);
                fprintf('  KII(theta_*)   = %.10e\n', Root.KII);
                fprintf('  KII/KI         = %.10e\n', Root.KII_over_KI);

            catch ME
                Root = struct();
                Root.success = false;
                Root.message = ME.message;

                fprintf('Root search failed:\n');
                fprintf('  %s\n', ME.message);
            end

        else
            fprintf('No sign-change bracket for KII was found.\n');
            fprintf('Using minimum |KII/KI| from the valid sweep as fallback.\n');

            ratio = abs(Sweep.KII_over_KI(idxValid));
            [~, jmin] = min(ratio);
            ibest = idxValid(jmin);

            Root = struct();
            Root.success      = false;
            Root.theta        = Sweep.theta(ibest);
            Root.thetaDeg     = Sweep.thetaDeg(ibest);
            Root.KI           = Sweep.KI(ibest);
            Root.KII          = Sweep.KII(ibest);
            Root.KII_over_KI  = Sweep.KII_over_KI(ibest);
            Root.message      = 'No sign-change bracket; fallback to min |KII/KI|.';

            fprintf('\nFallback angle:\n');
            fprintf('  theta_best     = %.6f deg\n', Root.thetaDeg);
            fprintf('  KI             = %.10e\n', Root.KI);
            fprintf('  KII            = %.10e\n', Root.KII);
            fprintf('  KII/KI         = %.10e\n', Root.KII_over_KI);
        end
    end

    %% ============================================================
    % 5. Optional: compute MTS angle for diagnostic only
    %% ============================================================
    % This is not the main rule in the current one-step appendix-angle study.
    % It is included only as a diagnostic and for later crack-growth work.

    if exist('KI0', 'var') && exist('KII0', 'var')
        [thetaMTS0, thetaMTS0Deg, sigttMTS0] = kink_angle_LEFM_MTS(KI0, KII0);

        fprintf('\nMTS diagnostic for normal appendix SIFs:\n');
        fprintf('  theta_MTS      = %.10f rad = %.6f deg\n', ...
            thetaMTS0, thetaMTS0Deg);
        fprintf('  sigma_tt_MTS   = %.10e\n', sigttMTS0);
    end

    %% ============================================================
    % 6. Store output
    %% ============================================================
    Results = struct();

    Results.C = C;

    Results.Stage1 = struct();
    Results.Stage1.G  = G;
    Results.Stage1.S1 = S1;
    Results.Stage1.B  = B;
    Results.Stage1.I  = I;

    Results.Normal = struct();
    Results.Normal.theta = theta0;
    Results.Normal.thetaDeg = rad2deg(theta0);
    Results.Normal.G2 = G2_0;
    Results.Normal.D  = D0;
    Results.Normal.M  = M0;
    Results.Normal.Mc = Mc0;
    Results.Normal.S2 = S2_0;
    Results.Normal.KI = KI0;
    Results.Normal.KII = KII0;
    Results.Normal.KII_over_KI = KII0 / KI0;
    Results.Normal.AuxSIF = Aux0;

    Results.Sweep = Sweep;
    Results.Root  = Root;

    fprintf('\n============================================================\n');
    fprintf('MAIN_CRACK_PATH_LEFM_FROM_HOLE finished.\n');
    fprintf('============================================================\n');
end


% ========================================================================
% Local helper: evaluate KI, KII for one appendix angle
% ========================================================================

function [KI, KII, Aux] = local_eval_SIF_for_theta(C, I, theta, rI_user)

    [G2, D, M, Mc] = build_stage2_cracked_mesh_for_theta( ...
        C, I, theta, ...
        'PlotGeom', false, ...
        'PlotMesh', false, ...
        'PlotCollapsed', false);

    S2 = solve_cracked_LEFM(C, Mc);

    [KI, KII, AuxSIF] = compute_SIF_for_stage2( ...
        C, G2, Mc, S2, ...
        'rI', rI_user);

    Aux = struct();
    Aux.G2 = G2;
    Aux.D  = D;
    Aux.M  = M;
    Aux.Mc = Mc;
    Aux.S2 = S2;
    Aux.SIF = AuxSIF;
end


% ========================================================================
% Local helper: evaluate only KII for root search
% ========================================================================

function KII = local_eval_KII_for_theta(C, I, theta, rI_user)

    [~, KII, ~] = local_eval_SIF_for_theta(C, I, theta, rI_user);
end


% ========================================================================
% Local helper: find sign-change bracket
% ========================================================================

function [bracket, hasBracket] = local_find_sign_change_bracket(theta, KII)

    bracket = [nan, nan];
    hasBracket = false;

    theta = theta(:);
    KII   = KII(:);

    good = isfinite(theta) & isfinite(KII);
    theta = theta(good);
    KII   = KII(good);

    if numel(theta) < 2
        return;
    end

    % Prefer a sign change closest to theta = 0.
    candidates = [];

    for i = 1:numel(theta)-1
        if KII(i) == 0
            bracket = [theta(i), theta(i)];
            hasBracket = true;
            return;
        end

        if KII(i) * KII(i+1) < 0
            midAbs = abs(0.5*(theta(i) + theta(i+1)));
            candidates = [candidates; i, midAbs]; %#ok<AGROW>
        end
    end

    if isempty(candidates)
        return;
    end

    [~, j] = min(candidates(:,2));
    iBest = candidates(j,1);

    bracket = [theta(iBest), theta(iBest+1)];
    hasBracket = true;
end


% ========================================================================
% Local plotting: Stage I boundary stress
% ========================================================================

function local_plot_stage1_boundary_stress(B, I)

    figure('Name', 'Stage I: hole-boundary stress', 'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    plot(rad2deg(B.phi), B.sig_tt_eff, 'LineWidth', 1.2);
    plot(rad2deg(I.phi_star), I.sig_tt_unit, ...
        'ro', 'MarkerSize', 7, 'LineWidth', 1.5);

    xlabel('\phi, deg');
    ylabel('\sigma_{\theta\theta}^{eff} at unit load');
    title('Stage I: tangential stress along hole boundary');
    xlim([0, 360]);
end


% ========================================================================
% Local plotting: sweep curves
% ========================================================================

function local_plot_sweep(Sweep)

    valid = Sweep.valid & isfinite(Sweep.KI) & isfinite(Sweep.KII);

    if ~any(valid)
        warning('No valid sweep points to plot.');
        return;
    end

    th = Sweep.thetaDeg(valid);

    figure('Name', 'LEFM appendix-angle sweep: SIFs', 'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    plot(th, Sweep.KI(valid),  'o-', 'LineWidth', 1.2);
    plot(th, Sweep.KII(valid), 's-', 'LineWidth', 1.2);
    yline(0, 'k--');

    xlabel('\theta, deg');
    ylabel('SIF at unit nominal load');
    legend('K_I', 'K_{II}', 'Location', 'best');
    title('SIFs versus appendix angle');

    figure('Name', 'LEFM appendix-angle sweep: KII/KI', 'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    plot(th, Sweep.KII_over_KI(valid), 'o-', 'LineWidth', 1.2);
    yline(0, 'k--');

    xlabel('\theta, deg');
    ylabel('K_{II}/K_I');
    title('Local-symmetry condition for appendix angle');
end