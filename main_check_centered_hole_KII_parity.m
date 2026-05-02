function Results = main_check_centered_hole_KII_parity()
%MAIN_CHECK_CENTERED_HOLE_KII_PARITY
% Parity diagnostic for KII in the centered-hole benchmark.
%
% Purpose:
%   1) Put the circular hole at the center of the plate.
%   2) Use a symmetric benchmark under remote vertical tension.
%   3) Build short cracks at angles +/-theta relative to the inward normal.
%   4) Compute KI,KII using different SIF debug variants.
%   5) Check the parity:
%
%          KII_odd(theta)  = 0.5*(KII(+theta)-KII(-theta))
%          KII_even(theta) = 0.5*(KII(+theta)+KII(-theta))
%
%      For a clean sign convention in a symmetric benchmark, the even part
%      should be small compared with the odd part.
%
% Notes:
%   - This driver does not choose a crack-growth angle.
%   - It is only a diagnostic test for the SIF extraction convention.
%   - It assumes that the following functions are on the MATLAB path:
%       cfg_hole_initiation
%       geom_hole_only
%       solve_hole_only
%       sample_hole_boundary_stress
%       find_hole_initiation_point
%       build_stage2_cracked_mesh_for_theta
%       solve_cracked_LEFM
%       SIF_LEFM_circle2_debug
%       SIF_LEFM_circle2_debug_swapPQ

    clc;
    close all;
    addpath(genpath(pwd));

    fprintf('\n============================================================\n');
    fprintf('MAIN_CHECK_CENTERED_HOLE_KII_PARITY\n');
    fprintf('Centered-hole parity check for KII\n');
    fprintf('============================================================\n');

    %% ============================================================
    % 0. Configuration
    %% ============================================================
    C0 = cfg_hole_initiation();

    % Center the hole.
    C0.hole.center = [0.5*C0.A, 0.0];
    C0.holes = {C0.hole};

    % Diagnostic settings, based on the stable centered-hole tests.
    a0_list        = 4*C0.a0;
    meshScaleList  = 0.4;
    thetaDegProbe  = [-10 -5 0 5 10];
    thetaProbe     = deg2rad(thetaDegProbe);

    % Relative contour radii.
    rI_over_a0_list = [0.5 0.6 0.7];

    % SIF variants to compare.
    variantList = {'original', 'swapPQ'};

    % Debug controls.
    nthet        = 200;
    debugVerbose = false;
    debugPlot    = false;

    % Plot controls.
    doPlotStage1 = true;
    doPlotMeshes = false;

    %% ============================================================
    % 1. Stage I: centered-hole solution and initiation point
    %% ============================================================
    fprintf('\n=== Stage I: centered-hole solution and initiation point ===\n');

    % Use baseline centered config for Stage I. The selected point is then
    % reused in Stage II. This follows the current two-stage workflow.
    G  = geom_hole_only(C0);
    S1 = solve_hole_only(C0, G, 'lambda', 1.0);
    B  = sample_hole_boundary_stress(C0, G, S1);
    I  = find_hole_initiation_point(C0, B);

    fprintf('\nCentered-hole initiation point selected by current rule:\n');
    fprintf('  phi_*          = %.8f rad = %.4f deg\n', ...
        I.phi_star, rad2deg(I.phi_star));
    fprintf('  x_*            = [%.10e, %.10e]\n', ...
        I.x_star(1), I.x_star(2));
    fprintf('  n_mat          = [%.10e, %.10e]\n', ...
        I.n_mat_star(1), I.n_mat_star(2));
    fprintf('  t_hat          = [%.10e, %.10e]\n', ...
        I.t_hat_star(1), I.t_hat_star(2));
    fprintf('  lambda_ini     = %.10e\n', I.lambda_ini);

    if doPlotStage1
        local_plot_stage1_boundary_stress(B, I);
    end

    %% ============================================================
    % 2. SIF parity loop
    %% ============================================================
    fprintf('\n=== Stage II: centered-hole KII parity loop ===\n');

    rows = [];

    for ia = 1:numel(a0_list)

        a0 = a0_list(ia);
        a0_over_base = a0 / C0.a0;

        for im = 1:numel(meshScaleList)

            meshScale = meshScaleList(im);

            C = local_scale_mesh_config(C0, meshScale);
            C.a0 = a0;

            fprintf('\n============================================================\n');
            fprintf('a0 = %.8e | a0/a0_base = %.3f | meshScale = %.3f\n', ...
                a0, a0_over_base, meshScale);
            fprintf('============================================================\n');

            for it = 1:numel(thetaProbe)

                theta = thetaProbe(it);
                thetaDeg = thetaDegProbe(it);

                fprintf('\n--- theta = %+8.3f deg ---\n', thetaDeg);

                try
                    [G2, D, M, Mc] = build_stage2_cracked_mesh_for_theta( ...
                        C, I, theta, ...
                        'PlotGeom', doPlotMeshes, ...
                        'PlotMesh', doPlotMeshes, ...
                        'PlotCollapsed', doPlotMeshes);

                    S2 = solve_cracked_LEFM(C, Mc);

                    nNode3 = size(S2.mesh.coord3, 1);
                    nElem3 = size(S2.mesh.connect3, 1);
                    nNode6 = size(S2.mesh.coord, 1);
                    nElem6 = size(S2.mesh.connect, 1);

                    fprintf('  mesh: T3 nodes=%d, T3 elems=%d | T6 nodes=%d, T6 elems=%d\n', ...
                        nNode3, nElem3, nNode6, nElem6);

                    for ir = 1:numel(rI_over_a0_list)

                        rfac = rI_over_a0_list(ir);
                        rI = rfac * G2.crack.a0;

                        for iv = 1:numel(variantList)

                            variant = variantList{iv};

                            try
                                [KI, KII, AuxSIF] = local_compute_SIF_variant( ...
                                    variant, C, G2, Mc, S2, rI, ...
                                    'nthet', nthet, ...
                                    'verbose', debugVerbose, ...
                                    'plot', debugPlot);

                                ratio = KII / KI;

                                [JII_over_absint, J2_sign_changes, ...
                                    nearP_1e4, nearQ_1e4, ...
                                    minBaryP, minBaryQ] = local_extract_debug_metrics(AuxSIF);

                                rows = [rows; ...
                                    a0, a0_over_base, meshScale, ...
                                    thetaDeg, rfac, rI, ...
                                    iv, ...
                                    KI, KII, ratio, ...
                                    nNode3, nElem3, nNode6, nElem6, ...
                                    JII_over_absint, J2_sign_changes, ...
                                    nearP_1e4, nearQ_1e4, minBaryP, minBaryQ, ...
                                    true]; %#ok<AGROW>

                                fprintf(['  %-8s | rI/a0 = %.3f | KI = %.8e | ', ...
                                         'KII = %+ .8e | KII/KI = %+ .8e'], ...
                                    variant, rfac, KI, KII, ratio);

                                if isfinite(JII_over_absint)
                                    fprintf(' | JII/int|J2| = %+ .3e | J2 changes = %g', ...
                                        JII_over_absint, J2_sign_changes);
                                end
                                fprintf('\n');

                            catch MEinner
                                rows = [rows; ...
                                    a0, a0_over_base, meshScale, ...
                                    thetaDeg, rfac, rI, ...
                                    iv, ...
                                    nan, nan, nan, ...
                                    nNode3, nElem3, nNode6, nElem6, ...
                                    nan, nan, nan, nan, nan, nan, ...
                                    false]; %#ok<AGROW>

                                fprintf('  %-8s | rI/a0 = %.3f FAILED: %s\n', ...
                                    variant, rfac, MEinner.message);
                            end
                        end
                    end

                catch MEouter
                    fprintf('  theta failed before SIF loop: %s\n', MEouter.message);

                    for ir = 1:numel(rI_over_a0_list)
                        rfac = rI_over_a0_list(ir);
                        rI = rfac * C.a0;

                        for iv = 1:numel(variantList)
                            rows = [rows; ...
                                a0, a0_over_base, meshScale, ...
                                thetaDeg, rfac, rI, ...
                                iv, ...
                                nan, nan, nan, ...
                                nan, nan, nan, nan, ...
                                nan, nan, nan, nan, nan, nan, ...
                                false]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end

    %% ============================================================
    % 3. Full table
    %% ============================================================
    T = array2table(rows, ...
        'VariableNames', { ...
            'a0', 'a0_over_base', 'meshScale', ...
            'thetaDeg', 'rI_over_a0', 'rI', ...
            'variantID', ...
            'KI', 'KII', 'KII_over_KI', ...
            'nNode3', 'nElem3', 'nNode6', 'nElem6', ...
            'JII_over_absint', 'J2_sign_changes', ...
            'nearP_1e4', 'nearQ_1e4', 'minBaryP', 'minBaryQ', ...
            'valid'});

    T.variant = strings(height(T),1);
    for iv = 1:numel(variantList)
        T.variant(T.variantID == iv) = string(variantList{iv});
    end
    T = movevars(T, 'variant', 'After', 'variantID');

    fprintf('\n=== Full centered-hole parity table ===\n');
    disp(T);

    %% ============================================================
    % 4. Summary over contour radii
    %% ============================================================
    Summary = local_make_variant_summary(T);

    fprintf('\n=== Summary over contour radii for each variant and theta ===\n');
    disp(Summary);

    %% ============================================================
    % 5. Parity table
    %% ============================================================
    Parity = local_make_parity_table(Summary);

    fprintf('\n=== KII parity table based on median KII/KI ===\n');
    disp(Parity);

    %% ============================================================
    % 6. Extra global-sign diagnostic for original variant
    %% ============================================================
    ParityFlip = Parity;
    if ~isempty(ParityFlip)
        idxOrig = strcmp(ParityFlip.variant, "original");
        ParityFlip.variant(idxOrig) = "original_flipSign";
        ParityFlip.Kminus(idxOrig) = -ParityFlip.Kminus(idxOrig);
        ParityFlip.Kplus(idxOrig)  = -ParityFlip.Kplus(idxOrig);
        ParityFlip.Keven(idxOrig)  = -ParityFlip.Keven(idxOrig);
        ParityFlip.Kodd(idxOrig)   = -ParityFlip.Kodd(idxOrig);
    end

    %% ============================================================
    % 7. Plots
    %% ============================================================
    local_plot_variant_radius(T);
    local_plot_variant_theta_summary(Summary);
    local_plot_parity(Parity);

    %% ============================================================
    % 8. Store output
    %% ============================================================
    Results = struct();
    Results.C0 = C0;
    Results.C = C;
    Results.Stage1.G = G;
    Results.Stage1.S1 = S1;
    Results.Stage1.B = B;
    Results.Stage1.I = I;

    Results.Table = T;
    Results.Summary = Summary;
    Results.Parity = Parity;
    Results.ParityFlip = ParityFlip;

    Results.variantList = variantList;
    Results.a0_list = a0_list;
    Results.meshScaleList = meshScaleList;
    Results.thetaDegProbe = thetaDegProbe;
    Results.rI_over_a0_list = rI_over_a0_list;

    fprintf('\n============================================================\n');
    fprintf('MAIN_CHECK_CENTERED_HOLE_KII_PARITY finished.\n');
    fprintf('============================================================\n');
end


% ========================================================================
% Compute SIF for selected variant
% ========================================================================

function [KI, KII, Aux] = local_compute_SIF_variant(variant, C, G2, Mc, S2, rI, varargin)

    ip = inputParser;
    addParameter(ip, 'verbose', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'plot', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'nthet', 200, @(x)isnumeric(x) && isscalar(x) && x >= 10);
    parse(ip, varargin{:});

    meshSIF = S2.mesh;

    matSIF = S2.mat;
    if ~isfield(matSIF, 'Dmat')
        if isfield(matSIF, 'D')
            matSIF.Dmat = matSIF.D;
        else
            error('local_compute_SIF_variant:MissingDmat', ...
                'Material struct must contain Dmat or D.');
        end
    end

    if isfield(Mc, 'crack') && isfield(Mc.crack, 'Pmid') && ~isempty(Mc.crack.Pmid)
        V = Mc.crack.Pmid;
    else
        V = G2.crack.polyline;
    end

    switch lower(strtrim(variant))
        case 'original'
            [KI, KII, Dbg] = SIF_LEFM_circle2_debug( ...
                meshSIF, S2.U, V, matSIF, rI, ...
                'nthet', ip.Results.nthet, ...
                'plot', logical(ip.Results.plot), ...
                'verbose', logical(ip.Results.verbose));

        case 'swappq'
            [KI, KII, Dbg] = SIF_LEFM_circle2_debug_swapPQ( ...
                meshSIF, S2.U, V, matSIF, rI, ...
                'nthet', ip.Results.nthet, ...
                'plot', logical(ip.Results.plot), ...
                'verbose', logical(ip.Results.verbose));

        otherwise
            error('local_compute_SIF_variant:UnknownVariant', ...
                'Unknown SIF variant "%s".', variant);
    end

    Aux = struct();
    Aux.meshSIF = meshSIF;
    Aux.matSIF = matSIF;
    Aux.V = V;
    Aux.rI = rI;
    Aux.Dbg = Dbg;
end


% ========================================================================
% Extract debug metrics
% ========================================================================

function [JII_over_absint, J2_sign_changes, nearP_1e4, nearQ_1e4, minBaryP, minBaryQ] = ...
    local_extract_debug_metrics(AuxSIF)

    JII_over_absint = nan;
    J2_sign_changes = nan;
    nearP_1e4 = nan;
    nearQ_1e4 = nan;
    minBaryP = nan;
    minBaryQ = nan;

    if ~isstruct(AuxSIF) || ~isfield(AuxSIF, 'Dbg')
        return;
    end

    Dbg = AuxSIF.Dbg;
    if ~isfield(Dbg, 'diagnostics')
        return;
    end

    DD = Dbg.diagnostics;

    JII_over_absint = getfield_if_exists(DD, 'JII_over_absint', nan);
    J2_sign_changes = getfield_if_exists(DD, 'J2_sign_changes', nan);
    nearP_1e4 = getfield_if_exists(DD, 'num_near_edge_P_1e4', nan);
    nearQ_1e4 = getfield_if_exists(DD, 'num_near_edge_Q_1e4', nan);
    minBaryP = getfield_if_exists(DD, 'min_baryP', nan);
    minBaryQ = getfield_if_exists(DD, 'min_baryQ', nan);
end


% ========================================================================
% Summary table for each variant and theta
% ========================================================================

function Summary = local_make_variant_summary(T)

    good = T.valid & isfinite(T.KI) & isfinite(T.KII_over_KI);

    aVals = unique(T.a0_over_base(good));
    mVals = unique(T.meshScale(good));
    vVals = unique(T.variant(good));
    thVals = unique(T.thetaDeg(good));

    rows = [];

    for ia = 1:numel(aVals)
        aa = aVals(ia);

        for im = 1:numel(mVals)
            ms = mVals(im);

            for iv = 1:numel(vVals)
                vv = vVals(iv);

                for it = 1:numel(thVals)
                    th = thVals(it);

                    idx = good ...
                        & abs(T.a0_over_base - aa) < 1e-12 ...
                        & abs(T.meshScale - ms) < 1e-12 ...
                        & T.variant == vv ...
                        & abs(T.thetaDeg - th) < 1e-12;

                    if ~any(idx)
                        continue;
                    end

                    valsKI = T.KI(idx);
                    valsRat = T.KII_over_KI(idx);

                    meanKI = mean(valsKI, 'omitnan');
                    stdKI  = std(valsKI, 'omitnan');

                    meanRat = mean(valsRat, 'omitnan');
                    stdRat  = std(valsRat, 'omitnan');
                    medRat  = median(valsRat, 'omitnan');
                    madRat  = median(abs(valsRat - medRat), 'omitnan');

                    minRat = min(valsRat);
                    maxRat = max(valsRat);
                    spreadRat = maxRat - minRat;

                    signStable = all(valsRat > 0) || all(valsRat < 0);

                    meanJIIabs = mean(T.JII_over_absint(idx), 'omitnan');
                    meanJ2chg  = mean(T.J2_sign_changes(idx), 'omitnan');

                    rows = [rows; ...
                        aa, ms, iv, th, ...
                        meanKI, stdKI, ...
                        meanRat, stdRat, medRat, madRat, ...
                        minRat, maxRat, spreadRat, ...
                        signStable, meanJIIabs, meanJ2chg]; %#ok<AGROW>
                end
            end
        end
    end

    Summary = array2table(rows, ...
        'VariableNames', { ...
            'a0_over_base', 'meshScale', 'variantID', 'thetaDeg', ...
            'KI_mean', 'KI_std', ...
            'KII_over_KI_mean', 'KII_over_KI_std', ...
            'KII_over_KI_median', 'KII_over_KI_MAD', ...
            'KII_over_KI_min', 'KII_over_KI_max', ...
            'KII_over_KI_spread', ...
            'signStable', ...
            'JII_over_absint_mean', 'J2_sign_changes_mean'});

    variantNames = unique(T(:, {'variantID','variant'}), 'rows');
    Summary.variant = strings(height(Summary),1);

    for i = 1:height(variantNames)
        Summary.variant(Summary.variantID == variantNames.variantID(i)) = variantNames.variant(i);
    end

    Summary = movevars(Summary, 'variant', 'After', 'variantID');
end


% ========================================================================
% Parity table from summary medians
% ========================================================================

function Parity = local_make_parity_table(Summary)

    if isempty(Summary)
        Parity = table();
        return;
    end

    rows = [];

    aVals = unique(Summary.a0_over_base);
    mVals = unique(Summary.meshScale);
    vVals = unique(Summary.variant);
    thetaAbsVals = unique(abs(Summary.thetaDeg));
    thetaAbsVals(thetaAbsVals == 0) = [];

    for ia = 1:numel(aVals)
        aa = aVals(ia);

        for im = 1:numel(mVals)
            ms = mVals(im);

            for iv = 1:numel(vVals)
                vv = vVals(iv);

                for it = 1:numel(thetaAbsVals)
                    th = thetaAbsVals(it);

                    idxP = abs(Summary.a0_over_base - aa) < 1e-12 ...
                        & abs(Summary.meshScale - ms) < 1e-12 ...
                        & Summary.variant == vv ...
                        & abs(Summary.thetaDeg - th) < 1e-12;

                    idxM = abs(Summary.a0_over_base - aa) < 1e-12 ...
                        & abs(Summary.meshScale - ms) < 1e-12 ...
                        & Summary.variant == vv ...
                        & abs(Summary.thetaDeg + th) < 1e-12;

                    if ~any(idxP) || ~any(idxM)
                        continue;
                    end

                    Kp = Summary.KII_over_KI_median(idxP);
                    Km = Summary.KII_over_KI_median(idxM);

                    Keven = 0.5*(Kp + Km);
                    Kodd  = 0.5*(Kp - Km);

                    parityRatio = abs(Keven) / max(abs(Kodd), eps);

                    spreadP = Summary.KII_over_KI_spread(idxP);
                    spreadM = Summary.KII_over_KI_spread(idxM);

                    rows = [rows; ...
                        aa, ms, iv, th, ...
                        Km, Kp, Keven, Kodd, parityRatio, ...
                        spreadM, spreadP]; %#ok<AGROW>
                end
            end
        end
    end

    Parity = array2table(rows, ...
        'VariableNames', { ...
            'a0_over_base', 'meshScale', 'variantID', 'thetaAbsDeg', ...
            'Kminus', 'Kplus', 'Keven', 'Kodd', 'absEven_over_absOdd', ...
            'spreadMinus', 'spreadPlus'});

    if isempty(Parity)
        return;
    end

    variantNames = unique(Summary(:, {'variantID','variant'}), 'rows');
    Parity.variant = strings(height(Parity),1);

    for i = 1:height(variantNames)
        Parity.variant(Parity.variantID == variantNames.variantID(i)) = variantNames.variant(i);
    end

    Parity = movevars(Parity, 'variant', 'After', 'variantID');
end


% ========================================================================
% Plot: KII/KI versus radius for each variant
% ========================================================================

function local_plot_variant_radius(T)

    good = T.valid & isfinite(T.KII_over_KI);

    if ~any(good)
        warning('No valid data for radius plot.');
        return;
    end

    variants = unique(T.variant(good));
    thetaVals = unique(T.thetaDeg(good));

    for iv = 1:numel(variants)
        vv = variants(iv);

        figure('Name', sprintf('Centered parity: radius sensitivity - %s', vv), ...
            'Color', 'w');
        clf;
        hold on;
        box on;
        grid on;

        for it = 1:numel(thetaVals)
            th = thetaVals(it);

            idx = good & T.variant == vv & abs(T.thetaDeg - th) < 1e-12;

            if ~any(idx)
                continue;
            end

            [rr, ord] = sort(T.rI_over_a0(idx));
            yy = T.KII_over_KI(idx);
            yy = yy(ord);

            plot(rr, yy, 'o-', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('\\theta = %.0f^\\circ', th));
        end

        yline(0, 'k--');

        xlabel('r_I/a_0');
        ylabel('K_{II}/K_I');
        title(sprintf('Centered-hole radius sensitivity, variant: %s', vv));
        legend('Location', 'best');
    end
end


% ========================================================================
% Plot: summary medians versus theta
% ========================================================================

function local_plot_variant_theta_summary(Summary)

    if isempty(Summary)
        warning('No summary data for theta plot.');
        return;
    end

    variants = unique(Summary.variant);

    figure('Name', 'Centered parity: median KII/KI versus theta', ...
        'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    for iv = 1:numel(variants)
        vv = variants(iv);

        idx = Summary.variant == vv;

        [th, ord] = sort(Summary.thetaDeg(idx));
        yy = Summary.KII_over_KI_median(idx);
        yy = yy(ord);

        plot(th, yy, 'o-', 'LineWidth', 1.4, ...
            'DisplayName', char(vv));
    end

    yline(0, 'k--');

    xlabel('\theta, deg');
    ylabel('median K_{II}/K_I');
    title('Centered-hole parity check');
    legend('Location', 'best');
end


% ========================================================================
% Plot: even and odd parts
% ========================================================================

function local_plot_parity(Parity)

    if isempty(Parity)
        warning('No parity data to plot.');
        return;
    end

    variants = unique(Parity.variant);

    figure('Name', 'Centered parity: even and odd parts', ...
        'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    for iv = 1:numel(variants)
        vv = variants(iv);

        idx = Parity.variant == vv;

        [th, ord] = sort(Parity.thetaAbsDeg(idx));
        Ke = Parity.Keven(idx);
        Ko = Parity.Kodd(idx);

        Ke = Ke(ord);
        Ko = Ko(ord);

        plot(th, Ke, 'o-', 'LineWidth', 1.3, ...
            'DisplayName', sprintf('%s: even', vv));
        plot(th, Ko, 's--', 'LineWidth', 1.3, ...
            'DisplayName', sprintf('%s: odd', vv));
    end

    yline(0, 'k--');

    xlabel('|\theta|, deg');
    ylabel('K_{II}/K_I decomposition');
    title('Even/odd decomposition of K_{II}/K_I');
    legend('Location', 'best');
end


% ========================================================================
% Stage-I plot
% ========================================================================

function local_plot_stage1_boundary_stress(B, I)

    figure('Name', 'Stage I: centered-hole boundary stress', 'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    plot(rad2deg(B.phi), B.sig_tt_eff, 'LineWidth', 1.2);
    plot(rad2deg(I.phi_star), I.sig_tt_unit, ...
        'ro', 'MarkerSize', 7, 'LineWidth', 1.5);

    xlabel('\phi, deg');
    ylabel('\sigma_{\theta\theta}^{eff} at unit load');
    title('Centered hole: tangential stress along hole boundary');
    xlim([0, 360]);
end


% ========================================================================
% Mesh scaling
% ========================================================================

function C = local_scale_mesh_config(C0, meshScale)

    C = C0;

    if isfield(C, 'mesh1')
        if isfield(C.mesh1, 'hmin') && ~isempty(C.mesh1.hmin)
            C.mesh1.hmin = meshScale * C0.mesh1.hmin;
        end
        if isfield(C.mesh1, 'hmax') && ~isempty(C.mesh1.hmax)
            C.mesh1.hmax = meshScale * C0.mesh1.hmax;
        end
        if isfield(C.mesh1, 'hhole') && ~isempty(C.mesh1.hhole)
            C.mesh1.hhole = meshScale * C0.mesh1.hhole;
        end
    end

    if isfield(C, 'mesh2')
        if isfield(C.mesh2, 'hmax') && ~isempty(C.mesh2.hmax)
            C.mesh2.hmax = meshScale * C0.mesh2.hmax;
        end
        if isfield(C.mesh2, 'hmin') && ~isempty(C.mesh2.hmin)
            C.mesh2.hmin = meshScale * C0.mesh2.hmin;
        end
        if isfield(C.mesh2, 'hhole') && ~isempty(C.mesh2.hhole)
            C.mesh2.hhole = meshScale * C0.mesh2.hhole;
        end
        if isfield(C.mesh2, 'hcrack') && ~isempty(C.mesh2.hcrack)
            C.mesh2.hcrack = meshScale * C0.mesh2.hcrack;
        end
    end
end


% ========================================================================
% Safe get field
% ========================================================================

function v = getfield_if_exists(S, field, default)

    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        v = S.(field);
    else
        v = default;
    end
end