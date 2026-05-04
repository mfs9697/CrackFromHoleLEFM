function Results = main_check_SIF_mesh_contour_sensitivity()
%MAIN_CHECK_SIF_MESH_CONTOUR_SENSITIVITY
% Diagnostic driver for checking the numerical stability of LEFM SIFs
% for the crack-from-hole problem.
%
% Purpose:
%   1) Solve the hole-only problem and determine the initiation point.
%   2) For selected appendix angles, build the collapsed cracked mesh.
%   3) For several mesh refinement levels and contour radii, compute KI,KII.
%   4) Assess whether KII/KI is stable enough to use KII(theta)=0.
%
% This driver intentionally does NOT run fzero for KII(theta)=0.
% It is a diagnostic check before using the local-symmetry condition.

clc;
close all;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('MAIN_CHECK_SIF_MESH_CONTOUR_SENSITIVITY\n');
fprintf('Mesh and contour-radius sensitivity of KI, KII\n');
fprintf('============================================================\n');

%% ============================================================
% 0. Configuration
%% ============================================================
C0 = cfg_hole_initiation();

% Plot controls
doPlotStage1 = true;
doPlotMeshes = false;  % keep false unless debugging one angle/mesh

%% ============================================================
% 1. Stage I: hole-only problem, once for the baseline configuration
%% ============================================================
fprintf('\n=== Stage I: hole-only solution and initiation point ===\n');

G  = geom_hole_only(C0);
S1 = solve_hole_only(C0, G, 'lambda', 1.0);
B  = sample_hole_boundary_stress(C0, G, S1);
I  = find_hole_initiation_point(C0, B);

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
% 2. Mesh + contour-radius / EDI-domain sensitivity loop
%% ============================================================
fprintf('\n=== Stage II: mesh, appendix-length, and SIF-domain sensitivity ===\n');

% ------------------------------------------------------------
% Study controls
% ------------------------------------------------------------

% SIF extraction method:
%   'circle2_debug'    : old circular mode-separation diagnostic
%   'interaction_EDI'  : new equivalent-domain interaction integral
sifMethod = 'interaction_EDI';

% For interaction_EDI:
% Existing rI/rI_over_a0 lists are interpreted as OUTER EDI-domain radii.
% The inner radius is chosen as:
%   r_inner = innerFactorEDI * r_outer
innerFactorEDI = 0.1;

% Appendix lengths to test
a0_list = [4*C0.a0];

% Mesh refinement
meshScaleList = [0.25];

% Angles
thetaDegProbe = -1:0.5:3;
thetaProbe    = deg2rad(thetaDegProbe);

% Radius mode:
%   'relative' : r_outer = rI_over_a0 * a0
%   'absolute' : r_outer prescribed directly
radiusMode = 'relative';

% Relative outer radii
rI_over_a0_list = 0.5:0.1:0.7;

% Absolute outer radii, used only if radiusMode = 'absolute'
rI_abs_list = [];

% Debug controls
debugVerbose = false;
debugPlot    = false;

% Mesh plotting
doPlotMeshes = false;

% ------------------------------------------------------------
% Storage
% ------------------------------------------------------------
rows = [];

DebugStore = struct();
DebugStore.items = {};
DebugStore.count = 0;

% ------------------------------------------------------------
% Main loop
% ------------------------------------------------------------
for ia = 1:numel(a0_list)

    a0 = a0_list(ia);
    a0_over_base = a0 / C0.a0;

    fprintf('\n============================================================\n');
    fprintf('Appendix length a0 = %.8e, a0/a0_base = %.3f\n', ...
        a0, a0_over_base);
    fprintf('============================================================\n');

    for im = 1:numel(meshScaleList)

        meshScale = meshScaleList(im);

        C = local_scale_mesh_config(C0, meshScale);
        C.a0 = a0;

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Mesh scale = %.3f\n', meshScale);
        fprintf('SIF method = %s\n', sifMethod);
        fprintf('------------------------------------------------------------\n');

        for it = 1:numel(thetaProbe)

            theta = thetaProbe(it);
            thetaDeg = thetaDegProbe(it);

            fprintf('\n--- a0/a0_base = %.3f | meshScale = %.3f | theta = %+8.3f deg ---\n', ...
                a0_over_base, meshScale, thetaDeg);

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

                % ----------------------------------------------------
                % Choose outer-radius list for this case
                % ----------------------------------------------------
                switch lower(strtrim(radiusMode))
                    case 'relative'
                        rfac_list = rI_over_a0_list(:).';
                        rI_list = rfac_list * G2.crack.a0;

                    case 'absolute'
                        rI_list = rI_abs_list(:).';
                        rfac_list = rI_list / G2.crack.a0;

                    otherwise
                        error('Unknown radiusMode = %s. Use relative or absolute.', radiusMode);
                end

                for ir = 1:numel(rI_list)

                    rI = rI_list(ir);
                    rfac = rfac_list(ir);

                    % For circle2_debug:
                    %   rI = circular contour radius.
                    %
                    % For interaction_EDI:
                    %   rI = r_outer, the outer EDI-domain radius.
                    r_outer = rI;
                    r_inner = innerFactorEDI * r_outer;

                    try
                        % ------------------------------------------------
                        % SIF calculation
                        % ------------------------------------------------
                        switch lower(strtrim(sifMethod))

                            case 'circle2_debug'
                                [KI, KII, AuxSIF] = local_compute_SIF_for_stage2_debug( ...
                                    C, G2, Mc, S2, rI, ...
                                    'verbose', debugVerbose, ...
                                    'plot', debugPlot);

                            case 'interaction_edi'
                                [KI, KII, AuxSIF] = compute_SIF_for_stage2_interaction( ...
                                    C, G2, Mc, S2, ...
                                    'r_inner', r_inner, ...
                                    'r_outer', r_outer, ...
                                    'Verbose', debugVerbose);

                            otherwise
                                error('Unknown sifMethod = %s.', sifMethod);
                        end

                        KII_over_KI = KII / KI;

                        % ------------------------------------------------
                        % Circular-method diagnostics
                        % ------------------------------------------------
                        diag_JII_over_absint = nan;
                        diag_J2_sign_changes = nan;
                        diag_nearP_1e4 = nan;
                        diag_nearQ_1e4 = nan;
                        diag_minBaryP = nan;
                        diag_minBaryQ = nan;

                        if isstruct(AuxSIF) && isfield(AuxSIF, 'Dbg')
                            Dbg = AuxSIF.Dbg;

                            if isfield(Dbg, 'diagnostics')
                                DD = Dbg.diagnostics;

                                if isfield(DD, 'JII_over_absint')
                                    diag_JII_over_absint = DD.JII_over_absint;
                                end
                                if isfield(DD, 'J2_sign_changes')
                                    diag_J2_sign_changes = DD.J2_sign_changes;
                                end
                                if isfield(DD, 'num_near_edge_P_1e4')
                                    diag_nearP_1e4 = DD.num_near_edge_P_1e4;
                                end
                                if isfield(DD, 'num_near_edge_Q_1e4')
                                    diag_nearQ_1e4 = DD.num_near_edge_Q_1e4;
                                end
                                if isfield(DD, 'min_baryP')
                                    diag_minBaryP = DD.min_baryP;
                                end
                                if isfield(DD, 'min_baryQ')
                                    diag_minBaryQ = DD.min_baryQ;
                                end
                            end
                        end

                        % ------------------------------------------------
                        % Interaction-EDI diagnostics
                        % ------------------------------------------------
                        diag_I_modeI = nan;
                        diag_I_modeII = nan;
                        diag_nGP_used = nan;
                        diag_nElem_used = nan;
                        diag_r_inner_over_a0 = nan;
                        diag_r_outer_over_a0 = nan;

                        if isstruct(AuxSIF)
                            if isfield(AuxSIF, 'I_modeI')
                                diag_I_modeI = AuxSIF.I_modeI;
                            end
                            if isfield(AuxSIF, 'I_modeII')
                                diag_I_modeII = AuxSIF.I_modeII;
                            end
                            if isfield(AuxSIF, 'nGP_used')
                                diag_nGP_used = AuxSIF.nGP_used;
                            end
                            if isfield(AuxSIF, 'nElem_used')
                                diag_nElem_used = AuxSIF.nElem_used;
                            end
                            if isfield(AuxSIF, 'r_inner_over_a0')
                                diag_r_inner_over_a0 = AuxSIF.r_inner_over_a0;
                            end
                            if isfield(AuxSIF, 'r_outer_over_a0')
                                diag_r_outer_over_a0 = AuxSIF.r_outer_over_a0;
                            end
                        end

                        rows = [rows; ...
                            a0, a0_over_base, ...
                            meshScale, thetaDeg, ...
                            rfac, rI, r_inner, r_outer, ...
                            KI, KII, KII_over_KI, ...
                            nNode3, nElem3, nNode6, nElem6, ...
                            diag_JII_over_absint, ...
                            diag_J2_sign_changes, ...
                            diag_nearP_1e4, diag_nearQ_1e4, ...
                            diag_minBaryP, diag_minBaryQ, ...
                            diag_I_modeI, diag_I_modeII, ...
                            diag_nGP_used, diag_nElem_used, ...
                            diag_r_inner_over_a0, diag_r_outer_over_a0, ...
                            true]; %#ok<AGROW>

                        fprintf(['  r_outer/a0 = %.4f | r_outer = %.4e | ', ...
                                 'r_inner = %.4e | KI = %.8e | KII = %+ .8e | ', ...
                                 'KII/KI = %+ .8e'], ...
                            rfac, r_outer, r_inner, KI, KII, KII_over_KI);

                        if strcmpi(sifMethod, 'circle2_debug') && isfinite(diag_JII_over_absint)
                            fprintf(' | JII/int|J2| = %+ .3e | J2 changes = %g', ...
                                diag_JII_over_absint, diag_J2_sign_changes);
                        end

                        if strcmpi(sifMethod, 'interaction_EDI') && isfinite(diag_nGP_used)
                            fprintf(' | I_I = %+ .3e | I_II = %+ .3e | GP = %g | Elem = %g', ...
                                diag_I_modeI, diag_I_modeII, diag_nGP_used, diag_nElem_used);
                        end

                        fprintf('\n');

                        % Store selected debug data only when available.
                        if isstruct(AuxSIF)
                            keepThisDebug = ...
                                abs(thetaDeg) < 1e-12 && ...
                                (abs(rI - rI_list(1)) < 1e-14 || abs(rI - rI_list(end)) < 1e-14);

                            if keepThisDebug
                                DebugStore.count = DebugStore.count + 1;
                                DebugStore.items{DebugStore.count} = struct( ...
                                    'sifMethod', sifMethod, ...
                                    'a0', a0, ...
                                    'a0_over_base', a0_over_base, ...
                                    'meshScale', meshScale, ...
                                    'thetaDeg', thetaDeg, ...
                                    'rI_over_a0', rfac, ...
                                    'rI', rI, ...
                                    'r_inner', r_inner, ...
                                    'r_outer', r_outer, ...
                                    'AuxSIF', AuxSIF);
                            end
                        end

                    catch MEinner
                        rows = [rows; ...
                            a0, a0_over_base, ...
                            meshScale, thetaDeg, ...
                            rfac, rI, r_inner, r_outer, ...
                            nan, nan, nan, ...
                            nNode3, nElem3, nNode6, nElem6, ...
                            nan, nan, nan, nan, nan, nan, ...
                            nan, nan, nan, nan, nan, nan, ...
                            false]; %#ok<AGROW>

                        fprintf('  r_outer/a0 = %.4f | r_outer = %.4e FAILED: %s\n', ...
                            rfac, r_outer, MEinner.message);
                    end
                end

            catch MEouter
                fprintf('  theta failed before radius/domain loop: %s\n', MEouter.message);

                switch lower(strtrim(radiusMode))
                    case 'relative'
                        rfac_list = rI_over_a0_list(:).';
                        rI_list = rfac_list * a0;
                    case 'absolute'
                        rI_list = rI_abs_list(:).';
                        rfac_list = rI_list / a0;
                end

                for ir = 1:numel(rI_list)
                    rI = rI_list(ir);
                    rfac = rfac_list(ir);

                    r_outer = rI;
                    r_inner = innerFactorEDI * r_outer;

                    rows = [rows; ...
                        a0, a0_over_base, ...
                        meshScale, thetaDeg, ...
                        rfac, rI, r_inner, r_outer, ...
                        nan, nan, nan, ...
                        nan, nan, nan, nan, ...
                        nan, nan, nan, nan, nan, nan, ...
                        nan, nan, nan, nan, nan, nan, ...
                        false]; %#ok<AGROW>
                end
            end
        end
    end
end

%% ============================================================
% 3. Convert to table and print
%% ============================================================

T = array2table(rows, ...
    'VariableNames', { ...
    'a0', 'a0_over_base', ...
    'meshScale', 'thetaDeg', ...
    'rI_over_a0', 'rI', 'r_inner', 'r_outer', ...
    'KI', 'KII', 'KII_over_KI', ...
    'nNode3', 'nElem3', 'nNode6', 'nElem6', ...
    'JII_over_absint', ...
    'J2_sign_changes', ...
    'nearP_1e4', 'nearQ_1e4', ...
    'minBaryP', 'minBaryQ', ...
    'I_modeI', 'I_modeII', ...
    'nGP_used', 'nElem_used', ...
    'r_inner_over_a0', 'r_outer_over_a0', ...
    'valid'});

fprintf('\n=== Full mesh-appendix-SIF-domain sensitivity table ===\n');
disp(T);
%% ============================================================
% 4. Summary statistics
%% ============================================================

Summary = local_make_summary_table_with_a0(T);

fprintf('\n=== Summary over contour radii for each a0, mesh scale, and theta ===\n');
disp(Summary);

%% ============================================================
% 5. Plots
%% ============================================================

local_plot_KII_over_KI_vs_radius_with_a0(T);
local_plot_KII_over_KI_vs_theta_with_a0(T);
local_plot_summary_spread_with_a0(Summary);

%% ============================================================
% 6. Store output
%% ============================================================

Results.Table = T;
Results.Summary = Summary;
Results.DebugStore = DebugStore;
Results.a0_list = a0_list;
Results.meshScaleList = meshScaleList;
Results.thetaDegProbe = thetaDegProbe;
Results.radiusMode = radiusMode;

if strcmpi(radiusMode, 'relative')
    Results.rI_over_a0_list = rI_over_a0_list;
else
    Results.rI_abs_list = rI_abs_list;
end
end


% ========================================================================
% Scale mesh-related fields in configuration
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
    if isfield(C.mesh1, 'hgrad') && ~isempty(C.mesh1.hgrad)
        % Usually keep hgrad unchanged. Uncomment only if needed.
        % C.mesh1.hgrad = C0.mesh1.hgrad;
    end
end

if isfield(C, 'mesh2')
    if isfield(C.mesh2, 'hmax') && ~isempty(C.mesh2.hmax)
        C.mesh2.hmax = meshScale * C0.mesh2.hmax;
    end
    if isfield(C.mesh2, 'hmin') && ~isempty(C.mesh2.hmin)
        C.mesh2.hmin = meshScale * C0.mesh2.hmin;
    end
    if isfield(C.mesh2, 'hcrack') && ~isempty(C.mesh2.hcrack)
        C.mesh2.hcrack = meshScale * C0.mesh2.hcrack;
    end
    if isfield(C.mesh2, 'hhole') && ~isempty(C.mesh2.hhole)
        C.mesh2.hhole = meshScale * C0.mesh2.hhole;
    end

    % Do not scale the channel width by default.
    % It is a geometric regularization parameter, not just a mesh size.
    % If needed later, test scaling C.mesh2.chw separately.
end
end


% ========================================================================
% Summary table
% ========================================================================

function Summary = local_make_summary_table(T)

good = T.valid & isfinite(T.KI) & isfinite(T.KII_over_KI);

meshVals  = unique(T.meshScale(good));
thetaVals = unique(T.thetaDeg(good));

rows = [];

for im = 1:numel(meshVals)
    ms = meshVals(im);

    for it = 1:numel(thetaVals)
        th = thetaVals(it);

        idx = good ...
            & abs(T.meshScale - ms) < 1e-12 ...
            & abs(T.thetaDeg - th) < 1e-12;

        if ~any(idx)
            continue;
        end

        valsKI  = T.KI(idx);
        valsRat = T.KII_over_KI(idx);

        meanKI = mean(valsKI, 'omitnan');
        stdKI  = std(valsKI,  'omitnan');

        meanRat = mean(valsRat, 'omitnan');
        stdRat  = std(valsRat,  'omitnan');
        minRat  = min(valsRat);
        maxRat  = max(valsRat);
        spreadRat = maxRat - minRat;

        signStable = all(valsRat > 0) || all(valsRat < 0);

        rows = [rows; ...
            ms, th, ...
            meanKI, stdKI, ...
            meanRat, stdRat, minRat, maxRat, spreadRat, ...
            signStable]; %#ok<AGROW>
    end
end

Summary = array2table(rows, ...
    'VariableNames', { ...
    'meshScale', 'thetaDeg', ...
    'KI_mean', 'KI_std', ...
    'KII_over_KI_mean', 'KII_over_KI_std', ...
    'KII_over_KI_min', 'KII_over_KI_max', ...
    'KII_over_KI_spread', ...
    'signStable'});
end


% ========================================================================
% Plot Stage-I boundary stress
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
% Plot KII/KI versus rI/a0 for each mesh scale
% ========================================================================

function local_plot_KII_over_KI_vs_radius(T)

good = T.valid & isfinite(T.KII_over_KI);

if ~any(good)
    warning('No valid data for radius-sensitivity plot.');
    return;
end

meshVals  = unique(T.meshScale(good));
thetaVals = unique(T.thetaDeg(good));

for im = 1:numel(meshVals)
    ms = meshVals(im);

    figure('Name', sprintf('KII/KI vs rI/a0, mesh scale %.2f', ms), ...
        'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    for it = 1:numel(thetaVals)
        th = thetaVals(it);

        idx = good ...
            & abs(T.meshScale - ms) < 1e-12 ...
            & abs(T.thetaDeg - th) < 1e-12;

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
    title(sprintf('Contour-radius sensitivity, mesh scale %.2f', ms));
    legend('Location', 'best');
end
end


% ========================================================================
% Plot KII/KI versus theta for each contour radius and mesh scale
% ========================================================================

function local_plot_KII_over_KI_vs_theta(T)

good = T.valid & isfinite(T.KII_over_KI);

if ~any(good)
    warning('No valid data for theta-sensitivity plot.');
    return;
end

meshVals = unique(T.meshScale(good));
rVals    = unique(T.rI_over_a0(good));

for im = 1:numel(meshVals)
    ms = meshVals(im);

    figure('Name', sprintf('KII/KI vs theta, mesh scale %.2f', ms), ...
        'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    for ir = 1:numel(rVals)
        rr = rVals(ir);

        idx = good ...
            & abs(T.meshScale - ms) < 1e-12 ...
            & abs(T.rI_over_a0 - rr) < 1e-12;

        if ~any(idx)
            continue;
        end

        [th, ord] = sort(T.thetaDeg(idx));
        yy = T.KII_over_KI(idx);
        yy = yy(ord);

        plot(th, yy, 'o-', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('r_I/a_0 = %.2f', rr));
    end

    yline(0, 'k--');

    xlabel('\theta, deg');
    ylabel('K_{II}/K_I');
    title(sprintf('Angular dependence, mesh scale %.2f', ms));
    legend('Location', 'best');
end
end


% ========================================================================
% Plot spread of KII/KI over contour radii
% ========================================================================

function local_plot_summary_spread(Summary)

if isempty(Summary)
    warning('No summary data to plot.');
    return;
end

meshVals = unique(Summary.meshScale);

figure('Name', 'Spread of KII/KI over contour radii', 'Color', 'w');
clf;
hold on;
box on;
grid on;

for im = 1:numel(meshVals)
    ms = meshVals(im);

    idx = abs(Summary.meshScale - ms) < 1e-12;

    [th, ord] = sort(Summary.thetaDeg(idx));
    spread = Summary.KII_over_KI_spread(idx);
    spread = spread(ord);

    plot(th, spread, 'o-', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('mesh scale %.2f', ms));
end

xlabel('\theta, deg');
ylabel('spread of K_{II}/K_I over r_I/a_0');
title('Contour sensitivity indicator');
legend('Location', 'best');
end

function [KI, KII, Aux] = local_compute_SIF_for_stage2_debug(C, G2, Mc, S2, rI, varargin)
%LOCAL_COMPUTE_SIF_FOR_STAGE2_DEBUG
% Adapter around SIF_LEFM_circle2_debug for the Stage-II hole-crack workflow.

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
        error('local_compute_SIF_for_stage2_debug:MissingDmat', ...
            'Material struct must contain Dmat or D.');
    end
end

if isfield(Mc, 'crack') && isfield(Mc.crack, 'Pmid') && ~isempty(Mc.crack.Pmid)
    V = Mc.crack.Pmid;
else
    V = G2.crack.polyline;
end

[KI, KII, Dbg] = SIF_LEFM_circle2_debug( ...
    meshSIF, S2.U, V, matSIF, rI, ...
    'nthet', ip.Results.nthet, ...
    'plot', logical(ip.Results.plot), ...
    'verbose', logical(ip.Results.verbose));

Aux = struct();
Aux.meshSIF = meshSIF;
Aux.matSIF = matSIF;
Aux.V = V;
Aux.rI = rI;
Aux.Dbg = Dbg;
end

function Summary = local_make_summary_table_with_a0(T)

good = T.valid & isfinite(T.KI) & isfinite(T.KII_over_KI);

a0Vals    = unique(T.a0_over_base(good));
meshVals  = unique(T.meshScale(good));
thetaVals = unique(T.thetaDeg(good));

rows = [];

for ia = 1:numel(a0Vals)
    aa = a0Vals(ia);

    for im = 1:numel(meshVals)
        ms = meshVals(im);

        for it = 1:numel(thetaVals)
            th = thetaVals(it);

            idx = good ...
                & abs(T.a0_over_base - aa) < 1e-12 ...
                & abs(T.meshScale - ms) < 1e-12 ...
                & abs(T.thetaDeg - th) < 1e-12;

            if ~any(idx)
                continue;
            end

            valsKI  = T.KI(idx);
            valsRat = T.KII_over_KI(idx);

            meanKI = mean(valsKI, 'omitnan');
            stdKI  = std(valsKI,  'omitnan');

            meanRat = mean(valsRat, 'omitnan');
            stdRat  = std(valsRat,  'omitnan');
            medRat  = median(valsRat, 'omitnan');
            madRat  = median(abs(valsRat - medRat), 'omitnan');

            minRat  = min(valsRat);
            maxRat  = max(valsRat);
            spreadRat = maxRat - minRat;

            signStable = all(valsRat > 0) || all(valsRat < 0);

            valsJIIabs = T.JII_over_absint(idx);
            meanJIIabs = mean(valsJIIabs, 'omitnan');

            valsJ2chg = T.J2_sign_changes(idx);
            meanJ2chg = mean(valsJ2chg, 'omitnan');

            rows = [rows; ...
                aa, ms, th, ...
                meanKI, stdKI, ...
                meanRat, stdRat, medRat, madRat, ...
                minRat, maxRat, spreadRat, ...
                signStable, ...
                meanJIIabs, meanJ2chg]; %#ok<AGROW>
        end
    end
end

Summary = array2table(rows, ...
    'VariableNames', { ...
    'a0_over_base', 'meshScale', 'thetaDeg', ...
    'KI_mean', 'KI_std', ...
    'KII_over_KI_mean', 'KII_over_KI_std', ...
    'KII_over_KI_median', 'KII_over_KI_MAD', ...
    'KII_over_KI_min', 'KII_over_KI_max', ...
    'KII_over_KI_spread', ...
    'signStable', ...
    'JII_over_absint_mean', 'J2_sign_changes_mean'});
end
function local_plot_KII_over_KI_vs_radius_with_a0(T)

good = T.valid & isfinite(T.KII_over_KI);

if ~any(good)
    warning('No valid data for radius-sensitivity plot.');
    return;
end

a0Vals = unique(T.a0_over_base(good));
meshVals = unique(T.meshScale(good));
thetaVals = unique(T.thetaDeg(good));

for ia = 1:numel(a0Vals)
    aa = a0Vals(ia);

    for im = 1:numel(meshVals)
        ms = meshVals(im);

        figure('Name', sprintf('KII/KI vs radius, a0/base %.2f, mesh %.2f', aa, ms), ...
            'Color', 'w');
        clf;
        hold on;
        box on;
        grid on;

        for it = 1:numel(thetaVals)
            th = thetaVals(it);

            idx = good ...
                & abs(T.a0_over_base - aa) < 1e-12 ...
                & abs(T.meshScale - ms) < 1e-12 ...
                & abs(T.thetaDeg - th) < 1e-12;

            if ~any(idx)
                continue;
            end

            [rr, ord] = sort(T.rI(idx));
            yy = T.KII_over_KI(idx);
            yy = yy(ord);

            plot(rr, yy, 'o-', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('\\theta = %.1f^\\circ', th));
        end

        yline(0, 'k--');

        xlabel('r_I');
        ylabel('K_{II}/K_I');
        title(sprintf('Radius sensitivity, a_0/a_{0,base}=%.2f, mesh scale %.2f', aa, ms));
        legend('Location', 'best');
    end
end
end

function local_plot_KII_over_KI_vs_theta_with_a0(T)

good = T.valid & isfinite(T.KII_over_KI);

if ~any(good)
    warning('No valid data for theta-sensitivity plot.');
    return;
end

a0Vals = unique(T.a0_over_base(good));
meshVals = unique(T.meshScale(good));
rVals = unique(T.rI(good));

for ia = 1:numel(a0Vals)
    aa = a0Vals(ia);

    for im = 1:numel(meshVals)
        ms = meshVals(im);

        figure('Name', sprintf('KII/KI vs theta, a0/base %.2f, mesh %.2f', aa, ms), ...
            'Color', 'w');
        clf;
        hold on;
        box on;
        grid on;

        for ir = 1:numel(rVals)
            rr = rVals(ir);

            idx = good ...
                & abs(T.a0_over_base - aa) < 1e-12 ...
                & abs(T.meshScale - ms) < 1e-12 ...
                & abs(T.rI - rr) < 1e-12;

            if ~any(idx)
                continue;
            end

            [th, ord] = sort(T.thetaDeg(idx));
            yy = T.KII_over_KI(idx);
            yy = yy(ord);

            plot(th, yy, 'o-', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('r_I = %.4g', rr));
        end

        yline(0, 'k--');

        xlabel('\theta, deg');
        ylabel('K_{II}/K_I');
        title(sprintf('Angular dependence, a_0/a_{0,base}=%.2f, mesh scale %.2f', aa, ms));
        legend('Location', 'best');
    end
end
end
function local_plot_summary_spread_with_a0(Summary)

if isempty(Summary)
    warning('No summary data to plot.');
    return;
end

a0Vals = unique(Summary.a0_over_base);
meshVals = unique(Summary.meshScale);

for im = 1:numel(meshVals)
    ms = meshVals(im);

    figure('Name', sprintf('Spread of KII/KI, mesh %.2f', ms), ...
        'Color', 'w');
    clf;
    hold on;
    box on;
    grid on;

    for ia = 1:numel(a0Vals)
        aa = a0Vals(ia);

        idx = abs(Summary.meshScale - ms) < 1e-12 ...
            & abs(Summary.a0_over_base - aa) < 1e-12;

        if ~any(idx)
            continue;
        end

        [th, ord] = sort(Summary.thetaDeg(idx));
        spread = Summary.KII_over_KI_spread(idx);
        spread = spread(ord);

        plot(th, spread, 'o-', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('a_0/a_{0,base}=%.2f', aa));
    end

    xlabel('\theta, deg');
    ylabel('spread of K_{II}/K_I over contour radii');
    title(sprintf('Contour sensitivity indicator, mesh scale %.2f', ms));
    legend('Location', 'best');
end
end