function Results = main_stage1_centered_hole_stress_check()
%MAIN_STAGE1_CENTERED_HOLE_STRESS_CHECK
% Stage-I diagnostic driver for a centered circular hole.
%
% Purpose:
%   1) Put the circular hole at the plate center.
%   2) Check the boundary tangential-stress distribution.
%   3) Check stress-concentration convergence under iterative local
%      refinement near the symmetry-equivalent maxima.

    clc;
    close all;

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(repoRoot));

    fprintf('\n============================================================\n');
    fprintf('MAIN_STAGE1_CENTERED_HOLE_STRESS_CHECK\n');
    fprintf('Centered-hole Stage-I stress-concentration diagnostic\n');
    fprintf('============================================================\n');

    %% Configuration
    C = cfg_hole_initiation();

    % Centered-hole benchmark
    C.hole.center = [0.5*C.A, 0.0];

    % Symmetry-friendly polygon and sampling
    C.hole.npoly = 360;
    C.mesh1.hmin  = 2*pi*C.hole.r / C.hole.npoly;
    C.mesh1.hhole = C.mesh1.hmin;
    C.mesh1.hmax  = 10*C.mesh1.hmin;
    C.stage1.nphi = 1440;
    C.holes = {C.hole};

    % Iterative local-refinement controls
    maxIter = 5;
    tolKtRel = 5e-3;
    baseHalfWidth = deg2rad(10);
    minHalfWidth = deg2rad(3);
    phiExpectedSet = [0, pi];
    exportOutputs = true;
    outputDir = fullfile(fileparts(mfilename('fullpath')), ...
        'outputs', 'centered_hole_iterative_refinement');

    rows = [];

    CaseStore = struct();
    CaseStore.items = {};
    CaseStore.count = 0;

    prevKt = NaN;

    for iter = 1:maxIter
        if iter == 1
            C.mesh1 = rmfield_safe(C.mesh1, { ...
                'refineHoleAngles', 'refineAngleHalfWidth', 'hrefine'});
            halfWidth = NaN;
            hrefine = NaN;
        else
            halfWidth = max(minHalfWidth, baseHalfWidth/(iter - 1));
            hrefine = C.mesh1.hhole/iter;

            % The centered benchmark has two symmetry-equivalent peaks.
            % Refining both avoids chasing one numerical tie-breaker.
            C.mesh1.refineHoleAngles = phiExpectedSet;
            C.mesh1.refineAngleHalfWidth = halfWidth;
            C.mesh1.hrefine = hrefine;
        end

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Stage-I local-refinement iteration %d/%d\n', iter, maxIter);
        if iter == 1
            fprintf('  refinement: off\n');
        else
            fprintf('  refinement: angles = [0, 180] deg, half-width = %.3f deg, hrefine = %.6e\n', ...
                rad2deg(halfWidth), hrefine);
        end
        fprintf('------------------------------------------------------------\n');

        G  = geom_hole_only(C);
        S1 = solve_hole_only(C, G, 'lambda', 1.0);
        B  = sample_hole_boundary_stress(C, G, S1);
        I  = find_hole_initiation_point(C, B);

        sig0 = C.load.sig0;
        Kt = I.sig_tt_pos_unit / sig0;
        if isfinite(prevKt)
            dKtRel = abs(Kt - prevKt) / max(abs(prevKt), eps);
        else
            dKtRel = NaN;
        end

        dphiSym = local_min_symmetry_angle_diff(I.phi_star, phiExpectedSet);

        nRefEdges = 0;
        nRefVertices = 0;
        refineActive = false;
        if isfield(G, 'meta') && isfield(G.meta, 'refinement')
            refineActive = G.meta.refinement.active;
            nRefEdges = numel(G.meta.refinement.edgeIDs);
            nRefVertices = numel(G.meta.refinement.vertexIDs);
        end

        rows = [rows; ...
            iter, ...
            double(refineActive), ...
            size(G.p,1), size(G.t,1), ...
            C.hole.npoly, C.stage1.nphi, ...
            hrefine, rad2deg(halfWidth), ...
            nRefEdges, nRefVertices, ...
            I.phi_star, rad2deg(I.phi_star), rad2deg(dphiSym), ...
            I.sig_tt_pos_unit, Kt, dKtRel, ...
            I.lambda_ini, B.eps_shift]; %#ok<AGROW>

        fprintf('  peak phi       = %.8f rad = %.5f deg\n', ...
            I.phi_star, rad2deg(I.phi_star));
        fprintf('  symmetry error = %.5f deg\n', rad2deg(dphiSym));
        fprintf('  Kt             = %.8f\n', Kt);
        if isfinite(dKtRel)
            fprintf('  relative dKt   = %.6e\n', dKtRel);
        end
        if refineActive
            fprintf('  refined IDs    = %d edges, %d vertices\n', nRefEdges, nRefVertices);
        end

        exportFiles = struct();
        if exportOutputs
            exportFiles = stage1_export_iteration_outputs(outputDir, iter, C, G, S1, B, I);
            fprintf('  outputs         = %s\n', outputDir);
            fprintf('  contour plot    = %s\n', exportFiles.holeContourPng);
            fprintf('  contour csv     = %s\n', exportFiles.holeContourCsv);
        end

        CaseStore.count = CaseStore.count + 1;
        CaseStore.items{CaseStore.count} = struct( ...
            'iter', iter, ...
            'C', C, 'G', G, 'S1', S1, 'B', B, 'I', I, ...
            'Kt', Kt, 'dKtRel', dKtRel, ...
            'dphiSym', dphiSym, ...
            'exportFiles', exportFiles);

        if iter > 1 && isfinite(dKtRel) && dKtRel < tolKtRel
            fprintf('  convergence reached: relative dKt < %.3e\n', tolKtRel);
            break;
        end

        prevKt = Kt;
    end

    T = array2table(rows, ...
        'VariableNames', { ...
            'iter', 'refineActive', ...
            'nNodeT3', 'nElemT3', ...
            'npoly', 'nphi', ...
            'hrefine', 'refineHalfWidthDeg', ...
            'nRefEdges', 'nRefVertices', ...
            'phi_peak_rad', 'phi_peak_deg', 'phi_error_to_symmetry_deg', ...
            'sig_tt_peak', 'Kt', 'dKt_rel', ...
            'lambda_ini', 'eps_shift'});

    fprintf('\n=== Stage-I centered-hole iterative-refinement table ===\n');
    disp(T);

    Results = struct();
    Results.Table = T;
    Results.CaseStore = CaseStore;
    Results.phiExpectedSet = phiExpectedSet;
    Results.tolKtRel = tolKtRel;
    Results.outputDir = outputDir;
    Results.PDEBenchmark = [];

    if CaseStore.count > 0
        last = CaseStore.items{CaseStore.count};
        try
            fprintf('\n=== PDE Toolbox circular-hole benchmark ===\n');
            P = stage1_pde_toolbox_circular_hole_benchmark(last.C, ...
                'Hmax', last.C.mesh1.hmax, ...
                'Hmin', last.C.mesh1.hmin, ...
                'Hgrad', last.C.mesh1.hgrad, ...
                'Nphi', last.C.stage1.nphi, ...
                'EpsShift', last.B.eps_shift);

            Results.PDEBenchmark = P;

            fprintf('  analysis type = %s\n', P.analysisType);
            fprintf('  peak phi      = %.8f rad = %.5f deg\n', ...
                P.phi_peak, P.phi_peak_deg);
            fprintf('  Kt benchmark  = %.8f\n', P.Kt);
            fprintf('  Kt custom     = %.8f\n', last.Kt);
            fprintf('  rel. diff     = %.6e\n', ...
                abs(last.Kt - P.Kt) / max(abs(P.Kt), eps));
        catch ME
            warning('main_stage1_centered_hole_stress_check:PDEBenchmarkFailed', ...
                'PDE Toolbox benchmark failed: %s', ME.message);
        end
    end

    % Plot the final iteration
    if CaseStore.count > 0
        last = CaseStore.items{CaseStore.count};
        Iexpected = stage1_force_initiation_phi(last.C, last.B, nearest_expected_phi(last.I.phi_star, phiExpectedSet));
        stage1_plot_boundary_stress(last.B, last.I, Iexpected);
    end
end


function Cmesh = rmfield_safe(Cmesh, fields)
%RMFIELD_SAFE Remove fields that may or may not exist.

    for k = 1:numel(fields)
        if isfield(Cmesh, fields{k})
            Cmesh = rmfield(Cmesh, fields{k});
        end
    end
end


function dphi = local_min_symmetry_angle_diff(phi, phiSet)
%LOCAL_MIN_SYMMETRY_ANGLE_DIFF Distance to nearest symmetry-equivalent angle.

    d = angle(exp(1i*(phi - phiSet(:))));
    dphi = min(abs(d));
end


function phi0 = nearest_expected_phi(phi, phiSet)
%NEAREST_EXPECTED_PHI Expected angle nearest to phi.

    d = angle(exp(1i*(phi - phiSet(:))));
    [~, idx] = min(abs(d));
    phi0 = phiSet(idx);
end
