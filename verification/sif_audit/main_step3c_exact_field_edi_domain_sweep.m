function Out = main_step3c_exact_field_edi_domain_sweep()
%MAIN_STEP3C_EXACT_FIELD_EDI_DOMAIN_SWEEP
% Step 3C of the asymmetric-mesh SIF audit.
%
% Purpose:
%   Determine whether the EDI annulus sensitivity observed on the numerical
%   FEM crack-tip field is intrinsic to the interaction-integral
%   implementation or is mainly caused by finite-element approximation of
%   the physical field.
%
% Strategy:
%   Reuse the exact Williams-field validator from Step 2. The production EDI
%   routine is evaluated on a fixed fine polar T6 mesh with known prescribed
%   SIFs while the EDI annulus is varied.
%
% Sweeps:
%   A) outer-radius sweep with fixed r_inner/r_outer = 0.20
%   B) inner-radius sweep with fixed r_outer = 0.12
%
% Exact fields already provided by validate_EDI_Williams_fields:
%   pure I       KI=1, KII=0
%   pure II      KI=0, KII=1
%   mixed I/II   KI=1, KII=0.35
%
% No physical FEM boundary-value problem is solved in this step.

    here = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(here));
    addpath(genpath(repoRoot));

    fprintf('\n============================================================\n');
    fprintf('SIF AUDIT STEP 3C: EXACT-FIELD EDI DOMAIN SWEEP\n');
    fprintf('============================================================\n');

    % Fixed exact-field interpolation mesh.
    Nr = 16;
    Nth = 128;
    rMeshInner = 0.005;
    rMeshOuter = 0.20;

    % ------------------------------------------------------------
    % A. Outer-radius sweep, fixed inner/outer ratio
    % ------------------------------------------------------------
    innerFactor = 0.20;
    outerList = [0.06 0.08 0.10 0.12 0.16];

    A = nan(numel(outerList), 10);

    for i = 1:numel(outerList)
        rout = outerList(i);
        rin = innerFactor*rout;

        R = validate_EDI_Williams_fields( ...
            'NrList', Nr, ...
            'NthList', Nth, ...
            'rMeshInner', rMeshInner, ...
            'rMeshOuter', rMeshOuter, ...
            'rInner', rin, ...
            'rOuter', rout, ...
            'Verbose', false, ...
            'AssertFine', false);

        T = R.table;

        pI = T(T.caseName=="pure_I",:);
        pII = T(T.caseName=="pure_II",:);
        mix = T(T.caseName=="mixed_I_II",:);

        A(i,:) = [ ...
            rout, rin, ...
            pI.KI_recovered_over_input, abs(pI.KII_recovered), ...
            pII.KII_recovered_over_input, abs(pII.KI_recovered), ...
            mix.KI_recovered_over_input, mix.KII_recovered_over_input, ...
            mix.KI_error_metric, mix.KII_error_metric];
    end

    Touter = array2table(A, 'VariableNames', { ...
        'r_outer','r_inner', ...
        'pureI_KI_ratio','pureI_abs_crossKII', ...
        'pureII_KII_ratio','pureII_abs_crossKI', ...
        'mixed_KI_ratio','mixed_KII_ratio', ...
        'mixed_KI_error','mixed_KII_error'});

    % ------------------------------------------------------------
    % B. Inner-radius sweep, fixed outer radius
    % ------------------------------------------------------------
    fixedOuter = 0.12;
    innerFactors = [0.10 0.20 0.30 0.40 0.50];

    B = nan(numel(innerFactors), 10);

    for i = 1:numel(innerFactors)
        fac = innerFactors(i);
        rin = fac*fixedOuter;

        R = validate_EDI_Williams_fields( ...
            'NrList', Nr, ...
            'NthList', Nth, ...
            'rMeshInner', rMeshInner, ...
            'rMeshOuter', rMeshOuter, ...
            'rInner', rin, ...
            'rOuter', fixedOuter, ...
            'Verbose', false, ...
            'AssertFine', false);

        T = R.table;

        pI = T(T.caseName=="pure_I",:);
        pII = T(T.caseName=="pure_II",:);
        mix = T(T.caseName=="mixed_I_II",:);

        B(i,:) = [ ...
            fac, rin, ...
            pI.KI_recovered_over_input, abs(pI.KII_recovered), ...
            pII.KII_recovered_over_input, abs(pII.KI_recovered), ...
            mix.KI_recovered_over_input, mix.KII_recovered_over_input, ...
            mix.KI_error_metric, mix.KII_error_metric];
    end

    Tinner = array2table(B, 'VariableNames', { ...
        'inner_over_outer','r_inner', ...
        'pureI_KI_ratio','pureI_abs_crossKII', ...
        'pureII_KII_ratio','pureII_abs_crossKI', ...
        'mixed_KI_ratio','mixed_KII_ratio', ...
        'mixed_KI_error','mixed_KII_error'});

    % ------------------------------------------------------------
    % Compact sensitivity summary
    % ------------------------------------------------------------
    S = struct();

    S.outerSweep = struct();
    S.outerSweep.pureI_ratio_range = range(Touter.pureI_KI_ratio);
    S.outerSweep.pureII_ratio_range = range(Touter.pureII_KII_ratio);
    S.outerSweep.mixed_KI_ratio_range = range(Touter.mixed_KI_ratio);
    S.outerSweep.mixed_KII_ratio_range = range(Touter.mixed_KII_ratio);

    S.innerSweep = struct();
    S.innerSweep.pureI_ratio_range = range(Tinner.pureI_KI_ratio);
    S.innerSweep.pureII_ratio_range = range(Tinner.pureII_KII_ratio);
    S.innerSweep.mixed_KI_ratio_range = range(Tinner.mixed_KI_ratio);
    S.innerSweep.mixed_KII_ratio_range = range(Tinner.mixed_KII_ratio);

    Out = struct();
    Out.outerSweep = Touter;
    Out.innerSweep = Tinner;
    Out.sensitivity = S;
    Out.settings = struct( ...
        'Nr',Nr, ...
        'Nth',Nth, ...
        'rMeshInner',rMeshInner, ...
        'rMeshOuter',rMeshOuter, ...
        'outerSweepInnerFactor',innerFactor, ...
        'innerSweepFixedOuter',fixedOuter);

    fprintf('\n--- Step 3C-A: exact field, outer-radius sweep ---\n');
    fprintf('fixed r_inner/r_outer = %.2f\n\n', innerFactor);
    disp(Touter);

    fprintf('\n--- Step 3C-B: exact field, inner-radius sweep ---\n');
    fprintf('fixed r_outer = %.4f\n\n', fixedOuter);
    disp(Tinner);

    fprintf('\nSensitivity ranges of recovered/input ratios:\n');
    fprintf('  outer sweep: pure I  = %.6e\n', S.outerSweep.pureI_ratio_range);
    fprintf('  outer sweep: pure II = %.6e\n', S.outerSweep.pureII_ratio_range);
    fprintf('  outer sweep: mixed KI  = %.6e\n', S.outerSweep.mixed_KI_ratio_range);
    fprintf('  outer sweep: mixed KII = %.6e\n', S.outerSweep.mixed_KII_ratio_range);
    fprintf('  inner sweep: pure I  = %.6e\n', S.innerSweep.pureI_ratio_range);
    fprintf('  inner sweep: pure II = %.6e\n', S.innerSweep.pureII_ratio_range);
    fprintf('  inner sweep: mixed KI  = %.6e\n', S.innerSweep.mixed_KI_ratio_range);
    fprintf('  inner sweep: mixed KII = %.6e\n', S.innerSweep.mixed_KII_ratio_range);

    fprintf('\nSTEP 3C completed.\n');
    fprintf(['Interpretation: if the exact-field ratios stay close to 1 and ', ...
        'vary only weakly with annulus choice, the much stronger KII ', ...
        'domain sensitivity seen in Step 3B is primarily a property of ', ...
        'the numerical FEM crack-tip field rather than of the EDI ', ...
        'formulation itself.\n']);
end
