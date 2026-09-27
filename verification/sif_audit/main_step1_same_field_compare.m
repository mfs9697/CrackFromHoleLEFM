function Results = main_step1_same_field_compare()
%MAIN_STEP1_SAME_FIELD_COMPARE
% First executable gate of the SIF asymmetric-mesh audit.
%
% It builds the two-leg Crack-Path control case, solves the elastic FEM
% problem once, and evaluates both SIF extractors on the identical field.
%
% Run from anywhere after checking out branch sif-asymmetric-mesh-audit:
%
%   Results = main_step1_same_field_compare();
%
% The returned values are preliminary method-to-method comparisons.  EDI is
% not yet treated as verified reference truth.

    here = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(here));
    addpath(genpath(repoRoot));

    fprintf('\n============================================================\n');
    fprintf('SIF AUDIT STEP 1: SAME-FIELD OLD-vs-EDI CONTROL\n');
    fprintf('============================================================\n');

    C = cfg_crack_path_two_leg_control(2.0, ...
        'plotGeom', false, ...
        'plotMesh', false);

    Results = run_crack_path_old_vs_edi(C, C.sigma0, ...
        'nthet', 100, ...
        'innerFactorEDI', 0.1, ...
        'Verbose', true);

    fprintf('\nSTEP 1 completed.\n');
    fprintf('Do not interpret the EDI column as reference truth yet.\n');
end
