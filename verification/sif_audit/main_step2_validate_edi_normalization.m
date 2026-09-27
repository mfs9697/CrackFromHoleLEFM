function Out = main_step2_validate_edi_normalization()
%MAIN_STEP2_VALIDATE_EDI_NORMALIZATION
% Second executable gate of the SIF asymmetric-mesh audit.
%
% Uses exact leading-order Williams displacement fields on an independently
% generated polar crack annulus. No physical FEM crack solve is involved.
%
% The production EDI implementation is not modified by this test.

    here = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(here));
    addpath(genpath(repoRoot));

    fprintf('\n============================================================\n');
    fprintf('SIF AUDIT STEP 2: EDI NORMALIZATION / SIGN\n');
    fprintf('============================================================\n');

    Out = validate_EDI_Williams_fields( ...
        'NrList', [8 16], ...
        'NthList', [64 128], ...
        'Verbose', true);

    fprintf('\nSTEP 2 regression completed.\n');
    fprintf('Fine-mesh Williams-field recovery must remain within the configured tolerance.\n');
    if isfield(Out,'checks') && isfield(Out.checks,'pass') && Out.checks.pass
        fprintf('STATUS: PASS (tolerance = %.3g)\n', Out.fineTolerance);
    end
end
