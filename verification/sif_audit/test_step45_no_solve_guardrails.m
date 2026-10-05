function result=test_step45_no_solve_guardrails()
%TEST_STEP45_NO_SOLVE_GUARDRAILS
% Verifies the default Step45 path NEVER silently creates a FEM result.
% This test itself needs NO meshing, finite-element solve or EDI.
tfile=[tempname '_nonexistent_step45.mat'];
assert(exist(tfile,'file')~=2);
blockedPrepare=false;
blockedPost=false;
try
    main_step45_prepare_symmetric_checkpoint([], ...
        'CheckpointFile',tfile);
catch ME
    blockedPrepare=strcmp(ME.identifier,'step45:NeedsPriorField');
    if ~blockedPrepare,rethrow(ME);end
end
try
    main_step45_symmetric_field_leakage(tfile);
catch ME
    blockedPost=strcmp(ME.identifier,'step45:CheckpointAbsent');
    if ~blockedPost,rethrow(ME);end
end
result=struct('passed',blockedPrepare&&blockedPost, ...
    'prepare_requires_reuse_or_explicit_solve',blockedPrepare, ...
    'postprocessor_rejects_missing_checkpoint',blockedPost);
disp(result);
assert(result.passed,'step45:GuardrailRegression', ...
    'Step45 no-solve safety guards failed.');
end
