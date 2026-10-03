% MAIN_STEP45_LOCAL
% One-command local entry point for the Step 45 symmetric FEM audit.
% Run from MATLAB after pulling branch sif-asymmetric-mesh-audit:
%   main_step45_local
%
% Uses, in order: (1) recognized saved Step-45 checkpoint,
% (2) the existing workspace O18 solved Step-18 results.
% If neither exists, STOPS WITHOUT SOLVING unless the investigator first
% explicitly types STEP45_ALLOW_NEW_SOLVE = true in the MATLAB Command Window.
% An explicitly authorized run constructs/solves precisely one centered,
% zero-angle Stage-II FEM field (Npoly=240), and immediately checkpoints it.
% This first batch performs COD ONLY, never an EDI integration.
%
% No dependence on the current working directory: this script lives at
% <repo>/verification/sif_audit/main_step45_local.m and writes compact
% results and the checkpoint under <repo>/verification.

step45ScriptDir=fileparts(mfilename('fullpath'));
step45RepoRoot=fileparts(fileparts(step45ScriptDir));
addpath(genpath(step45RepoRoot));
step45DataDir=fullfile(step45RepoRoot,'verification');
step45Checkpoint=fullfile(step45DataDir,'step45_symmetric_theta0_solved.mat');
step45OutputPrefix=fullfile(step45DataDir,'step45_symmetric_field_leakage');

fprintf('\n============================================================\n');
fprintf('STEP 45 LOCAL — PHASE 1: ACTUAL SYMMETRIC FEM, NATIVE COD\n');
fprintf('============================================================\n');
test_step45_no_solve_guardrails();

if exist('STEP45_ALLOW_NEW_SOLVE','var')~=1
    STEP45_ALLOW_NEW_SOLVE=false;
end
if ~(islogical(STEP45_ALLOW_NEW_SOLVE) && isscalar(STEP45_ALLOW_NEW_SOLVE))
    error('step45:AllowSolveFlag', ...
        'STEP45_ALLOW_NEW_SOLVE must be explicitly set to true or false.');
end

if exist(step45Checkpoint,'file')==2
    % This helper recognizes only an identified Step-45 checkpoint and
    % refuses to overwrite or silently recompute a solved field.
    fprintf('\nReusing existing Step-45 symmetric FEM checkpoint.\n');
    P45=main_step45_prepare_symmetric_checkpoint([], ...
        'CheckpointFile',step45Checkpoint);
elseif exist('O18','var')==1 && isstruct(O18) && ~isempty(O18)
    fprintf('\nReusing the existing workspace O18 theta=0 FEM field.\n');
    P45=main_step45_prepare_symmetric_checkpoint(O18, ...
        'CheckpointFile',step45Checkpoint);
elseif STEP45_ALLOW_NEW_SOLVE
    fprintf(['\nExplicit authorization received. ', ...
        'Building exactly ONE centered theta=0 Stage-II solution.\n']);
    P45=main_step45_prepare_symmetric_checkpoint([], ...
        'AllowSolve',true,'Npoly',240, ...
        'CheckpointFile',step45Checkpoint);
else
    fprintf(['\nNO SAVED STEP-18 FIELD OR STEP-45 CHECKPOINT FOUND.\n', ...
        'No FEM solve was started. If the original Step-18 MAT exists,\n', ...
        'load O18 into the MATLAB workspace and run main_step45_local.\n', ...
        'Otherwise authorize exactly ONE symmetric FEM solve by typing:\n', ...
        '    STEP45_ALLOW_NEW_SOLVE = true;\n', ...
        '    main_step45_local\n', ...
        'That call will checkpoint the solved field automatically.\n']);
    return
end

% The first step deliberately stops after inexpensive COD diagnostics.
% Review results before running matched, checkpointed 16-point EDI.
O45=main_step45_symmetric_field_leakage(P45, ...
    'SavePrefix',step45OutputPrefix);

fprintf('\n============================================================\n');
fprintf('STEP 45 PHASE 1 RESULTS TO RETURN\n');
fprintf('============================================================\n');
fprintf('\nTip-adjacent T3 topology:\n');
disp(O45.tipTopology);
fprintf('\nOpposite-face abscissa mismatch [m]:\n');
disp(O45.faceGridMismatch);
fprintf('\nPointwise COD bands:\n');
disp(O45.rawBands);
fprintf('\nLinear/quadratic COD extrapolation:\n');
disp(O45.fitTable);
fprintf('\nProposed COD numerical symmetry gate:\n');
disp(O45.CODgates);
fprintf(['\nPHASE 1 COMPLETE. Send this output back for interpretation.\n', ...
    'No EDI integration is started by main_step45_local.\n']);
