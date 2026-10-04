function R63r=main_step63r_recover_lost_physical_field(varargin)
%MAIN_STEP63R_RECOVER_LOST_PHYSICAL_FIELD
% Guarded replacement of the LOST Step63 physical checkpoint.
%
% This is NOT a new mesh experiment. It reproduces the exact selected C03
% mesh and repeats the already-authorized physical problem only because the
% first solved checkpoint was lost locally after it had completed.
%
% DEFAULT IS SAFE:
%   AllowRecoverySolve = false
%
% With explicit investigator authorization:
%   R63r=main_step63r_recover_lost_physical_field( ...
%       'AllowRecoverySolve',true);
%
% The recovery:
%   - uses main_step63_calibrated_asymmetric_physical_solve;
%   - inherits its deterministic C03 reconstruction and all provenance gates;
%   - performs at most one replacement physical solve;
%   - reruns native COD;
%   - compares the recovered COD fingerprint against the recorded first-run
%     Step63 console values before declaring the field recovered.
%
% NO EDI is run here. Step64 remains a separate postprocessing command.
%
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');

ip=inputParser;
addParameter(ip,'AllowRecoverySolve',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'CheckpointFile', ...
    fullfile(vdir,'step63_calibrated_asymmetric_physical_solved.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'ReportFile', ...
    fullfile(vdir,'step63r_recovered_physical_field_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;

addpath(genpath(root));
assert_recovery_branch(root);

cp=char(opt.CheckpointFile);
reportFile=char(opt.ReportFile);

fprintf('\n============================================================\n');
fprintf('STEP 63R: RECOVER LOST STEP63 PHYSICAL FIELD\n');
fprintf('============================================================\n');
fprintf('  This is a replacement of the lost first Step63 checkpoint.\n');
fprintf('  Same C03 mesh, same material/loading/anchors, same solver.\n');
fprintf('  NO EDI in this recovery step.\n');

if exist(cp,'file')==2
    fprintf('  A Step63 checkpoint already exists. No replacement solve will run.\n');
else
    if ~opt.AllowRecoverySolve
        error('step63r:ExplicitRecoveryAuthorizationRequired', ...
            ['The original Step63 field is confirmed lost. A replacement ', ...
             'physical solve requires explicit authorization. Rerun with ', ...
             '''AllowRecoverySolve'',true only if authorized.']);
    end
    fprintf('  No Step63 checkpoint exists.\n');
    fprintf('  Starting ONE explicitly authorized replacement solve.\n');
end

[P63,O63]=main_step63_calibrated_asymmetric_physical_solve( ...
    'AllowSolve',opt.AllowRecoverySolve, ...
    'CheckpointFile',cp);

% -------------------------------------------------------------------------
% Compare recovered COD with the first successful Step63 run recorded in
% the investigator console log. These reference numbers are intentionally
% only the displayed values, so tolerances are looser than machine precision.
% -------------------------------------------------------------------------
F=O63.fitTable;
R=O63.rawTable;
if height(F)~=8 || height(R)~=5
    error('step63r:CODShape','Recovered COD tables have unexpected size.');
end

expectedWindows=[ ...
    .04 .20 1 38;
    .04 .20 2 38;
    .04 .30 1 55;
    .04 .30 2 55;
    .08 .30 1 44;
    .08 .30 2 44;
    .12 .30 1 34;
    .12 .30 2 34];
actualWindows=[F.lower_r_over_a0,F.upper_r_over_a0,F.degree,F.n_native];
if max(abs(actualWindows(:)-expectedWindows(:)))>1e-12
    error('step63r:CODWindows','Recovered COD windows/sampling changed.');
end

expectedFitRatio=[ ...
    1.0526e-4;
    1.0562e-4;
    1.0489e-4;
    1.0571e-4;
    1.0456e-4;
    1.0582e-4;
    1.0416e-4;
    1.0583e-4];
expectedRawMedian=[ ...
    1.0365e-4;
    1.0178e-4;
    9.9274e-5;
    9.5709e-5;
    9.0841e-5];

fitDelta=F.ratio_COD-expectedFitRatio;
rawDelta=R.median_raw_ratio-expectedRawMedian;

% Printed first-run values were rounded. 5e-8 absolute is ~0.05% of the
% target ratio and comfortably larger than display-rounding uncertainty.
tolRatio=5e-8;
fitFingerprintPass=max(abs(fitDelta))<=tolRatio;
rawFingerprintPass=max(abs(rawDelta))<=tolRatio;
positiveSignalPass=all(F.ratio_COD>0) && all(R.median_raw_ratio>0);
rangePass=min(F.ratio_COD)>1.03e-4 && max(F.ratio_COD)<1.07e-4;

fingerprintTable=table( ...
    expectedFitRatio,F.ratio_COD,fitDelta, ...
    'VariableNames',{'firstRunDisplayedRatio','recoveredRatio','difference'});

rawFingerprintTable=table( ...
    expectedRawMedian,R.median_raw_ratio,rawDelta, ...
    'VariableNames',{'firstRunDisplayedMedian','recoveredMedian','difference'});

gates=struct( ...
    'fitFingerprintPass',fitFingerprintPass, ...
    'rawFingerprintPass',rawFingerprintPass, ...
    'positiveSignalPass',positiveSignalPass, ...
    'rangePass',rangePass);
gates.recoveryPass=all(structfun(@logical,gates));

fprintf('\nRECOVERY COD FINGERPRINT — FIRST RUN vs REPLACEMENT\n');
disp(fingerprintTable);
fprintf('\nRECOVERY RAW-BAND FINGERPRINT\n');
disp(rawFingerprintTable);
disp(gates);

if ~gates.recoveryPass
    error('step63r:RecoveryFingerprintMismatch', ...
        ['Replacement field does not reproduce the recorded first-run ', ...
         'Step63 COD fingerprint within the predeclared tolerance. ', ...
         'Do not run Step64.']);
end

cpHash=sha256(cp);
R63r=struct( ...
    'checkpointPath',cp, ...
    'checkpointSHA256',cpHash, ...
    'replacementSolvePerformed',logical(P63.newSolve), ...
    'firstRunCheckpointWasLost',true, ...
    'fitFingerprintTable',fingerprintTable, ...
    'rawFingerprintTable',rawFingerprintTable, ...
    'gates',gates, ...
    'tolRatio',tolRatio, ...
    'O63',O63, ...
    'noEDI',true, ...
    'interpretation',['Exact C03 physical problem repeated solely to replace ', ...
      'a lost solved checkpoint; recovered COD agrees with the recorded ', ...
      'first successful Step63 fingerprint.']);

save(reportFile,'R63r','-v7');
fprintf('  Recovered checkpoint SHA-256: %s\n',cpHash);
fprintf('  Recovery report saved: %s\n',reportFile);
fprintf('STEP63R PASS. Keep this branch checked out; Step64 may now use the field.\n');
end

function h=sha256(path)
f=fopen(path,'rb');
if f<0,error('step63r:HashOpen','Cannot open %s.',path);end
guard=onCleanup(@()fclose(f));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(f)
    bytes=fread(f,1024*1024,'*uint8');
    md.update(typecast(bytes,'int8'));
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
clear guard
end

function assert_recovery_branch(root)
[status,b]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0&&strcmp(strtrim(b),'audit/step63r-recover-lost-physical-field'), ...
    'step63r:Branch', ...
    'Run Step63R only on audit/step63r-recover-lost-physical-field.');
end
