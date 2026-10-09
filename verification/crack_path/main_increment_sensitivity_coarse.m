function M = main_increment_sensitivity_coarse(varargin)
%MAIN_INCREMENT_SENSITIVITY_COARSE
% Independent crack-path run for a controlled crack-increment study.
%
% The experiment compares Delta a = 4 mm, 2 mm, and 1 mm using the SAME
% dimensionless coarse numerical family:
%
%   structured tip core : CoreScale = 2
%   hTip/Delta a        : 2*0.00675308135
%   r_i/Delta a         : 0.10
%   r_o/Delta a         : 0.65
%   r_c/Delta a         : 0.75
%
%   exterior M1 family:
%     far cap/Delta a   : 1.25
%     transition/Delta a: 1.00
%     far slope         : 0.15
%     boundary growth   : 0.35
%
% Hence halving Delta a halves the local/core/EDI/exterior length scales.
% This is a convergence test of the fixed-increment numerical scheme as a
% dimensionless family; it is NOT an experiment holding absolute mesh sizes
% fixed while changing only Delta a.
%
% DEFAULT IS SAFE: qualification only, no physical solve.
%
% Examples
% --------
% 4-mm preflight:
%   M4q = main_increment_sensitivity_coarse('IncrementMM',4);
%
% 2-mm preflight:
%   M2q = main_increment_sensitivity_coarse('IncrementMM',2);
%
% 1-mm preflight:
%   M1q = main_increment_sensitivity_coarse('IncrementMM',1);
%
% Full 4-mm trajectory through 92 mm:
%   M4 = main_increment_sensitivity_coarse( ...
%       'IncrementMM',4,'TargetLengthMM',92,'AllowPhysicalSolves',true);
%
% Full 2-mm trajectory through 92 mm:
%   M2 = main_increment_sensitivity_coarse( ...
%       'IncrementMM',2,'TargetLengthMM',92,'AllowPhysicalSolves',true);
%
% Full 1-mm trajectory through 92 mm:
%   M1 = main_increment_sensitivity_coarse( ...
%       'IncrementMM',1,'TargetLengthMM',92,'AllowPhysicalSolves',true);

    ip=inputParser;
    addParameter(ip,'IncrementMM',4,@(x)isnumeric(x)&&isscalar(x)&& ...
        isfinite(x)&&any(abs(x-[1 2 4])<=1e-12));
    addParameter(ip,'TargetLengthMM',92,@(x)isnumeric(x)&&isscalar(x)&& ...
        isfinite(x)&&x>0);
    addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FastEDI',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ReuseCandidates',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'Resume',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PlotEachStep',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(root));

    da=1e-3*opt.IncrementMM;
    target=1e-3*opt.TargetLengthMM;
    nSeg=target/da;
    assert(abs(nSeg-round(nSeg))<=1e-12,'incstudy:TargetNotDivisible', ...
        'TargetLengthMM must be an integer multiple of IncrementMM.');
    nSeg=round(nSeg);
    assert(nSeg>=2,'incstudy:TargetTooShort', ...
        'TargetLengthMM must contain at least two increments for the trajectory driver.');

    frozenFile=fullfile(root,'paper','data','accepted_stage1_source.mat');
    assert(exist(frozenFile,'file')==2,'incstudy:MissingFrozen', ...
        'Accepted Stage-I source is missing: %s',frozenFile);
    d=load(frozenFile,'R0');
    assert(isfield(d,'R0')&&isstruct(d.R0)&&d.R0.stage1Pass, ...
        'incstudy:BadFrozen','Expected a passed accepted Stage-I R0.');
    R0=d.R0;

    % Stage I selects the initiation point/frame.  Its reserved 4-mm first
    % segment is a Stage-II numerical choice, not part of the hole-only
    % physical solve.  Override only that reserved increment in memory.
    assert(isfield(R0,'summary')&&istable(R0.summary)&&height(R0.summary)==1, ...
        'incstudy:FrozenSummary','Frozen R0.summary is required.');
    assert(ismember('a0_reserved_m',R0.summary.Properties.VariableNames), ...
        'incstudy:IncrementField','Frozen summary lacks a0_reserved_m.');
    originalReserved=R0.summary.a0_reserved_m;
    R0.summary.a0_reserved_m=da;

    % Coarse family used for BOTH increments.
    coreScale=2;
    farCap=1.25;
    transition=1.00;
    cal=struct('farSlope',0.15,'boundaryMetricGrowth',0.35);

    outDir=char(opt.OutputDir);
    if isempty(outDir)
        tag=strrep(sprintf('da_%gmm',opt.IncrementMM),'.','p');
        outDir=fullfile(root,'verification','crack_path', ...
            'increment_sensitivity_coarse',tag);
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end
    p1Dir=fullfile(outDir,'p1_seed');
    pathDir=fullfile(outDir,'trajectory');
    if exist(p1Dir,'dir')~=7,mkdir(p1Dir);end
    if exist(pathDir,'dir')~=7,mkdir(pathDir);end

    fprintf('\n============================================================\n');
    fprintf('CRACK-INCREMENT SENSITIVITY: COARSE DIMENSIONLESS FAMILY\n');
    fprintf('============================================================\n');
    fprintf('  Delta a              = %.6f mm\n',opt.IncrementMM);
    fprintf('  target crack length  = %.6f mm\n',opt.TargetLengthMM);
    fprintf('  target segments      = %d\n',nSeg);
    fprintf('  physical solves      = %d\n',logical(opt.AllowPhysicalSolves));
    fprintf('  Stage-I point/frame  = accepted, unchanged\n');
    fprintf('  reserved increment   = %.6f -> %.6f mm (in-memory only)\n', ...
        1e3*originalReserved,opt.IncrementMM);
    fprintf('  core scale           = 2 (coarse)\n');
    fprintf('  hTip                 = %.10f mm\n', ...
        1e3*coreScale*0.00675308135*da);
    fprintf('  EDI r_i / r_o        = %.6f / %.6f mm\n', ...
        1e3*.10*da,1e3*.65*da);
    fprintf('  core radius          = %.6f mm\n',1e3*.75*da);
    fprintf('  M1 far cap           = %.6f mm = 1.25 Delta a\n',1e3*farCap*da);
    fprintf('  M1 transition        = %.6f mm = 1.00 Delta a\n',1e3*transition*da);
    fprintf('  M1 far slope         = %.3f\n',cal.farSlope);
    fprintf('  M1 boundary growth   = %.3f\n',cal.boundaryMetricGrowth);

    % --------------------------------------------------------------
    % A. P1 qualification on the selected increment and coarse family.
    % --------------------------------------------------------------
    p1CandidateFile=fullfile(p1Dir,'P1_candidate.mat');
    p1QualFile=fullfile(p1Dir,'P1_qualification_small.mat');

    F1=main_stage2_embed_scaled_core_full_domain_theta0( ...
        'FrozenState',R0, ...
        'CoreScale',coreScale, ...
        'ExteriorVerbose',true, ...
        'ExteriorFarCapOverA0',farCap, ...
        'ExteriorTransitionOverA0',transition, ...
        'ExteriorCalibration',cal, ...
        'SaveCandidate',true, ...
        'CandidateFile',p1CandidateFile, ...
        'SaveCompact',true, ...
        'CompactFile',p1QualFile, ...
        'Plot',false);

    assert(F1.pass,'incstudy:P1Qualification', ...
        'P1 did not qualify for Delta a = %.6g mm.',opt.IncrementMM);
    assert(abs(F1.candidate.coreMeshControls.scale-coreScale)<=1e-14, ...
        'incstudy:P1CoreScale','P1 candidate does not use CoreScale=2.');
    assert(F1.candidate.coreMeshControls.expectedCoreT3==3318 && ...
           F1.candidate.coreMeshControls.expectedEDIElements==2976 && ...
           F1.candidate.coreMeshControls.expectedOptimizedSupport==2700 && ...
           isequal(F1.candidate.coreMeshControls.expectedNativeSamples(:), ...
                   [19;28;23;18]), ...
        'incstudy:P1Fingerprint', ...
        'P1 candidate does not carry the qualified CoreScale=2 fingerprints.');
    assert(isfield(F1.candidate,'exteriorMeshControls')&& ...
        ~logical(F1.candidate.exteriorMeshControls.isReferenceProductionExterior), ...
        'incstudy:P1Exterior','P1 candidate does not carry M1 exterior provenance.');

    M=struct();
    M.study='crack_increment_sensitivity_coarse_family';
    M.increment_mm=opt.IncrementMM;
    M.targetLength_mm=opt.TargetLengthMM;
    M.targetSegments=nSeg;
    M.frozenStateFile=frozenFile;
    M.originalReservedIncrement_mm=1e3*originalReserved;
    M.controls=struct( ...
        'coreScale',coreScale, ...
        'hTipOverIncrement',coreScale*0.00675308135, ...
        'rInnerOverIncrement',.10, ...
        'rOuterOverIncrement',.65, ...
        'rCoreOverIncrement',.75, ...
        'farCapOverIncrement',farCap, ...
        'transitionOverIncrement',transition, ...
        'farSlope',cal.farSlope, ...
        'boundaryMetricGrowth',cal.boundaryMetricGrowth);
    M.outputDir=outDir;
    M.P1Qualification=F1.summary;
    M.P1Physical=[];
    M.Path=[];

    if ~opt.AllowPhysicalSolves
        save(fullfile(outDir,'preflight.mat'),'M','-v7');
        fprintf('\nCOARSE INCREMENT-STUDY PREFLIGHT PASS.\n');
        fprintf('  No physical solve was performed.\n');
        return
    end

    % --------------------------------------------------------------
    % B. P1 physical seed. Reuse an accepted compact result if present.
    % --------------------------------------------------------------
    p1Checkpoint=fullfile(p1Dir,'P1_physical_solved.mat');
    p1Small=fullfile(p1Dir,'P1_physical_small.mat');

    reuseP1=false;
    if opt.Resume && exist(p1Small,'file')==2
        z=load(p1Small,'R');
        assert(isfield(z,'R')&&isstruct(z.R),'incstudy:P1ReuseMismatch', ...
            'Cached P1 file does not contain an accepted R.');
        validate_increment_study_p1_cache(z.R,F1.candidate,R0,p1Checkpoint);
        R1=z.R;
        reuseP1=true;
        fprintf('\nReusing accepted compatible P1 physical seed: %s\n',p1Small);
    end

    if ~reuseP1
        R1=main_stage2_theta0_physical_solve( ...
            'FrozenState',R0, ...
            'Candidate',F1.candidate, ...
            'AlternativeQualifiedCandidate',true, ...
            'AllowSolve',true, ...
            'CheckpointFile',p1Checkpoint, ...
            'SaveFile',p1Small);
    end
    assert(R1.pass,'incstudy:P1Physical','P1 physical seed did not pass.');

    KI1=R1.EDI.KI_unit;
    KII1=R1.EDI.KII_unit;
    [~,theta2Deg]=kink_angle_LEFM_MTS(KI1,KII1);

    fprintf('\nCOARSE-FAMILY P1 SEED\n');
    fprintf('  KI(P1)     = %.15g MPa*sqrt(m)\n',KI1);
    fprintf('  KII(P1)    = %+.15g MPa*sqrt(m)\n',KII1);
    fprintf('  KII/KI     = %+.15g\n',KII1/KI1);
    fprintf('  theta_2    = %+.15g deg\n',theta2Deg);

    % --------------------------------------------------------------
    % C. Independent trajectory.
    % --------------------------------------------------------------
    runArgs={ ...
        'FrozenState',R0, ...
        'MaxSegments',nSeg, ...
        'AllowPhysicalSolves',true, ...
        'RunSynthetic',opt.RunSynthetic, ...
        'FastEDI',opt.FastEDI, ...
        'CoreScale',coreScale, ...
        'ReuseCandidates',opt.ReuseCandidates, ...
        'PlotEachStep',opt.PlotEachStep, ...
        'RegressionGates',false, ...
        'StopAtCoreClearance',true, ...
        'OutputDir',pathDir, ...
        'SeedKI',KI1, ...
        'SeedKII',KII1, ...
        'ExteriorFarCapOverIncrement',farCap, ...
        'ExteriorTransitionOverIncrement',transition, ...
        'ExteriorCalibration',cal, ...
        'MeshFamilyLabel',sprintf('increment_%gmm_coarse_M1_2h0',opt.IncrementMM)};

    stateFile=fullfile(pathDir,'path_run_state.mat');
    if opt.Resume && exist(stateFile,'file')==2
        runArgs=[runArgs,{'ResumeStateFile',stateFile,'ResumeSourceDir',pathDir}]; %#ok<AGROW>
        fprintf('  trajectory resume    = %s\n',stateFile);
    end

    Path=run_incremental_crack_path(runArgs{:});

    M.P1Physical=R1.summary;
    M.P1EDI=R1.EDI;
    M.P1Theta2Deg=theta2Deg;
    M.Path=Path;

    % Portable result table.
    T=Path.stepTable;
    T.crack_length_mm=1e3*da*T.segment;
    csvFile=fullfile(outDir,'states.csv');
    writetable(T,csvFile);

    V=Path.vertices;
    vertexIndex=(0:size(V,1)-1).';
    Vtab=table(vertexIndex,1e3*da*vertexIndex,V(:,1),V(:,2), ...
        'VariableNames',{'vertex','crack_length_mm','x_m','y_m'});
    verticesFile=fullfile(outDir,'vertices.csv');
    writetable(Vtab,verticesFile);

    M.statesCsv=csvFile;
    M.verticesCsv=verticesFile;
    save(fullfile(outDir,'increment_study_summary.mat'),'M','-v7');

    fprintf('\n============================================================\n');
    fprintf('COARSE INCREMENT-STUDY RUN COMPLETE\n');
    fprintf('============================================================\n');
    fprintf('  Delta a             = %.6f mm\n',opt.IncrementMM);
    fprintf('  solved through      = P%d (%.6f mm)\n', ...
        T.segment(end),T.crack_length_mm(end));
    fprintf('  stop reason         = %s\n',Path.stopReason);
    fprintf('  states CSV          = %s\n',csvFile);
    fprintf('  vertices CSV        = %s\n',verticesFile);
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
