function M = main_m1_independent_trajectory(varargin)
%MAIN_M1_INDEPENDENT_TRAJECTORY
% Fully independent LEFM trajectory on the qualified M1 exterior mesh family.
%
% Scientific rule:
%   - Stage-I initiation point/frame are the frozen accepted values.
%   - theta_1=0 and Delta a=4 mm remain prescribed exactly as before.
%   - P1 is rebuilt and solved on M1 to obtain an M1-specific KI,KII seed.
%   - theta_2 is computed from that M1 P1 field.
%   - every P2...Pk geometry is then generated recursively from M1 SIFs/MTS.
%   - no reference trajectory vertex after P1 is imposed.
%
% M1 exterior controls:
%   far cap               = 1.25 Delta a = 5 mm
%   transition            = 1.00 Delta a = 4 mm
%   far slope             = 0.15
%   boundary metric growth= 0.35
%
% The reflection-paired tip core, EDI radii, COD windows, material, loading,
% MTS selector, solver tolerances, and acceptance gates remain unchanged.
%
% DEFAULT IS SAFE: AllowPhysicalSolves=false.
%
% Qualification-only preflight:
%   M = main_m1_independent_trajectory();
%
% Full independent P1--P23 run:
%   M = main_m1_independent_trajectory('AllowPhysicalSolves',true);

    ip=inputParser;
    addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'MaxSegments',23,@(x)isnumeric(x)&&isscalar(x)&& ...
        isfinite(x)&&x>=2&&x==round(x));
    addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FastEDI',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ReuseCandidates',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PlotEachStep',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(root));

    frozenFile=fullfile(root,'paper','data','accepted_stage1_source.mat');
    assert(exist(frozenFile,'file')==2,'m1path:MissingFrozen', ...
        'Recovered accepted Stage-I source is missing: %s',frozenFile);
    d=load(frozenFile,'R0');
    assert(isfield(d,'R0')&&isstruct(d.R0)&&d.R0.stage1Pass, ...
        'm1path:BadFrozen','Expected a passed R0 in the recovered Stage-I source.');
    R0=d.R0;

    outDir=char(opt.OutputDir);
    if isempty(outDir)
        outDir=fullfile(root,'verification','crack_path','m1_independent_run');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    p1Dir=fullfile(outDir,'p1_seed');
    pathDir=fullfile(outDir,'trajectory');
    if exist(p1Dir,'dir')~=7,mkdir(p1Dir);end
    if exist(pathDir,'dir')~=7,mkdir(pathDir);end

    cal=struct('farSlope',0.15,'boundaryMetricGrowth',0.35);
    farCap=1.25;
    transition=1.0;

    fprintf('\n============================================================\n');
    fprintf('M1 FULLY INDEPENDENT LEFM TRAJECTORY\n');
    fprintf('============================================================\n');
    fprintf('  frozen Stage I       = %s\n',frozenFile);
    fprintf('  max physical segment = P%d\n',opt.MaxSegments);
    fprintf('  physical solves      = %d\n',logical(opt.AllowPhysicalSolves));
    fprintf('  M1 far cap           = %.3f * DeltaA\n',farCap);
    fprintf('  M1 transition        = %.3f * DeltaA\n',transition);
    fprintf('  M1 far slope         = %.3f\n',cal.farSlope);
    fprintf('  M1 boundary growth   = %.3f\n',cal.boundaryMetricGrowth);
    fprintf('  P1 seed              = RECOMPUTED ON M1\n');
    fprintf('  P2 onward            = GENERATED ONLY FROM M1 HISTORY\n');
    fprintf('  reference vertices   = NOT imposed after P1\n');

    % --------------------------------------------------------------
    % Step A. Build the prescribed first 4-mm segment on M1.
    % --------------------------------------------------------------
    p1CandidateFile=fullfile(p1Dir,'M1_P1_candidate.mat');
    p1QualFile=fullfile(p1Dir,'M1_P1_qualification_small.mat');

    F1=main_stage2_embed_scaled_core_full_domain_theta0( ...
        'FrozenState',R0, ...
        'ExteriorVerbose',true, ...
        'ExteriorFarCapOverA0',farCap, ...
        'ExteriorTransitionOverA0',transition, ...
        'ExteriorCalibration',cal, ...
        'SaveCandidate',true, ...
        'CandidateFile',p1CandidateFile, ...
        'SaveCompact',true, ...
        'CompactFile',p1QualFile, ...
        'Plot',false);

    assert(F1.pass,'m1path:P1Qualification', ...
        'M1 P1 candidate did not pass. No trajectory solve is allowed.');

    M=struct();
    M.meshFamily='M1';
    M.controls=struct('farCapOverIncrement',farCap, ...
        'transitionOverIncrement',transition, ...
        'farSlope',cal.farSlope, ...
        'boundaryMetricGrowth',cal.boundaryMetricGrowth);
    M.frozenStateFile=frozenFile;
    M.outputDir=outDir;
    M.P1Qualification=F1.summary;
    M.P1Physical=[];
    M.Path=[];

    if ~opt.AllowPhysicalSolves
        save(fullfile(outDir,'M1_independent_preflight.mat'),'M','-v7');
        fprintf('\nM1 INDEPENDENT TRAJECTORY PREFLIGHT PASS.\n');
        fprintf('  P1 M1 mesh is qualified. No physical solve was performed.\n');
        fprintf('  Rerun with ''AllowPhysicalSolves'',true for the independent path.\n');
        return
    end

    % --------------------------------------------------------------
    % Step B. Solve P1 on M1. This removes the final reference-mesh seed.
    % --------------------------------------------------------------
    p1Checkpoint=fullfile(p1Dir,'M1_P1_physical_solved.mat');
    p1Small=fullfile(p1Dir,'M1_P1_physical_small.mat');

    R1=main_stage2_theta0_physical_solve( ...
        'FrozenState',R0, ...
        'Candidate',F1.candidate, ...
        'AlternativeQualifiedCandidate',true, ...
        'AllowSolve',true, ...
        'CheckpointFile',p1Checkpoint, ...
        'SaveFile',p1Small);

    assert(R1.pass,'m1path:P1Physical','M1 P1 physical solve did not pass.');

    KI1=R1.EDI.KI_unit;
    KII1=R1.EDI.KII_unit;
    [delta2,theta2Deg]=kink_angle_LEFM_MTS(KI1,KII1); %#ok<ASGLU>

    fprintf('\nM1 P1 SEED\n');
    fprintf('  KI(P1)     = %.15g MPa*sqrt(m)\n',KI1);
    fprintf('  KII(P1)    = %+.15g MPa*sqrt(m)\n',KII1);
    fprintf('  KII/KI     = %+.15g\n',KII1/KI1);
    fprintf('  theta_2    = %+.15g deg\n',theta2Deg);

    % --------------------------------------------------------------
    % Step C. Propagate P2...Pmax using only M1 physical results.
    % Reference exact-value regression gates are intentionally disabled:
    % this is a different qualified mesh family, not a reproduction run.
    % All structural, synthetic, solver, residual, EDI, and MTS gates remain.
    % --------------------------------------------------------------
    Path=run_incremental_crack_path( ...
        'FrozenState',R0, ...
        'MaxSegments',opt.MaxSegments, ...
        'AllowPhysicalSolves',true, ...
        'RunSynthetic',opt.RunSynthetic, ...
        'FastEDI',opt.FastEDI, ...
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
        'MeshFamilyLabel','M1');

    M.P1Physical=R1.summary;
    M.P1EDI=R1.EDI;
    M.P1Theta2Deg=theta2Deg;
    M.Path=Path;

    % --------------------------------------------------------------
    % Step D. Compare independent M1 history with the accepted reference
    % only after the M1 trajectory has been generated.
    % --------------------------------------------------------------
    refStateFile=fullfile(root,'verification','crack_path', ...
        'final_clean_run','path_run_state.mat');
    M.referenceStateFile=refStateFile;
    M.comparison=table();

    if exist(refStateFile,'file')==2
        q=load(refStateFile,'State');
        if isfield(q,'State')&&isstruct(q.State)
            Rref=q.State;
            names=cellstr(string(Rref.rowVariableNames));
            Tref=array2table(Rref.rowsThroughCompleted,'VariableNames',names);
            Tm1=Path.stepTable;

            common=intersect(Tref.segment,Tm1.segment,'stable');
            rows=nan(numel(common),13);
            for jj=1:numel(common)
                k=common(jj);
                a=Tref(Tref.segment==k,:);
                b=Tm1(Tm1.segment==k,:);
                rows(jj,:)=[ ...
                    k,4*k, ...
                    1e3*(b.tip_x_m-a.tip_x_m), ...
                    1e3*(b.tip_y_m-a.tip_y_m), ...
                    b.theta_deg-a.theta_deg, ...
                    b.KI_unit-a.KI_unit, ...
                    b.KII_unit-a.KII_unit, ...
                    b.KII_over_KI-a.KII_over_KI, ...
                    b.delta_theta_next_deg-a.delta_theta_next_deg, ...
                    b.KI_unit, ...
                    b.KII_unit, ...
                    b.KII_over_KI, ...
                    b.delta_theta_next_deg];
            end
            M.comparison=array2table(rows,'VariableNames',{ ...
                'segment','crack_length_mm', ...
                'dx_tip_mm_M1_minus_ref','dy_tip_mm_M1_minus_ref', ...
                'dtheta_deg_M1_minus_ref', ...
                'dKI_M1_minus_ref','dKII_M1_minus_ref', ...
                'dq_M1_minus_ref','dturn_deg_M1_minus_ref', ...
                'KI_M1','KII_M1','q_M1','turn_M1_deg'});
            writetable(M.comparison,fullfile(outDir,'M1_vs_reference.csv'));
        end
    end

    save(fullfile(outDir,'M1_independent_run_summary.mat'),'M','-v7');

    fprintf('\n============================================================\n');
    fprintf('M1 INDEPENDENT TRAJECTORY RUN COMPLETE THROUGH P%d\n', ...
        Path.stepTable.segment(end));
    fprintf('============================================================\n');
    fprintf('  M1 P1 was recomputed on M1.\n');
    fprintf('  No reference path vertex after P1 was imposed.\n');
    fprintf('  Reference comparison was generated only after propagation.\n');
    fprintf('  output = %s\n',outDir);
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
