function M = main_tip2h0_independent_trajectory(varargin)
%MAIN_TIP2H0_INDEPENDENT_TRAJECTORY
% Fully independent LEFM trajectory with the qualified 2*h0 tip core.
%
% Scientific isolation:
%   - frozen Stage-I initiation point/frame unchanged;
%   - theta_1=0 and Delta a=4 mm unchanged;
%   - reference exterior mesh law unchanged;
%   - EDI radii, material, loading, MTS rule, solver tolerances and gates
%     unchanged;
%   - only the structured near-tip scale is changed:
%         CoreScale = 2
%         hTip      = 0.0540246508 mm.
%
% Independence:
%   - P1 is rebuilt and physically solved with the 2*h0 core;
%   - its own KI(P1),KII(P1) determine theta_2;
%   - every later vertex is generated recursively from the 2*h0 history;
%   - no accepted reference vertex after P1 is imposed.
%
% DEFAULT IS SAFE: AllowPhysicalSolves=false.
%
% Qualification-only preflight:
%   M0 = main_tip2h0_independent_trajectory();
%
% Full independent P1--P23 run:
%   M = main_tip2h0_independent_trajectory( ...
%       'AllowPhysicalSolves',true,'MaxSegments',23);

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
    assert(exist(frozenFile,'file')==2,'tip2h0:MissingFrozen', ...
        'Recovered accepted Stage-I source is missing: %s',frozenFile);
    d=load(frozenFile,'R0');
    assert(isfield(d,'R0')&&isstruct(d.R0)&&d.R0.stage1Pass, ...
        'tip2h0:BadFrozen','Expected a passed R0 in the accepted Stage-I source.');
    R0=d.R0;

    outDir=char(opt.OutputDir);
    if isempty(outDir)
        outDir=fullfile(root,'verification','crack_path','tip_2h0_independent_run');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    p1Dir=fullfile(outDir,'p1_seed');
    pathDir=fullfile(outDir,'trajectory');
    if exist(p1Dir,'dir')~=7,mkdir(p1Dir);end
    if exist(pathDir,'dir')~=7,mkdir(pathDir);end

    coreScale=2;
    farCap=0.625;
    transition=1.0;
    cal=struct();

    fprintf('\n============================================================\n');
    fprintf('2*h0 FULLY INDEPENDENT LEFM TRAJECTORY\n');
    fprintf('============================================================\n');
    fprintf('  frozen Stage I       = %s\n',frozenFile);
    fprintf('  max physical segment = P%d\n',opt.MaxSegments);
    fprintf('  physical solves      = %d\n',logical(opt.AllowPhysicalSolves));
    fprintf('  core scale           = %.3f * h0\n',coreScale);
    fprintf('  hTip                 = %.10f mm\n',0.0540246508);
    fprintf('  exterior mesh law    = REFERENCE, UNCHANGED\n');
    fprintf('  EDI radii            = REFERENCE, UNCHANGED\n');
    fprintf('  P1 seed              = RECOMPUTED WITH 2*h0 CORE\n');
    fprintf('  P2 onward            = GENERATED ONLY FROM 2*h0 HISTORY\n');
    fprintf('  reference vertices   = NOT imposed after P1\n');

    % --------------------------------------------------------------
    % A. Rebuild/qualify the prescribed first segment using 2*h0.
    % --------------------------------------------------------------
    p1CandidateFile=fullfile(p1Dir,'tip2h0_P1_candidate.mat');
    p1QualFile=fullfile(p1Dir,'tip2h0_P1_qualification_small.mat');

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

    assert(F1.pass,'tip2h0:P1Qualification', ...
        '2*h0 P1 candidate did not pass. No trajectory solve is allowed.');
    assert(abs(F1.candidate.coreMeshControls.scale-coreScale)<=1e-14, ...
        'tip2h0:P1CoreScale','P1 candidate does not carry CoreScale=2.');
    assert(logical(F1.candidate.exteriorMeshControls.isReferenceProductionExterior), ...
        'tip2h0:P1Exterior','P1 exterior is not the reference production law.');

    M=struct();
    M.meshFamily='tip_2h0';
    M.controls=struct( ...
        'coreScale',coreScale, ...
        'hTip_mm',1e3*F1.candidate.coreMeshControls.hTip_m, ...
        'farCapOverIncrement',farCap, ...
        'transitionOverIncrement',transition, ...
        'referenceExterior',true);
    M.frozenStateFile=frozenFile;
    M.outputDir=outDir;
    M.P1Qualification=F1.summary;
    M.P1Physical=[];
    M.Path=[];

    if ~opt.AllowPhysicalSolves
        save(fullfile(outDir,'tip2h0_independent_preflight.mat'),'M','-v7');
        fprintf('\n2*h0 INDEPENDENT TRAJECTORY PREFLIGHT PASS.\n');
        fprintf('  P1 2*h0 mesh is qualified with the reference exterior law.\n');
        fprintf('  No physical solve was performed.\n');
        fprintf('  Inspect P1 structural/Williams gates before authorizing the path.\n');
        return
    end

    % --------------------------------------------------------------
    % B. Solve P1 on 2*h0. Its own field supplies theta_2.
    % --------------------------------------------------------------
    p1Checkpoint=fullfile(p1Dir,'tip2h0_P1_physical_solved.mat');
    p1Small=fullfile(p1Dir,'tip2h0_P1_physical_small.mat');

    R1=main_stage2_theta0_physical_solve( ...
        'FrozenState',R0, ...
        'Candidate',F1.candidate, ...
        'AlternativeQualifiedCandidate',true, ...
        'AllowSolve',true, ...
        'CheckpointFile',p1Checkpoint, ...
        'SaveFile',p1Small);

    assert(R1.pass,'tip2h0:P1Physical', ...
        '2*h0 P1 physical solve did not pass.');

    KI1=R1.EDI.KI_unit;
    KII1=R1.EDI.KII_unit;
    [~,theta2Deg]=kink_angle_LEFM_MTS(KI1,KII1);

    fprintf('\n2*h0 P1 SEED\n');
    fprintf('  KI(P1)     = %.15g MPa*sqrt(m)\n',KI1);
    fprintf('  KII(P1)    = %+.15g MPa*sqrt(m)\n',KII1);
    fprintf('  KII/KI     = %+.15g\n',KII1/KI1);
    fprintf('  theta_2    = %+.15g deg\n',theta2Deg);

    % --------------------------------------------------------------
    % C. Propagate with only 2*h0 physical results.
    % Exact reference-value regression gates are disabled because this is
    % a different qualified local mesh family. All ordinary qualification,
    % synthetic, solver, residual, EDI and MTS gates remain active.
    % --------------------------------------------------------------
    Path=run_incremental_crack_path( ...
        'FrozenState',R0, ...
        'MaxSegments',opt.MaxSegments, ...
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
        'MeshFamilyLabel','tip_2h0');

    M.P1Physical=R1.summary;
    M.P1EDI=R1.EDI;
    M.P1Theta2Deg=theta2Deg;
    M.Path=Path;

    % --------------------------------------------------------------
    % D. Reference comparison is post hoc only.
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
            Talt=Path.stepTable;

            common=intersect(Tref.segment,Talt.segment,'stable');
            rows=nan(numel(common),14);
            for jj=1:numel(common)
                k=common(jj);
                a=Tref(Tref.segment==k,:);
                b=Talt(Talt.segment==k,:);
                dx=1e3*(b.tip_x_m-a.tip_x_m);
                dy=1e3*(b.tip_y_m-a.tip_y_m);
                rows(jj,:)=[ ...
                    k,4*k,dx,dy,hypot(dx,dy), ...
                    b.theta_deg-a.theta_deg, ...
                    b.KI_unit-a.KI_unit, ...
                    b.KII_unit-a.KII_unit, ...
                    b.KII_over_KI-a.KII_over_KI, ...
                    b.delta_theta_next_deg-a.delta_theta_next_deg, ...
                    b.KI_unit,b.KII_unit,b.KII_over_KI,b.delta_theta_next_deg];
            end
            M.comparison=array2table(rows,'VariableNames',{ ...
                'segment','crack_length_mm', ...
                'dx_tip_mm_2h0_minus_h0','dy_tip_mm_2h0_minus_h0', ...
                'tip_separation_mm', ...
                'dtheta_deg_2h0_minus_h0', ...
                'dKI_2h0_minus_h0','dKII_2h0_minus_h0', ...
                'dq_2h0_minus_h0','dturn_deg_2h0_minus_h0', ...
                'KI_2h0','KII_2h0','q_2h0','turn_2h0_deg'});
            writetable(M.comparison, ...
                fullfile(outDir,'tip2h0_vs_reference.csv'));
        end
    end

    save(fullfile(outDir,'tip2h0_independent_run_summary.mat'),'M','-v7');

    fprintf('\n============================================================\n');
    fprintf('2*h0 INDEPENDENT TRAJECTORY RUN COMPLETE THROUGH P%d\n', ...
        Path.stepTable.segment(end));
    fprintf('============================================================\n');
    fprintf('  P1 was recomputed with the 2*h0 core.\n');
    fprintf('  The reference exterior law was retained.\n');
    fprintf('  No reference path vertex after P1 was imposed.\n');
    fprintf('  Reference comparison was generated only after propagation.\n');
    fprintf('  output = %s\n',outDir);
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
