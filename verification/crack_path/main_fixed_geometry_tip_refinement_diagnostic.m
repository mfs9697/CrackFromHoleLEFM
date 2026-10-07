function D = main_fixed_geometry_tip_refinement_diagnostic(varargin)
%MAIN_FIXED_GEOMETRY_TIP_REFINEMENT_DIAGNOSTIC
% Fixed-geometry sensitivity test for the audited h/2 near-tip core.
%
% Scientific isolation:
%   - accepted crack geometry is held exactly fixed at selected states;
%   - exterior grading remains the reference production law;
%   - rInner=0.10 DeltaA, rOuter=0.65 DeltaA, rCore=0.75 DeltaA unchanged;
%   - only structured core scale changes from L0=1 to L1=0.5;
%   - material, loading, solver, EDI, MTS, and all validity gates unchanged.
%
% Default fixed states:
%   P17 : positive KII/KI maximum;
%   P22 : first accepted negative KII/KI state.
%
% Safe default performs qualification/synthetic replay only:
%   D = main_fixed_geometry_tip_refinement_diagnostic();
%
% After inspecting qualification, explicitly authorize the two physical solves:
%   D = main_fixed_geometry_tip_refinement_diagnostic( ...
%       'AllowPhysicalSolves',true);

    ip=inputParser;
    addParameter(ip,'Segments',[17 22],@(x)isnumeric(x)&&isvector(x)&& ...
        all(isfinite(x))&&all(x==round(x))&&all(x>=2));
    addParameter(ip,'CoreScale',0.5,@(x)isnumeric(x)&&isscalar(x)&& ...
        isfinite(x)&&abs(x-.5)<=1e-14);
    addParameter(ip,'ReferenceEvidenceFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'FrozenStateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FastEDI',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PlotCandidates',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'SaveMeshFigures',true,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(root));

    evidenceFile=char(opt.ReferenceEvidenceFile);
    if isempty(evidenceFile)
        evidenceFile=fullfile(root,'paper','data','evidence_exact.mat');
    elseif ~local_is_absolute_path(evidenceFile)
        evidenceFile=fullfile(root,evidenceFile);
    end

    frozenFile=char(opt.FrozenStateFile);
    if isempty(frozenFile)
        paperCopy=fullfile(root,'paper','data','accepted_stage1_source.mat');
        standardCopy=fullfile(root,'verification','crack_path','stage1_starting_state.mat');
        if exist(paperCopy,'file')==2
            frozenFile=paperCopy;
        else
            frozenFile=standardCopy;
        end
    elseif ~local_is_absolute_path(frozenFile)
        frozenFile=fullfile(root,frozenFile);
    end

    outDir=char(opt.OutputDir);
    if isempty(outDir)
        outDir=fullfile(root,'verification','crack_path', ...
            'fixed_geometry_tip_refinement');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    assert(exist(evidenceFile,'file')==2,'tipref:MissingEvidence', ...
        'Committed audited evidence not found: %s',evidenceFile);
    assert(exist(frozenFile,'file')==2,'tipref:MissingFrozen', ...
        'Frozen Stage-I state not found: %s',frozenFile);

    ev=load(evidenceFile,'E','T','Q');
    assert(isfield(ev,'E')&&isstruct(ev.E)&& ...
        isfield(ev,'T')&&istable(ev.T)&& ...
        isfield(ev,'Q')&&istable(ev.Q), ...
        'tipref:BadEvidence', ...
        'evidence_exact.mat must contain E, T, and Q.');
    Eref=ev.E;
    Ref=ev.T;
    RefQual=ev.Q;
    assert(height(Ref)==23 && isequal(Ref.segment,(1:23)'), ...
        'tipref:ReferenceRows','Expected accepted reference rows P1--P23.');
    assert(size(Eref.vertices_m,1)==25 && size(Eref.vertices_m,2)==2, ...
        'tipref:ReferenceVertices', ...
        'Expected audited vertices P0--P24 in committed evidence.');
    assert(Eref.audit.completedPhysicalSegments==23 && ...
        Eref.audit.qualifiedUnsolvedSegment==24, ...
        'tipref:ReferenceClassification', ...
        'Committed evidence no longer has the accepted P23/P24 classification.');
    assert(all(ismember(opt.Segments,Ref.segment)), ...
        'tipref:IncompleteEvidence','Requested state is not physically accepted.');

    d=load(frozenFile,'R0');
    assert(isfield(d,'R0')&&isstruct(d.R0)&&d.R0.stage1Pass, ...
        'tipref:BadFrozen','Frozen MAT must contain passed R0.');
    R0=d.R0;

    mouthRef=Eref.vertices_m(1,:);
    mouthFrozen=[R0.summary.x_star_m,R0.summary.y_star_m];
    assert(norm(mouthRef-mouthFrozen)<=2e-12, ...
        'tipref:FrozenEvidenceMismatch', ...
        'Committed path evidence and frozen Stage-I mouth differ.');

    segments=unique(opt.Segments(:).','stable');

    fprintf('\n============================================================\n');
    fprintf('FIXED-GEOMETRY NEAR-TIP h/2 DIAGNOSTIC\n');
    fprintf('============================================================\n');
    fprintf('  states                 = %s\n',mat2str(segments));
    fprintf('  reference evidence     = %s\n',evidenceFile);
    fprintf('  frozen Stage I         = %s\n',frozenFile);
    fprintf('  physical solves        = %d\n',logical(opt.AllowPhysicalSolves));
    fprintf('  core scale L0 -> L1    = 1 -> %.3f\n',opt.CoreScale);
    fprintf('  hTip                    = 0.0270123254 -> 0.0135061627 mm\n');
    fprintf('  rInner / rOuter / rCore= 0.4 / 2.6 / 3.0 mm UNCHANGED\n');
    fprintf('  exterior mesh law      = REFERENCE, UNCHANGED\n');
    fprintf('  EDI formulation/radii  = UNCHANGED\n');
    fprintf('  all physical settings  = UNCHANGED\n');

    qRows=[];
    pRows=[];
    Qcell=cell(numel(segments),1);
    Rcell=cell(numel(segments),1);

    for ii=1:numel(segments)
        k=segments(ii);
        path=Eref.vertices_m(1:k+1,:);

        refState=Ref(Ref.segment==k,:);
        refQS=RefQual(RefQual.segment==k,:);
        assert(height(refState)==1 && height(refQS)==1 && logical(refQS.pass), ...
            'tipref:BadReference','Committed reference P%d is incomplete.',k);
        assert(abs(refState.crack_length_mm-4*k)<=2e-10, ...
            'tipref:ReferenceLength','Reference P%d crack length changed.',k);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('P%d: SAME accepted crack geometry, L1 h/2 tip core\n',k);
        fprintf('------------------------------------------------------------\n');

        candidateFile=fullfile(outDir,sprintf('L1_step_%03d_candidate.mat',k));
        qualFile=fullfile(outDir,sprintf('L1_step_%03d_qualification_small.mat',k));

        Q=qualify_incremental_crack_candidate(path, ...
            'FrozenState',R0, ...
            'CoreScale',opt.CoreScale, ...
            'RunSynthetic',opt.RunSynthetic, ...
            'FastEDI',opt.FastEDI, ...
            'ExteriorVerbose',true, ...
            'SaveCandidate',true, ...
            'CandidateFile',candidateFile, ...
            'SaveCompact',true, ...
            'CompactFile',qualFile, ...
            'Plot',opt.PlotCandidates);

        assert(Q.pass,'tipref:QualificationFailed', ...
            'L1 qualification failed at P%d. Do not run a physical solve.',k);
        Qcell{ii}=Q;

        if opt.PlotCandidates && opt.SaveMeshFigures
            fig=gcf;
            png=fullfile(outDir,sprintf('L1_step_%03d_mesh.png',k));
            try
                exportgraphics(fig,png,'Resolution',220);
                fprintf('  Mesh preview saved: %s\n',png);
            catch ME
                warning('tipref:FigureExport', ...
                    'Could not export P%d mesh preview: %s',k,ME.message);
            end
        end

        altQS=Q.summary;
        sc=Q.sampleCounts.nativePoints;
        qRows(end+1,:)=[ ...
            k, ...
            1e3*refQS.hTip_m,1e3*altQS.hTip_m,altQS.hTip_m/refQS.hTip_m, ...
            refQS.core_T3_elements,altQS.core_T3_elements, ...
            refQS.T3_elements,altQS.T3_elements, ...
            100*(altQS.T3_elements/refQS.T3_elements-1), ...
            refQS.EDI_elements,altQS.EDI_elements, ...
            refQS.min_angle_deg,altQS.min_angle_deg, ...
            refQS.max_neighbor_ratio,altQS.max_neighbor_ratio, ...
            sc(1),sc(2),sc(3),sc(4)]; %#ok<AGROW>

        if ~opt.AllowPhysicalSolves
            fprintf('  Qualification PASS at P%d. Physical solve remains guarded.\n',k);
            continue
        end

        cp=fullfile(outDir,sprintf('L1_step_%03d_physical_solved.mat',k));
        sf=fullfile(outDir,sprintf('L1_step_%03d_physical_small.mat',k));

        R=solve_incremental_crack_tip(Q.candidate, ...
            'FrozenState',R0, ...
            'AllowSolve',true, ...
            'FastEDI',opt.FastEDI, ...
            'CheckpointFile',cp, ...
            'SaveFile',sf);
        assert(R.pass,'tipref:PhysicalFailed', ...
            'L1 physical result failed at P%d.',k);
        Rcell{ii}=R;

        refKI=refState.KI_unit;
        refKII=refState.KII_unit;
        refQmix=refState.KII_over_KI;
        refTurn=refState.delta_theta_next_deg;

        altKI=R.EDI.KI_unit;
        altKII=R.EDI.KII_unit;
        altQmix=R.EDI.KII_over_KI;
        altTurn=R.deltaThetaNextDeg;

        pRows(end+1,:)=[ ...
            k, ...
            refKI,altKI,100*(altKI/refKI-1), ...
            refKII,altKII,100*(altKII/refKII-1), ...
            refQmix,altQmix,altQmix-refQmix,100*(altQmix/refQmix-1), ...
            refTurn,altTurn,altTurn-refTurn, ...
            refState.PCG_iterations,R.solverInfo.iter, ...
            R.solverInfo.relres,R.solverInfo.trueRelResidual]; %#ok<AGROW>
    end

    Qualification=array2table(qRows,'VariableNames',{ ...
        'segment', ...
        'hTip_ref_mm','hTip_L1_mm','hTip_ratio', ...
        'core_T3_ref','core_T3_L1', ...
        'T3_elements_ref','T3_elements_L1','T3_element_change_percent', ...
        'EDI_elements_ref','EDI_elements_L1', ...
        'min_angle_ref_deg','min_angle_L1_deg', ...
        'max_neighbor_ratio_ref','max_neighbor_ratio_L1', ...
        'native_w1','native_w2','native_w3','native_w4'});

    if isempty(pRows)
        Physical=table();
    else
        Physical=array2table(pRows,'VariableNames',{ ...
            'segment', ...
            'KI_ref','KI_L1','KI_change_percent', ...
            'KII_ref','KII_L1','KII_change_percent', ...
            'q_ref','q_L1','q_change_absolute','q_change_percent', ...
            'turn_ref_deg','turn_L1_deg','turn_change_deg', ...
            'PCG_iterations_ref','PCG_iterations_L1', ...
            'PCG_relres_L1','true_rel_residual_L1'});
    end

    fprintf('\nQUALIFICATION COMPARISON\n');
    disp(Qualification);
    if ~isempty(Physical)
        fprintf('\nPHYSICAL FIXED-GEOMETRY h/2 COMPARISON\n');
        disp(Physical);
    end

    D=struct();
    D.profile='L1_tip_h2_fixed_geometry';
    D.referenceEvidenceFile=evidenceFile;
    D.referenceSourceCommit=Eref.sourceCommit;
    D.frozenStateFile=frozenFile;
    D.outputDir=outDir;
    D.segments=segments;
    D.controls=struct( ...
        'coreScale',opt.CoreScale, ...
        'hTipRatio',opt.CoreScale, ...
        'pairedCoreRadiusChanged',false, ...
        'EDIRadiiChanged',false, ...
        'exteriorMeshLawChanged',false, ...
        'crackGeometryChanged',false, ...
        'physicalSolvesPerformed',logical(opt.AllowPhysicalSolves));
    D.qualification=Qualification;
    D.physical=Physical;
    D.qualifications=Qcell;
    D.physicalResults=Rcell;

    save(fullfile(outDir,'fixed_geometry_tip_refinement_diagnostic.mat'),'D','-v7');
    writetable(Qualification,fullfile(outDir,'qualification_comparison.csv'));
    if ~isempty(Physical)
        writetable(Physical,fullfile(outDir,'physical_comparison.csv'));
    end

    if ~opt.AllowPhysicalSolves
        fprintf('\nTIP-REFINEMENT DIAGNOSTIC PHASE 1 COMPLETE.\n');
        fprintf('  L1 h/2 meshes qualified on both fixed crack states.\n');
        fprintf('  No physical solve was performed.\n');
        fprintf('  Inspect mesh previews/tables before explicit physical authorization.\n');
    else
        fprintf('\nFIXED-GEOMETRY NEAR-TIP h/2 DIAGNOSTIC COMPLETE.\n');
        fprintf('  Crack geometries were held fixed exactly at all compared states.\n');
        fprintf('  Differences therefore isolate near-tip discretization sensitivity.\n');
    end
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
