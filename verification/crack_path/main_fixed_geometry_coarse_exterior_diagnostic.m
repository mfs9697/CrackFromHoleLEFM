function D = main_fixed_geometry_coarse_exterior_diagnostic(varargin)
%MAIN_FIXED_GEOMETRY_COARSE_EXTERIOR_DIAGNOSTIC
% Fixed-geometry sensitivity test for a coarser EXTERIOR mesh.
%
% The accepted crack geometry is held exactly fixed at selected states.
% The paired current-tip core, EDI radii, material, loading, solver, MTS
% rule, and COD windows are unchanged. Only the exterior target-size law is
% changed. This separates local SIF sensitivity to the exterior mesh from
% recursive trajectory sensitivity.
%
% Default states:
%   P17 : accepted positive KII/KI maximum
%   P22 : accepted first negative KII/KI state
%
% Safe default:
%   AllowPhysicalSolves = false
%
% Typical use:
%   Q = main_fixed_geometry_coarse_exterior_diagnostic();
%   % inspect qualification and mesh figures
%   D = main_fixed_geometry_coarse_exterior_diagnostic( ...
%       'AllowPhysicalSolves',true);
%
% The trial M1 exterior is deliberately NOT called accepted until all
% unchanged structural/synthetic/physical gates pass.

    ip=inputParser;
    addParameter(ip,'Segments',[17 22],@(x)isnumeric(x)&&isvector(x)&& ...
        all(isfinite(x))&&all(x==round(x))&&all(x>=2));
    addParameter(ip,'RunDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'FrozenStateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FastEDI',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PlotCandidates',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'SaveMeshFigures',true,@(x)islogical(x)&&isscalar(x));

    % Trial M1 exterior. The reference values are FarCap=0.625*DeltaA,
    % farSlope=0.10, boundaryMetricGrowth=0.25.
    addParameter(ip,'FarCapOverIncrement',1.25, ...
        @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'TransitionOverIncrement',1.0, ...
        @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'FarSlope',0.15, ...
        @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'BoundaryMetricGrowth',0.35, ...
        @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(root));

    runDir=char(opt.RunDir);
    if isempty(runDir)
        runDir=fullfile(root,'verification','crack_path','final_clean_run');
    elseif ~local_is_absolute_path(runDir)
        runDir=fullfile(root,runDir);
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
            'fixed_geometry_coarse_exterior');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    stateFile=fullfile(runDir,'path_run_state.mat');
    assert(exist(stateFile,'file')==2,'coarsemesh:MissingState', ...
        'Authoritative state not found: %s',stateFile);
    assert(exist(frozenFile,'file')==2,'coarsemesh:MissingFrozen', ...
        'Frozen Stage-I state not found: %s',frozenFile);

    s=load(stateFile,'State');
    assert(isfield(s,'State')&&isstruct(s.State), ...
        'coarsemesh:BadState','path_run_state.mat must contain State.');
    State=s.State;
    assert(State.completedPhysicalSegments>=max(opt.Segments), ...
        'coarsemesh:IncompleteState','Requested state is not physically accepted.');

    d=load(frozenFile,'R0');
    assert(isfield(d,'R0')&&isstruct(d.R0)&&d.R0.stage1Pass, ...
        'coarsemesh:BadFrozen','Frozen MAT must contain passed R0.');
    R0=d.R0;

    rowNames=cellstr(string(State.rowVariableNames));
    Ref=array2table(State.rowsThroughCompleted,'VariableNames',rowNames);
    assert(isequal(Ref.segment,(1:height(Ref))'), ...
        'coarsemesh:RowOrder','State rows are not segment ordered.');

    segments=unique(opt.Segments(:).','stable');
    calibration=struct( ...
        'farSlope',opt.FarSlope, ...
        'boundaryMetricGrowth',opt.BoundaryMetricGrowth);

    fprintf('\n============================================================\n');
    fprintf('FIXED-GEOMETRY COARSE-EXTERIOR DIAGNOSTIC\n');
    fprintf('============================================================\n');
    fprintf('  states                 = %s\n',mat2str(segments));
    fprintf('  reference path         = %s\n',stateFile);
    fprintf('  frozen Stage I         = %s\n',frozenFile);
    fprintf('  physical solves        = %d\n',logical(opt.AllowPhysicalSolves));
    fprintf('  tip core               = UNCHANGED\n');
    fprintf('  EDI radii              = UNCHANGED\n');
    fprintf('  exterior far cap       = %.6g * DeltaA\n',opt.FarCapOverIncrement);
    fprintf('  exterior transition    = %.6g * DeltaA\n',opt.TransitionOverIncrement);
    fprintf('  exterior far slope     = %.6g\n',opt.FarSlope);
    fprintf('  boundary metric growth = %.6g\n',opt.BoundaryMetricGrowth);
    fprintf('  all acceptance gates   = UNCHANGED\n');

    qRows=[];
    pRows=[];
    Qcell=cell(numel(segments),1);
    Rcell=cell(numel(segments),1);

    for ii=1:numel(segments)
        k=segments(ii);
        path=State.vertices(1:k+1,:);

        refQFile=fullfile(runDir,sprintf('step_%03d_qualification_small.mat',k));
        refPFile=fullfile(runDir,sprintf('step_%03d_physical_small.mat',k));
        assert(exist(refQFile,'file')==2,'coarsemesh:MissingReferenceQ', ...
            'Missing reference qualification: %s',refQFile);
        assert(exist(refPFile,'file')==2,'coarsemesh:MissingReferenceP', ...
            'Missing reference physical result: %s',refPFile);

        rq=load(refQFile,'Small'); refSmall=rq.Small;
        rp=load(refPFile,'R'); refR=rp.R;
        assert(refSmall.pass&&refR.pass,'coarsemesh:BadReference', ...
            'Reference P%d did not pass.',k);
        assert(norm(refR.pathFixed-path,'fro')<=2e-12, ...
            'coarsemesh:ReferencePath','Reference P%d path mismatch.',k);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('P%d: qualifying the SAME accepted crack geometry on trial M1\n',k);
        fprintf('------------------------------------------------------------\n');

        candidateFile=fullfile(outDir,sprintf('M1_step_%03d_candidate.mat',k));
        qualFile=fullfile(outDir,sprintf('M1_step_%03d_qualification_small.mat',k));

        Q=qualify_incremental_crack_candidate(path, ...
            'FrozenState',R0, ...
            'RunSynthetic',opt.RunSynthetic, ...
            'FastEDI',opt.FastEDI, ...
            'ExteriorVerbose',true, ...
            'ExteriorFarCapOverIncrement',opt.FarCapOverIncrement, ...
            'ExteriorTransitionOverIncrement',opt.TransitionOverIncrement, ...
            'ExteriorCalibration',calibration, ...
            'SaveCandidate',true, ...
            'CandidateFile',candidateFile, ...
            'SaveCompact',true, ...
            'CompactFile',qualFile, ...
            'Plot',opt.PlotCandidates);

        assert(Q.pass,'coarsemesh:QualificationFailed', ...
            'Trial M1 qualification failed at P%d. Do not run a physical solve.',k);
        Qcell{ii}=Q;

        if opt.PlotCandidates && opt.SaveMeshFigures
            fig=gcf;
            png=fullfile(outDir,sprintf('M1_step_%03d_mesh.png',k));
            try
                exportgraphics(fig,png,'Resolution',220);
                fprintf('  Mesh preview saved: %s\n',png);
            catch ME
                warning('coarsemesh:FigureExport', ...
                    'Could not export P%d mesh preview: %s',k,ME.message);
            end
        end

        refQS=refSmall.summary;
        altQS=Q.summary;
        qRows(end+1,:)=[ ...
            k, ...
            refQS.T3_nodes,altQS.T3_nodes, ...
            refQS.T3_elements,altQS.T3_elements, ...
            refQS.T6_nodes,altQS.T6_nodes, ...
            100*(altQS.T3_elements/refQS.T3_elements-1), ...
            refQS.min_angle_deg,altQS.min_angle_deg, ...
            refQS.max_neighbor_ratio,altQS.max_neighbor_ratio, ...
            refQS.EDI_elements,altQS.EDI_elements, ...
            refQS.physical_clearance_m*1e3,altQS.physical_clearance_m*1e3]; %#ok<AGROW>

        if ~opt.AllowPhysicalSolves
            fprintf('  Qualification PASS at P%d. Physical solve remains guarded.\n',k);
            continue
        end

        cp=fullfile(outDir,sprintf('M1_step_%03d_physical_solved.mat',k));
        sf=fullfile(outDir,sprintf('M1_step_%03d_physical_small.mat',k));

        R=solve_incremental_crack_tip(Q.candidate, ...
            'FrozenState',R0, ...
            'AllowSolve',true, ...
            'FastEDI',opt.FastEDI, ...
            'CheckpointFile',cp, ...
            'SaveFile',sf);
        assert(R.pass,'coarsemesh:PhysicalFailed', ...
            'Trial M1 physical result failed at P%d.',k);
        Rcell{ii}=R;

        refKI=refR.EDI.KI_unit;
        refKII=refR.EDI.KII_unit;
        refQmix=refR.EDI.KII_over_KI;
        refTurn=refR.deltaThetaNextDeg;

        altKI=R.EDI.KI_unit;
        altKII=R.EDI.KII_unit;
        altQmix=R.EDI.KII_over_KI;
        altTurn=R.deltaThetaNextDeg;

        pRows(end+1,:)=[ ...
            k, ...
            refKI,altKI,100*(altKI/refKI-1), ...
            refKII,altKII,100*(altKII/refKII-1), ...
            refQmix,altQmix,altQmix-refQmix, ...
            refTurn,altTurn,altTurn-refTurn, ...
            refR.solverInfo.iter,R.solverInfo.iter, ...
            R.solverInfo.relres,R.solverInfo.trueRelResidual]; %#ok<AGROW>
    end

    Qualification=array2table(qRows,'VariableNames',{ ...
        'segment', ...
        'T3_nodes_ref','T3_nodes_M1', ...
        'T3_elements_ref','T3_elements_M1', ...
        'T6_nodes_ref','T6_nodes_M1', ...
        'T3_element_change_percent', ...
        'min_angle_ref_deg','min_angle_M1_deg', ...
        'max_neighbor_ratio_ref','max_neighbor_ratio_M1', ...
        'EDI_elements_ref','EDI_elements_M1', ...
        'physical_clearance_ref_mm','physical_clearance_M1_mm'});

    if isempty(pRows)
        Physical=table();
    else
        Physical=array2table(pRows,'VariableNames',{ ...
            'segment', ...
            'KI_ref','KI_M1','KI_change_percent', ...
            'KII_ref','KII_M1','KII_change_percent', ...
            'q_ref','q_M1','q_change_absolute', ...
            'turn_ref_deg','turn_M1_deg','turn_change_deg', ...
            'PCG_iterations_ref','PCG_iterations_M1', ...
            'PCG_relres_M1','true_rel_residual_M1'});
    end

    fprintf('\nQUALIFICATION COMPARISON\n');
    disp(Qualification);
    if ~isempty(Physical)
        fprintf('\nPHYSICAL FIXED-GEOMETRY COMPARISON\n');
        disp(Physical);
    end

    D=struct();
    D.profile='M1_coarse_exterior_fixed_geometry';
    D.referenceRunDir=runDir;
    D.frozenStateFile=frozenFile;
    D.outputDir=outDir;
    D.segments=segments;
    D.controls=struct( ...
        'farCapOverIncrement',opt.FarCapOverIncrement, ...
        'transitionOverIncrement',opt.TransitionOverIncrement, ...
        'farSlope',opt.FarSlope, ...
        'boundaryMetricGrowth',opt.BoundaryMetricGrowth, ...
        'pairedCoreChanged',false, ...
        'EDIRadiiChanged',false, ...
        'crackGeometryChanged',false, ...
        'physicalSolvesPerformed',logical(opt.AllowPhysicalSolves));
    D.qualification=Qualification;
    D.physical=Physical;
    D.qualifications=Qcell;
    D.physicalResults=Rcell;

    save(fullfile(outDir,'fixed_geometry_coarse_exterior_diagnostic.mat'),'D','-v7');
    writetable(Qualification,fullfile(outDir,'qualification_comparison.csv'));
    if ~isempty(Physical)
        writetable(Physical,fullfile(outDir,'physical_comparison.csv'));
    end

    if ~opt.AllowPhysicalSolves
        fprintf('\nDIAGNOSTIC PHASE 1 COMPLETE.\n');
        fprintf('  Both trial meshes qualified. Inspect the mesh previews and table.\n');
        fprintf('  No physical solve was performed.\n');
        fprintf('  Rerun with ''AllowPhysicalSolves'',true only if the qualification is accepted.\n');
    else
        fprintf('\nFIXED-GEOMETRY COARSE-EXTERIOR DIAGNOSTIC COMPLETE.\n');
        fprintf('  Crack geometries were held fixed exactly at all compared states.\n');
        fprintf('  Differences therefore measure exterior-mesh sensitivity, not path divergence.\n');
    end
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
