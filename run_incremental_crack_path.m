function Path = run_incremental_crack_path(varargin)
%RUN_INCREMENTAL_CRACK_PATH
% Production incremental LEFM crack-path driver.
%
% Initialization:
%   theta_1 = 0 by prescription.
%   Accepted P1 SIFs seed theta_2 through MTS.
%
% Recurrence for k>=2:
%   qualify fixed path P0->...->Pk
%   solve exactly one physical state at Pk
%   Delta theta_{k+1} = MTS(KI_k,KII_k)
%   theta_{k+1} = theta_k + Delta theta_{k+1}
%   append one frozen-length segment if continuation is admissible.
%
% DEFAULT IS SAFE:
%   AllowPhysicalSolves = false
%
% Optional resume:
%   ResumeState     : in-memory State struct written by this driver
%   ResumeStateFile : MAT file containing State
%   ResumeSourceDir : optional directory containing prior compact physical
%                     results. If present, surviving rows are recovered.
%   ResumeAcceptedResult     : accepted in-memory R struct for the first
%                     unsolved resumed segment (for example R17u).
%   ResumeAcceptedResultFile : MAT file containing that struct as R.
%
% A resumed State is validated against the frozen Stage-I geometry, the
% prescribed 4-mm increments, stored segment angles, and accepted physical
% history before any new candidate or physical solve is attempted.
%
% The default MaxSegments=5 is intended as the first production regression
% run. Increase only after the short run passes.

    ip=inputParser;
    addParameter(ip,'FrozenState',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'MaxSegments',5,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=2&&x==round(x));
    addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FastEDI',false,@(x)islogical(x)&&isscalar(x));
    % Exterior-only controls. Defaults reproduce the accepted production
    % family; nondefault values define an alternative propagated mesh family.
    addParameter(ip,'ExteriorFarCapOverIncrement',0.625, ...
        @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'ExteriorTransitionOverIncrement',1.0, ...
        @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'ExteriorCalibration',struct(), ...
        @(x)isstruct(x)&&isscalar(x));
    addParameter(ip,'MeshFamilyLabel','reference', ...
        @(x)ischar(x)||isstring(x));
    addParameter(ip,'ReuseCandidates',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PlotEachStep',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'RegressionGates',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'StopAtCoreClearance',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'ResumeState',[],@(x)isempty(x)||(isstruct(x)&&isscalar(x)));
    addParameter(ip,'ResumeStateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'ResumeSourceDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'ResumeAcceptedResult',[],@(x)isempty(x)||(isstruct(x)&&isscalar(x)));
    addParameter(ip,'ResumeAcceptedResultFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'SeedKI',0.366479612185,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    addParameter(ip,'SeedKII',4.19826648316e-6,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    meshControls=struct( ...
        'farCapOverIncrement',opt.ExteriorFarCapOverIncrement, ...
        'transitionOverIncrement',opt.ExteriorTransitionOverIncrement, ...
        'calibrationOverride',opt.ExteriorCalibration, ...
        'label',char(opt.MeshFamilyLabel));
    meshControls.isReferenceProductionExterior = ...
        abs(opt.ExteriorTransitionOverIncrement-1.0)<=10*eps && ...
        abs(opt.ExteriorFarCapOverIncrement-0.625)<=10*eps && ...
        isempty(fieldnames(opt.ExteriorCalibration));

    incremental_profile_clock('begin','run',0);
    profileCleanup=onCleanup(@()incremental_profile_clock('end','run'));
    root=fileparts(mfilename('fullpath'));
    R0=local_load_frozen(root,opt.FrozenState);
    C=R0.C;
    S0=R0.summary(1,:);

    req={'a0_reserved_m','x_star_m','y_star_m', ...
        'nmat_x','nmat_y','that_x','that_y','stage1_pass'};
    local_require_table_variables(S0,req);
    assert(logical(S0.stage1_pass),'pathrun:Stage1NotPassed', ...
        'Frozen Stage-I state did not pass.');

    increment=S0.a0_reserved_m;
    mouth=[S0.x_star_m,S0.y_star_m];
    nMat=[S0.nmat_x,S0.nmat_y];nMat=nMat/norm(nMat);
    tHat=[S0.that_x,S0.that_y];tHat=tHat/norm(tHat);

    outDir=char(opt.OutputDir);
    if isempty(outDir)
        outDir=fullfile(root,'verification','crack_path','incremental_run');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    % --------------------------------------------------------------
    % Deterministic seed and accepted Stage III-D regression constants.
    % --------------------------------------------------------------
    theta1Deg=0;
    p0=mouth;
    p1=p0+increment*nMat;

    [deltaTheta2,theta2Deg]=kink_angle_LEFM_MTS(opt.SeedKI,opt.SeedKII);
    theta2=deltaTheta2; % theta_1=0, so absolute theta_2 equals first turn.
    e2=cos(theta2)*nMat+sin(theta2)*tHat;
    e2=e2/norm(e2);
    p2=p1+increment*e2;

    rows=nan(opt.MaxSegments,14);
    rows(1,:)=[1,p1,theta1Deg,opt.SeedKI,opt.SeedKII,opt.SeedKII/opt.SeedKI, ...
        theta2Deg,theta2Deg,NaN,NaN,NaN,NaN,0];
    seedRow=rows(1,:);

    regression=struct();
    regression.theta2_expected_deg=-0.00131272214162;
    regression.KI2_expected=0.437864119925;
    regression.KII2_expected=4.10524020214e-05;
    regression.deltaTheta3_expected_deg=-0.0107436495453;
    regression.theta3_expected_deg=-0.0120563716869;
    regression.theta2_pass=abs(theta2Deg-regression.theta2_expected_deg)<=5e-11;

    if opt.RegressionGates && ~regression.theta2_pass
        error('pathrun:SeedRegression', ...
            'Seed MTS theta_2 changed: got %.15g deg.',theta2Deg);
    end

    stepResults=cell(opt.MaxSegments,1);
    qualification=cell(opt.MaxSegments,1);

    [resumeState,resumeLabel,resumeSourceDir]=local_load_resume_state( ...
        root,opt.ResumeState,char(opt.ResumeStateFile),char(opt.ResumeSourceDir));

    if isempty(resumeState)
        vertices=[p0;p1;p2];
        thetaDeg=[theta1Deg;theta2Deg];
        startK=2;
        completedPhysicalSegments=1;
        resumed=false;
    else
        [vertices,thetaDeg,rows,stepResults,regression, ...
            completedPhysicalSegments,startK,resumeHistoryMode]=local_prepare_resume( ...
            resumeState,resumeSourceDir,opt.MaxSegments,rows,seedRow,stepResults, ...
            regression,p0,p1,p2,theta2Deg,increment,nMat,tHat,C, ...
            opt.RegressionGates,opt.StopAtCoreClearance,meshControls);
        resumed=true;
    end
    if ~resumed,resumeHistoryMode='fresh';end

    [resumeAcceptedResult,resumeAcceptedLabel]=local_load_resume_accepted_result( ...
        root,opt.ResumeAcceptedResult,char(opt.ResumeAcceptedResultFile));

    promotedAcceptedSegment=0;
    promotedNextThetaDeg=NaN;
    if ~isempty(resumeAcceptedResult)
        if ~resumed
            error('pathrun:AcceptedResultWithoutResume', ...
                'ResumeAcceptedResult requires ResumeState or ResumeStateFile.');
        end
        [vertices,thetaDeg,rows,stepResults,completedPhysicalSegments,startK, ...
            promotedAcceptedSegment,promotedNextThetaDeg]= ...
            local_promote_accepted_result( ...
                resumeAcceptedResult,opt.MaxSegments,vertices,thetaDeg,rows, ...
                stepResults,completedPhysicalSegments,startK,increment,nMat,tHat,C, ...
                opt.StopAtCoreClearance);
    end

    fprintf('\n============================================================\n');
    fprintf('GENERAL INCREMENTAL CRACK-PATH DRIVER\n');
    fprintf('============================================================\n');
    fprintf('  max segments      = %d\n',opt.MaxSegments);
    fprintf('  increment         = %.9f mm\n',1e3*increment);
    fprintf('  theta_1           = %+g deg (prescribed)\n',theta1Deg);
    fprintf('  accepted P1 KI    = %.12g MPa*sqrt(m)\n',opt.SeedKI);
    fprintf('  accepted P1 KII   = %+.12g MPa*sqrt(m)\n',opt.SeedKII);
    fprintf('  theta_2           = %+.12g deg (MTS from P1)\n',theta2Deg);
    fprintf('  physical solves   = %d\n',logical(opt.AllowPhysicalSolves));
    fprintf('  mesh family       = %s\n',meshControls.label);
    fprintf('  exterior far cap  = %.6g * DeltaA\n',meshControls.farCapOverIncrement);
    if ~isempty(fieldnames(meshControls.calibrationOverride))
        fprintf('  exterior override = enabled\n');
    end
    fprintf('  output directory  = %s\n',outDir);
    fprintf('  resume mode       = %d\n',resumed);
    if resumed
        fprintf('  resume state      = %s\n',resumeLabel);
        fprintf('  completed physical= %d segments\n',completedPhysicalSegments);
        fprintf('  first resumed step= %d\n',startK);
        fprintf('  resume history    = %s\n',resumeHistoryMode);
        if promotedAcceptedSegment>0
            fprintf('  promoted result   = P%d from %s\n', ...
                promotedAcceptedSegment,resumeAcceptedLabel);
            if startK<=opt.MaxSegments
                fprintf('  first new step    = %d\n',startK);
            else
                fprintf('  first new step    = none (target already accepted)\n');
            end
        end
    end

    stopReason='max_segments_reached';

    if promotedAcceptedSegment>0 && startK>opt.MaxSegments
        State=local_make_resume_state(vertices,thetaDeg,completedPhysicalSegments, ...
            promotedNextThetaDeg,regression,'max_segments_reached', ...
            rows,opt.FastEDI,outDir,meshControls);
        local_atomic_save_state(outDir,State);
    end

    % --------------------------------------------------------------
    % k = number of already existing finite crack segments.
    % A fresh run starts at k=2. A resumed run starts at the first
    % validated segment that does not yet have an accepted physical state.
    % --------------------------------------------------------------
    for k=startK:opt.MaxSegments
        incremental_profile_clock('begin','step',k);
        stepProfileCleanup=onCleanup(@()incremental_profile_clock('end','step'));
        pathNow=vertices(1:k+1,:);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('INCREMENTAL STEP k=%d / %d\n',k,opt.MaxSegments);
        fprintf('  current theta_k = %+.12g deg\n',thetaDeg(k));
        fprintf('  current tip     = [%.12g, %.12g] m\n', ...
            pathNow(end,1),pathNow(end,2));

        candidateFile=fullfile(outDir,sprintf('step_%03d_candidate.mat',k));
        qualFile=fullfile(outDir,sprintf('step_%03d_qualification_small.mat',k));

        incremental_profile_clock('phase','step','qualification_or_candidate_reuse');
        useSaved=false;
        if opt.ReuseCandidates && exist(candidateFile,'file')==2
            d=load(candidateFile,'candidate');
            if isfield(d,'candidate') && isstruct(d.candidate) && ...
                    isfield(d.candidate,'path') && ...
                    size(d.candidate.path,1)==size(pathNow,1) && ...
                    norm(d.candidate.path-pathNow,'fro')<=2e-12 && ...
                    isfield(d.candidate,'scientificallyReadyForIncrementalPhysicalSolve') && ...
                    logical(d.candidate.scientificallyReadyForIncrementalPhysicalSolve) && ...
                    local_candidate_mesh_controls_match(d.candidate,meshControls)
                candidate=d.candidate;
                useSaved=true;
                fprintf('  Reusing qualified candidate: %s\n',candidateFile);
            end
        end

        if ~useSaved
            Q=qualify_incremental_crack_candidate(pathNow, ...
                'FrozenState',R0, ...
                'RunSynthetic',opt.RunSynthetic, ...
                'FastEDI',opt.FastEDI, ...
                'ExteriorFarCapOverIncrement',opt.ExteriorFarCapOverIncrement, ...
                'ExteriorTransitionOverIncrement',opt.ExteriorTransitionOverIncrement, ...
                'ExteriorCalibration',opt.ExteriorCalibration, ...
                'SaveCandidate',true, ...
                'CandidateFile',candidateFile, ...
                'SaveCompact',true, ...
                'CompactFile',qualFile, ...
                'Plot',opt.PlotEachStep);
            assert(Q.pass,'pathrun:QualificationFailed', ...
                'Candidate qualification failed at segment %d.',k);
            candidate=Q.candidate;
            qualification{k}=Q.summary;
        else
            qualification{k}=table();
        end

        cp=fullfile(outDir,sprintf('step_%03d_physical_solved.mat',k));
        sf=fullfile(outDir,sprintf('step_%03d_physical_small.mat',k));

        incremental_profile_clock('phase','step','physical');
        R=solve_incremental_crack_tip(candidate, ...
            'FrozenState',R0, ...
            'AllowSolve',opt.AllowPhysicalSolves, ...
            'FastEDI',opt.FastEDI, ...
            'CheckpointFile',cp, ...
            'SaveFile',sf);

        if ~R.pass
            error('pathrun:PhysicalStepFailed', ...
                'Physical solve/postprocessing failed at segment %d.',k);
        end
        stepResults{k}=R;

        incremental_profile_clock('phase','step','recurrence_and_checkpoint');
        KI=R.EDI.KI_unit;
        KII=R.EDI.KII_unit;
        ratio=R.EDI.KII_over_KI;
        deltaNextDeg=R.deltaThetaNextDeg;
        thetaNextDeg=R.thetaNextDeg;

        rows(k,:)=[k,pathNow(end,:),thetaDeg(k),KI,KII,ratio, ...
            deltaNextDeg,thetaNextDeg,R.solverInfo.iter, ...
            R.solverInfo.relres,R.solverInfo.trueRelResidual, ...
            R.EDI.EDI_elements,R.newSolve];

        % Exact regression bridge to the accepted Stage III-D calculation.
        if k==2
            regression.KI2_pass=abs(KI-regression.KI2_expected)<=5e-10;
            regression.KII2_pass=abs(KII-regression.KII2_expected)<=5e-10;
            regression.deltaTheta3_pass= ...
                abs(deltaNextDeg-regression.deltaTheta3_expected_deg)<=5e-9;
            regression.theta3_pass= ...
                abs(thetaNextDeg-regression.theta3_expected_deg)<=5e-9;
            regression.step2_pass=regression.KI2_pass&&regression.KII2_pass&& ...
                regression.deltaTheta3_pass&&regression.theta3_pass;
            if opt.RegressionGates && ~regression.step2_pass
                error('pathrun:Stage3DRegression', ...
                    ['Generic driver did not reproduce accepted Stage III-D. ', ...
                     'KI=%g, KII=%g, Delta=%g deg, thetaNext=%g deg.'], ...
                     KI,KII,deltaNextDeg,thetaNextDeg);
            end
        end

        % Driver stops after recording the physical state at MaxSegments.
        % Save a solved-end resume state as well. On a later continuation
        % the next segment is reconstructed from thetaNextDeg without
        % repeating this accepted physical solve.
        if k==opt.MaxSegments
            State=local_make_resume_state(vertices,thetaDeg,k,thetaNextDeg, ...
                regression,'max_segments_reached',rows,opt.FastEDI,outDir,meshControls);
            local_atomic_save_state(outDir,State);
            clear stepProfileCleanup
            break
        end

        % Append exactly one new segment using the MTS-predicted ABSOLUTE
        % local direction in the frozen Stage-I frame.
        eNext=cosd(thetaNextDeg)*nMat+sind(thetaNextDeg)*tHat;
        eNext=eNext/norm(eNext);
        pNext=pathNow(end,:)+increment*eNext;

        if pNext(1)<0 || pNext(1)>C.A || pNext(2)<-C.B || pNext(2)>C.B
            stopReason='next_tip_outside_plate';
            fprintf('  STOP: proposed next tip is outside the plate.\n');
            State=local_make_resume_state(vertices,thetaDeg,k,thetaNextDeg, ...
                regression,stopReason,rows,opt.FastEDI,outDir,meshControls);
            local_atomic_save_state(outDir,State);
            clear stepProfileCleanup
            break
        end

        if opt.StopAtCoreClearance
            rCore=.75*increment;
            clearance=local_physical_clearance(pNext,C);
            if clearance<=rCore
                stopReason=sprintf('next_tip_core_clearance_%.6g_m',clearance);
                fprintf('  STOP: proposed next tip clearance %.6f mm <= rCore %.6f mm.\n', ...
                    1e3*clearance,1e3*rCore);
                State=local_make_resume_state(vertices,thetaDeg,k,thetaNextDeg, ...
                    regression,stopReason,rows,opt.FastEDI,outDir,meshControls);
                local_atomic_save_state(outDir,State);
                clear stepProfileCleanup
                break
            end
        end

        vertices(end+1,:)=pNext; %#ok<AGROW>
        thetaDeg(end+1,1)=thetaNextDeg; %#ok<AGROW>

        % Atomic compact path-state checkpoint after each accepted append.
        State=local_make_resume_state(vertices,thetaDeg,k,thetaNextDeg, ...
            regression,'running',rows,opt.FastEDI,outDir,meshControls);
        local_atomic_save_state(outDir,State);
        clear stepProfileCleanup
    end

    used=find(isfinite(rows(:,1)));
    StepTable=array2table(rows(used,:), ...
        'VariableNames',local_row_variable_names());
    StepTable.pass=true(height(StepTable),1);

    nExisting=size(vertices,1)-1;
    Path=struct();
    Path.vertices=vertices;
    Path.thetaDeg=thetaDeg;
    Path.increment=increment;
    Path.nSegments=nExisting;
    Path.stepTable=StepTable;
    Path.stepResults=stepResults;
    Path.qualification=qualification;
    Path.regression=regression;
    Path.stopReason=stopReason;
    Path.outputDir=outDir;
    Path.fastEDI=opt.FastEDI;
    Path.exteriorMeshControls=meshControls;
    Path.meshFamilyLabel=meshControls.label;
    Path.resumed=resumed;
    Path.resumeStateSource=resumeLabel;
    Path.resumeSourceDir=resumeSourceDir;
    Path.resumeStartSegment=startK;
    Path.resumeHistoryMode=resumeHistoryMode;
    Path.promotedAcceptedSegment=promotedAcceptedSegment;
    Path.resumeAcceptedResultSource=resumeAcceptedLabel;
    Path.completedPhysicalSegmentsAtStart=completedPhysicalSegments;
    Path.thirdAndLaterGenerated=nExisting>=3;
    Path.complete=strcmp(stopReason,'max_segments_reached') || ...
        startsWith(stopReason,'next_tip_');

    save(fullfile(outDir,'incremental_path_result.mat'),'Path','-v7');

    fprintf('\n============================================================\n');
    fprintf('INCREMENTAL PATH RUN COMPLETE\n');
    fprintf('============================================================\n');
    fprintf('  segments represented = %d\n',Path.nSegments);
    fprintf('  stop reason          = %s\n',Path.stopReason);
    if height(StepTable)>=2
        fprintf('  last solved theta    = %+.12g deg\n',StepTable.theta_deg(end));
        fprintf('  last KI/KII          = %.12g / %+.12g\n', ...
            StepTable.KI_unit(end),StepTable.KII_unit(end));
        fprintf('  next predicted theta = %+.12g deg\n',StepTable.theta_next_deg(end));
    end
end

function [State,label,sourceDir]=local_load_resume_state(root,Rin,fileIn,sourceIn)
    if ~isempty(Rin) && ~isempty(strtrim(fileIn))
        error('pathrun:ResumeSourceConflict', ...
            'Pass either ResumeState or ResumeStateFile, not both.');
    end

    State=[];
    label='<none>';
    sourceDir='';

    if ~isempty(Rin)
        State=Rin;
        label='<in-memory ResumeState>';
    elseif ~isempty(strtrim(fileIn))
        f=char(fileIn);
        if ~local_is_absolute_path(f),f=fullfile(root,f);end
        if exist(f,'file')~=2
            error('pathrun:MissingResumeState','ResumeStateFile not found: %s',f);
        end
        d=load(f,'State');
        if ~isfield(d,'State')||~isstruct(d.State)||~isscalar(d.State)
            error('pathrun:BadResumeStateFile', ...
                'ResumeStateFile must contain one scalar struct named State.');
        end
        State=d.State;
        label=f;
        sourceDir=fileparts(f);
    end

    if isempty(State),return,end

    if isempty(sourceDir) && isfield(State,'sourceOutputDir') && ...
            (ischar(State.sourceOutputDir)||isstring(State.sourceOutputDir))
        sourceDir=char(State.sourceOutputDir);
    end

    if ~isempty(strtrim(sourceIn))
        sourceDir=char(sourceIn);
        if ~local_is_absolute_path(sourceDir),sourceDir=fullfile(root,sourceDir);end
    end
end

function [R,label]=local_load_resume_accepted_result(root,Rin,fileIn)
    if ~isempty(Rin) && ~isempty(strtrim(fileIn))
        error('pathrun:AcceptedResultSourceConflict', ...
            'Pass either ResumeAcceptedResult or ResumeAcceptedResultFile, not both.');
    end
    R=[];
    label='<none>';
    if ~isempty(Rin)
        R=Rin;
        label='<in-memory accepted physical result>';
        return
    end
    if isempty(strtrim(fileIn)),return,end
    f=char(fileIn);
    if ~local_is_absolute_path(f),f=fullfile(root,f);end
    if exist(f,'file')~=2
        error('pathrun:MissingAcceptedResult', ...
            'ResumeAcceptedResultFile not found: %s',f);
    end
    d=load(f,'R');
    if ~isfield(d,'R')||~isstruct(d.R)||~isscalar(d.R)
        error('pathrun:BadAcceptedResultFile', ...
            'ResumeAcceptedResultFile must contain one scalar struct named R.');
    end
    R=d.R;
    label=f;
end

function [vertices,thetaDeg,rows,stepResults,kDone,startK,promotedK,nextThetaDeg]= ...
        local_promote_accepted_result(R,maxSegments,vertices,thetaDeg,rows, ...
        stepResults,kDone,startK,increment,nMat,tHat,C,stopAtCoreClearance)

    promotedK=startK;
    if promotedK~=kDone+1
        error('pathrun:AcceptedResultOrder', ...
            'Accepted result must correspond to the first unsolved resumed segment.');
    end
    if promotedK>maxSegments
        error('pathrun:AcceptedResultBeyondTarget', ...
            'Accepted result segment %d exceeds MaxSegments=%d.',promotedK,maxSegments);
    end
    if size(vertices,1)-1~=promotedK || numel(thetaDeg)~=promotedK
        error('pathrun:AcceptedResultPathShape', ...
            'Resume path must contain exactly the promoted unsolved segment.');
    end

    local_validate_resume_result(R,promotedK,vertices,thetaDeg);
    rows(promotedK,:)=local_row_from_result(R,promotedK);
    stepResults{promotedK}=R;
    kDone=promotedK;
    nextThetaDeg=R.thetaNextDeg;

    if kDone==maxSegments
        % Keep exactly the solved path. The schema-2 state records the next
        % MTS angle but does not append a segment beyond the requested target.
        startK=maxSegments+1;
        return
    end

    if ~isfinite(nextThetaDeg)
        error('pathrun:AcceptedResultNextAngle', ...
            'Accepted result has a nonfinite next MTS angle.');
    end

    eNext=cosd(nextThetaDeg)*nMat+sind(nextThetaDeg)*tHat;
    eNext=eNext/norm(eNext);
    pNext=vertices(end,:)+increment*eNext;

    if pNext(1)<0 || pNext(1)>C.A || pNext(2)<-C.B || pNext(2)>C.B
        error('pathrun:AcceptedResultNextTipOutside', ...
            'Accepted result predicts a next tip outside the plate.');
    end
    if stopAtCoreClearance
        rCore=.75*increment;
        clearance=local_physical_clearance(pNext,C);
        if clearance<=rCore
            error('pathrun:AcceptedResultNextTipClearance', ...
                'Accepted result predicts a next tip inside the core-clearance stop.');
        end
    end

    vertices(end+1,:)=pNext;
    thetaDeg(end+1,1)=nextThetaDeg;
    startK=kDone+1;
end

function [vertices,thetaDeg,rows,stepResults,regression,kDone,startK,historyMode]= ...
        local_prepare_resume(State,sourceDir,maxSegments,rows,seedRow,stepResults, ...
        regression,p0,p1,p2,theta2Deg,increment,nMat,tHat,C, ...
        regressionGates,stopAtCoreClearance,meshControls)

    req={'vertices','thetaDeg','completedPhysicalSegments','nextThetaDeg'};
    for j=1:numel(req)
        if ~isfield(State,req{j})||isempty(State.(req{j}))
            error('pathrun:ResumeField','Resume State missing %s.',req{j});
        end
    end

    if isfield(State,'exteriorMeshControls')
        if ~local_mesh_controls_equal(State.exteriorMeshControls,meshControls)
            error('pathrun:ResumeMeshFamily', ...
                'Resume State was created with different exterior mesh controls.');
        end
    elseif ~meshControls.isReferenceProductionExterior
        error('pathrun:ResumeMeshFamilyMissing', ...
            'Alternative mesh-family resume requires stored exteriorMeshControls.');
    end

    vertices=State.vertices;
    thetaDeg=State.thetaDeg(:);
    validateattributes(vertices,{'numeric'},{'2d','ncols',2,'finite'});
    validateattributes(thetaDeg,{'numeric'},{'column','finite'});

    nSeg=size(vertices,1)-1;
    kDone=State.completedPhysicalSegments;
    if ~isscalar(kDone)||~isfinite(kDone)||kDone~=round(kDone)||kDone<1
        error('pathrun:ResumeCompleted','Invalid completedPhysicalSegments.');
    end
    if numel(thetaDeg)~=nSeg
        error('pathrun:ResumeAngles','thetaDeg count does not match resumed segments.');
    end
    if ~(nSeg==kDone || nSeg==kDone+1)
        error('pathrun:ResumeShape', ...
            ['Resume State must contain either the solved path (N=kDone) ', ...
             'or exactly one already-appended unsolved segment (N=kDone+1).']);
    end
    if maxSegments<max(2,kDone+1)
        error('pathrun:ResumeMaxSegments', ...
            'MaxSegments=%d cannot continue after completed segment %d.', ...
            maxSegments,kDone);
    end

    seg=diff(vertices,1,1);
    segLength=vecnorm(seg,2,2);
    if max(abs(segLength-increment))>2e-12
        error('pathrun:ResumeIncrement','Resumed path changed the frozen increment.');
    end

    directions=seg./segLength;
    thetaGeomDeg=atan2d(directions*tHat(:),directions*nMat(:));
    if max(abs(thetaGeomDeg-thetaDeg))>1e-10
        error('pathrun:ResumeGeometryAngles', ...
            'Stored resume angles do not match the path geometry.');
    end

    if norm(vertices(1,:)-p0)>2e-12 || ...
            norm(vertices(2,:)-p1)>2e-12 || ...
            abs(thetaDeg(1))>5e-11
        error('pathrun:ResumeSeedGeometry', ...
            'Resumed mouth/first segment differs from the frozen prescription.');
    end
    if nSeg>=2 && (norm(vertices(3,:)-p2)>2e-12 || ...
            abs(thetaDeg(2)-theta2Deg)>5e-11)
        error('pathrun:ResumeTheta2', ...
            'Resumed second segment differs from the accepted P1 MTS seed.');
    end

    if abs(State.nextThetaDeg-thetaDeg(end))>1e-10 && nSeg==kDone+1
        error('pathrun:ResumeNextAngle', ...
            'Legacy appended resume state has inconsistent nextThetaDeg.');
    end

    % Recover as much accepted numeric history as is available. A legacy
    % path_run_state.mat is itself an atomic accepted-step checkpoint, so
    % missing old per-step result files must not make the path unusable.
    historyMode='checkpoint_only';
    hasEmbeddedShape=isfield(State,'rowsThroughCompleted') && ...
        isnumeric(State.rowsThroughCompleted) && ...
        size(State.rowsThroughCompleted,2)==size(rows,2) && ...
        size(State.rowsThroughCompleted,1)>=kDone;

    if hasEmbeddedShape
        hist=State.rowsThroughCompleted(1:kDone,:);
        for kk=2:kDone
            if isfinite(hist(kk,1))
                if hist(kk,1)~=kk
                    error('pathrun:ResumeRows', ...
                        'Embedded resume row index %g is inconsistent at segment %d.', ...
                        hist(kk,1),kk);
                end
                rows(kk,:)=hist(kk,:);
            end
        end
        historyMode='embedded_partial';
    end

    % Opportunistically hydrate any missing historical rows from surviving
    % compact physical results. ResumeSourceDir is an aid, not a dependency.
    if kDone>=2 && ~isempty(sourceDir) && exist(sourceDir,'dir')==7
        hydratedAny=false;
        for kk=2:kDone
            if local_row_complete(rows(kk,:)),continue,end
            f=fullfile(sourceDir,sprintf('step_%03d_physical_small.mat',kk));
            if exist(f,'file')~=2,continue,end
            d=load(f,'R');
            if ~isfield(d,'R')||~isstruct(d.R)
                error('pathrun:ResumeHistoryBad','%s does not contain struct R.',f);
            end
            R=d.R;
            local_validate_resume_result(R,kk,vertices,thetaDeg);
            rows(kk,:)=local_row_from_result(R,kk);
            stepResults{kk}=R;
            hydratedAny=true;
        end
        if hydratedAny,historyMode='mixed_or_external';end
    end

    % Never trust a stored seed row over the current frozen P1 seed.
    rows(1,:)=seedRow;

    historyComplete=true;
    if kDone>=2
        historyComplete=all(arrayfun(@(kk)local_row_complete(rows(kk,:)),2:kDone));
    end
    if historyComplete
        historyMode='complete';
    end

    % P2 is the immutable regression bridge. Prefer its numeric row when
    % available. Otherwise validate the regression record stored in the
    % atomic legacy checkpoint.
    if kDone>=2 && local_row_complete(rows(2,:))
        KI=rows(2,5);KII=rows(2,6);
        deltaNextDeg=rows(2,8);thetaNextDeg=rows(2,9);
        regression.KI2_pass=abs(KI-regression.KI2_expected)<=5e-10;
        regression.KII2_pass=abs(KII-regression.KII2_expected)<=5e-10;
        regression.deltaTheta3_pass= ...
            abs(deltaNextDeg-regression.deltaTheta3_expected_deg)<=5e-9;
        regression.theta3_pass= ...
            abs(thetaNextDeg-regression.theta3_expected_deg)<=5e-9;
        regression.step2_pass=regression.KI2_pass&&regression.KII2_pass&& ...
            regression.deltaTheta3_pass&&regression.theta3_pass;
    elseif kDone>=2
        local_validate_legacy_regression(State,regression);
        regression.theta2_pass=true;
        regression.KI2_pass=true;
        regression.KII2_pass=true;
        regression.deltaTheta3_pass=true;
        regression.theta3_pass=true;
        regression.step2_pass=true;
        historyMode='checkpoint_only';
    end

    if kDone>=2 && regressionGates && ~regression.step2_pass
        error('pathrun:ResumeStage3DRegression', ...
            'Resume history does not reproduce accepted Stage III-D.');
    end

    % New solved-end states contain no next leg yet. Reconstruct exactly one
    % proposed segment from the stored MTS angle before entering the loop.
    if nSeg==kDone
        thetaNextDeg=State.nextThetaDeg;
        if ~isfinite(thetaNextDeg)
            error('pathrun:ResumeNextAngle','Solved-end nextThetaDeg is not finite.');
        end
        eNext=cosd(thetaNextDeg)*nMat+sind(thetaNextDeg)*tHat;
        eNext=eNext/norm(eNext);
        pNext=vertices(end,:)+increment*eNext;

        if pNext(1)<0 || pNext(1)>C.A || pNext(2)<-C.B || pNext(2)>C.B
            error('pathrun:ResumeNextTipOutside', ...
                'Saved solved-end state predicts a next tip outside the plate.');
        end
        if stopAtCoreClearance
            rCore=.75*increment;
            clearance=local_physical_clearance(pNext,C);
            if clearance<=rCore
                error('pathrun:ResumeNextTipClearance', ...
                    'Saved solved-end next tip violates the current core-clearance gate.');
            end
        end

        vertices(end+1,:)=pNext;
        thetaDeg(end+1,1)=thetaNextDeg;
        nSeg=nSeg+1;
    end

    startK=kDone+1;
    if nSeg~=startK
        error('pathrun:ResumeInternal', ...
            'Validated resume state did not produce exactly one unsolved segment.');
    end
end

function tf=local_row_complete(row)
    tf=isnumeric(row)&&isvector(row)&&numel(row)==14 && ...
        all(isfinite(row));
end

function local_validate_legacy_regression(State,current)
    if ~isfield(State,'regression')||~isstruct(State.regression)|| ...
            ~isscalar(State.regression)
        error('pathrun:ResumeRegressionMissing', ...
            ['Legacy checkpoint has no complete numeric P2 history and no ', ...
             'stored regression record. It cannot be resumed safely.']);
    end

    old=State.regression;
    req={'theta2_expected_deg','KI2_expected','KII2_expected', ...
        'deltaTheta3_expected_deg','theta3_expected_deg','step2_pass'};
    for j=1:numel(req)
        if ~isfield(old,req{j})||isempty(old.(req{j}))
            error('pathrun:ResumeRegressionMissing', ...
                'Legacy regression record is missing %s.',req{j});
        end
    end

    if abs(old.theta2_expected_deg-current.theta2_expected_deg)>5e-11 || ...
            abs(old.KI2_expected-current.KI2_expected)>5e-10 || ...
            abs(old.KII2_expected-current.KII2_expected)>5e-10 || ...
            abs(old.deltaTheta3_expected_deg-current.deltaTheta3_expected_deg)>5e-9 || ...
            abs(old.theta3_expected_deg-current.theta3_expected_deg)>5e-9 || ...
            ~logical(old.step2_pass)
        error('pathrun:ResumeRegressionMismatch', ...
            'Legacy checkpoint regression does not match the accepted Stage III-D bridge.');
    end

    passNames={'theta2_pass','KI2_pass','KII2_pass', ...
        'deltaTheta3_pass','theta3_pass'};
    for j=1:numel(passNames)
        if isfield(old,passNames{j}) && ~logical(old.(passNames{j}))
            error('pathrun:ResumeRegressionMismatch', ...
                'Legacy checkpoint records failed regression gate %s.',passNames{j});
        end
    end
end

function local_validate_resume_result(R,k,vertices,thetaDeg)
    req={'pass','nSegments','pathFixed','thetaCurrentDeg','EDI', ...
        'solverInfo','gates','newSolve','deltaThetaNextDeg','thetaNextDeg'};
    for j=1:numel(req)
        if ~isfield(R,req{j})||isempty(R.(req{j}))
            error('pathrun:ResumeResultField', ...
                'Accepted compact result for segment %d is missing %s.',k,req{j});
        end
    end
    if ~logical(R.pass)||R.nSegments~=k || ...
            norm(R.pathFixed-vertices(1:k+1,:),'fro')>2e-12 || ...
            abs(R.thetaCurrentDeg-thetaDeg(k))>1e-10 || ...
            ~all(structfun(@logical,R.gates)) || ...
            R.EDI.EDI_elements(1)~=11316
        error('pathrun:ResumeResultMismatch', ...
            'Accepted compact physical result for segment %d does not match resume path.',k);
    end
end

function row=local_row_from_result(R,k)
    row=[k,R.pathFixed(end,:),R.thetaCurrentDeg, ...
        R.EDI.KI_unit(1),R.EDI.KII_unit(1),R.EDI.KII_over_KI(1), ...
        R.deltaThetaNextDeg,R.thetaNextDeg,R.solverInfo.iter, ...
        R.solverInfo.relres,R.solverInfo.trueRelResidual, ...
        R.EDI.EDI_elements(1),logical(R.newSolve)];
end

function State=local_make_resume_state(vertices,thetaDeg,kDone,nextThetaDeg, ...
        regression,stopReason,rows,fastEDI,outDir,meshControls)
    State=struct();
    State.schemaVersion=2;
    State.vertices=vertices;
    State.thetaDeg=thetaDeg;
    State.completedPhysicalSegments=kDone;
    State.nextThetaDeg=nextThetaDeg;
    State.regression=regression;
    State.stopReason=stopReason;
    State.rowsThroughCompleted=rows(1:kDone,:);
    State.rowVariableNames=local_row_variable_names();
    if kDone<=1
        State.rowHistoryComplete=true;
    else
        State.rowHistoryComplete=all(arrayfun( ...
            @(kk)local_row_complete(rows(kk,:)),2:kDone));
    end
    State.fastEDI=logical(fastEDI);
    State.sourceOutputDir=outDir;
    State.exteriorMeshControls=meshControls;
    State.meshFamilyLabel=meshControls.label;
end

function local_atomic_save_state(outDir,State)
    tmp=fullfile(outDir,'path_run_state.incomplete.mat');
    dst=fullfile(outDir,'path_run_state.mat');
    save(tmp,'State','-v7');
    [ok,msg]=movefile(tmp,dst,'f');
    if ~ok,error('pathrun:StateSave','%s',msg);end
end

function names=local_row_variable_names()
    names={'segment','tip_x_m','tip_y_m','theta_deg', ...
        'KI_unit','KII_unit','KII_over_KI','delta_theta_next_deg', ...
        'theta_next_deg','PCG_iterations','PCG_relres','true_rel_residual', ...
        'EDI_elements','newSolve'};
end

function R0=local_load_frozen(root,Rin)
    if ~isempty(Rin),R0=Rin;return,end
    f=fullfile(root,'verification','crack_path','stage1_starting_state.mat');
    if exist(f,'file')~=2
        error('pathrun:MissingFrozenState', ...
            'Frozen Stage-I state is missing; pass ''FrozenState'',R0.');
    end
    d=load(f);
    assert(isfield(d,'R0')&&isstruct(d.R0),'pathrun:BadFrozenState');
    R0=d.R0;
end

function d=local_physical_clearance(x,C)
    d=min([x(1),C.A-x(1),x(2)+C.B,C.B-x(2)]);
    if isfield(C,'holes')&&~isempty(C.holes)
        holes=C.holes;
        if isstruct(holes),holes=num2cell(holes);end
        for k=1:numel(holes)
            h=holes{k};
            if strcmpi(strtrim(h.type),'circle')
                d=min(d,abs(norm(x-h.center(:).')-h.r));
            end
        end
    end
end

function local_require_table_variables(T,names)
    miss=names(~ismember(names,T.Properties.VariableNames));
    if ~isempty(miss)
        error('pathrun:FrozenFields','Missing frozen fields: %s',strjoin(miss,', '));
    end
end

function tf=local_candidate_mesh_controls_match(candidate,requested)
    tf=false;
    if ~isfield(candidate,'exteriorMeshControls') || ...
            ~isstruct(candidate.exteriorMeshControls)
        % Historical reference candidates predate explicit control metadata.
        tf=requested.isReferenceProductionExterior;
        return
    end
    tf=local_mesh_controls_equal(candidate.exteriorMeshControls,requested);
end

function tf=local_mesh_controls_equal(a,b)
    tf=false;
    if ~isstruct(a)||~isstruct(b),return,end

    if isfield(a,'farCapOverIncrement')
        fa=a.farCapOverIncrement;
    elseif isfield(a,'farCapOverA0')
        fa=a.farCapOverA0;
    else
        return
    end
    if isfield(a,'transitionOverIncrement')
        ta=a.transitionOverIncrement;
    elseif isfield(a,'transitionOverA0')
        ta=a.transitionOverA0;
    else
        return
    end

    if ~isfield(b,'farCapOverIncrement')||~isfield(b,'transitionOverIncrement')
        return
    end

    if abs(fa-b.farCapOverIncrement)>1e-14 || ...
            abs(ta-b.transitionOverIncrement)>1e-14
        return
    end

    if isfield(a,'calibrationOverride')
        ca=a.calibrationOverride;
    else
        ca=struct();
    end
    if isfield(b,'calibrationOverride')
        cb=b.calibrationOverride;
    else
        cb=struct();
    end
    tf=isequaln(ca,cb);
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
