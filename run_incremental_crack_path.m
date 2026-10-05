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
% The default MaxSegments=5 is intended as the first production regression
% run. Increase only after the short run passes.

    ip=inputParser;
    addParameter(ip,'FrozenState',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'MaxSegments',5,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=2&&x==round(x));
    addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'RunSynthetic',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ReuseCandidates',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PlotEachStep',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'RegressionGates',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'StopAtCoreClearance',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'SeedKI',0.366479612185,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    addParameter(ip,'SeedKII',4.19826648316e-6,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(mfilename('fullpath'));
    local_assert_branch(root);
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
    % Initial finite segment: prescribed radial/material-normal direction.
    % --------------------------------------------------------------
    theta1Deg=0;
    p0=mouth;
    p1=p0+increment*nMat;

    [deltaTheta2,theta2Deg]=kink_angle_LEFM_MTS(opt.SeedKI,opt.SeedKII);
    theta2=deltaTheta2; % theta_1=0, so absolute theta_2 equals first turn.
    e2=cos(theta2)*nMat+sin(theta2)*tHat;
    e2=e2/norm(e2);
    p2=p1+increment*e2;

    vertices=[p0;p1;p2];
    thetaDeg=[theta1Deg;theta2Deg];

    % Accepted seed is itself a physical tip state at P1.
    rows=nan(opt.MaxSegments,14);
    rows(1,:)=[1,p1,theta1Deg,opt.SeedKI,opt.SeedKII,opt.SeedKII/opt.SeedKI, ...
        theta2Deg,theta2Deg,NaN,NaN,NaN,NaN,0,1];

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
    fprintf('  output directory  = %s\n',outDir);

    stopReason='max_segments_reached';
    stepResults=cell(opt.MaxSegments,1);
    qualification=cell(opt.MaxSegments,1);

    % --------------------------------------------------------------
    % k = number of already existing finite crack segments.
    % Start at k=2 because the accepted P1 state seeded segment 2.
    % --------------------------------------------------------------
    for k=2:opt.MaxSegments
        pathNow=vertices(1:k+1,:);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('INCREMENTAL STEP k=%d / %d\n',k,opt.MaxSegments);
        fprintf('  current theta_k = %+.12g deg\n',thetaDeg(k));
        fprintf('  current tip     = [%.12g, %.12g] m\n', ...
            pathNow(end,1),pathNow(end,2));

        candidateFile=fullfile(outDir,sprintf('step_%03d_candidate.mat',k));
        qualFile=fullfile(outDir,sprintf('step_%03d_qualification_small.mat',k));

        useSaved=false;
        if opt.ReuseCandidates && exist(candidateFile,'file')==2
            d=load(candidateFile,'candidate');
            if isfield(d,'candidate') && isstruct(d.candidate) && ...
                    isfield(d.candidate,'path') && ...
                    size(d.candidate.path,1)==size(pathNow,1) && ...
                    norm(d.candidate.path-pathNow,'fro')<=2e-12 && ...
                    isfield(d.candidate,'scientificallyReadyForIncrementalPhysicalSolve') && ...
                    logical(d.candidate.scientificallyReadyForIncrementalPhysicalSolve)
                candidate=d.candidate;
                useSaved=true;
                fprintf('  Reusing qualified candidate: %s\n',candidateFile);
            end
        end

        if ~useSaved
            Q=qualify_incremental_crack_candidate(pathNow, ...
                'FrozenState',R0, ...
                'RunSynthetic',opt.RunSynthetic, ...
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

        R=solve_incremental_crack_tip(candidate, ...
            'FrozenState',R0, ...
            'AllowSolve',opt.AllowPhysicalSolves, ...
            'CheckpointFile',cp, ...
            'SaveFile',sf);

        if ~R.pass
            error('pathrun:PhysicalStepFailed', ...
                'Physical solve/postprocessing failed at segment %d.',k);
        end
        stepResults{k}=R;

        KI=R.EDI.KI_unit;
        KII=R.EDI.KII_unit;
        ratio=R.EDI.KII_over_KI;
        deltaNextDeg=R.deltaThetaNextDeg;
        thetaNextDeg=R.thetaNextDeg;

        rows(k,:)=[k,pathNow(end,:),thetaDeg(k),KI,KII,ratio, ...
            deltaNextDeg,thetaNextDeg,R.solverInfo.iter, ...
            R.solverInfo.relres,R.solverInfo.trueRelResidual, ...
            R.EDI.EDI_elements,R.newSolve,R.pass];

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
        if k==opt.MaxSegments
            break
        end

        % Append exactly one new segment using the MTS-predicted ABSOLUTE
        % local direction in the frozen Stage-I frame.
        eNext=cosd(thetaNextDeg)*nMat+sind(thetaNextDeg)*tHat;
        eNext=eNext/norm(eNext);
        pNext=pathNow(end,:)+increment*eNext;

        if opt.StopAtCoreClearance
            rCore=.75*increment;
            clearance=local_physical_clearance(pNext,C);
            if clearance<=rCore
                stopReason=sprintf('next_tip_core_clearance_%.6g_m',clearance);
                fprintf('  STOP: proposed next tip clearance %.6f mm <= rCore %.6f mm.\n', ...
                    1e3*clearance,1e3*rCore);
                break
            end
        end

        if pNext(1)<0 || pNext(1)>C.A || pNext(2)<-C.B || pNext(2)>C.B
            stopReason='next_tip_outside_plate';
            fprintf('  STOP: proposed next tip is outside the plate.\n');
            break
        end

        vertices(end+1,:)=pNext; %#ok<AGROW>
        thetaDeg(end+1,1)=thetaNextDeg; %#ok<AGROW>

        % Atomic compact path-state checkpoint after each accepted append.
        State=struct('vertices',vertices,'thetaDeg',thetaDeg, ...
            'completedPhysicalSegments',k,'nextThetaDeg',thetaNextDeg, ...
            'regression',regression,'stopReason','running');
        tmp=fullfile(outDir,'path_run_state.incomplete.mat');
        dst=fullfile(outDir,'path_run_state.mat');
        save(tmp,'State','-v7');
        [ok,msg]=movefile(tmp,dst,'f');
        if ~ok,error('pathrun:StateSave','%s',msg);end
    end

    used=find(isfinite(rows(:,1)));
    StepTable=array2table(rows(used,:), ...
        'VariableNames',{'segment','tip_x_m','tip_y_m','theta_deg', ...
        'KI_unit','KII_unit','KII_over_KI','delta_theta_next_deg', ...
        'theta_next_deg','PCG_iterations','PCG_relres','true_rel_residual', ...
        'EDI_elements','newSolve'});
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

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end

function local_assert_branch(root)
    [status,b]=system(sprintf('git -C "%s" branch --show-current',root));
    assert(status==0&&strcmp(strtrim(b),'incremental-general-crack-path'), ...
        'pathrun:Branch', ...
        'Run the general path driver only on incremental-general-crack-path.');
end
