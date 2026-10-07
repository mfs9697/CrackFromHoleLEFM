function R = main_stage2_theta0_physical_solve(varargin)
%MAIN_STAGE2_THETA0_PHYSICAL_SOLVE
% One guarded physical LEFM solve for the qualified full-domain theta_1=0
% candidate. This is NOT an angle sweep.
%
% DEFAULT IS SAFE:
%   AllowSolve = false
%
% Required scientific lineage:
%   frozen Stage-I initiation state
%     -> qualified 4-mm theta_1=0 crack geometry
%     -> qualified a0-scaled paired core
%     -> qualified full-domain core embedding
%     -> THIS DRIVER: exactly one unit-load physical solve
%
% Solver:
%   free-DOF SPD formulation
%   symamd permutation
%   parameter-free symmetric Gauss-Seidel preconditioner
%   pcg tolerance 1e-10, maxit 5000
%
% Postprocessing happens ONLY after the physical field is checkpointed:
%   - native COD on the four predeclared r/a0 windows;
%   - exactly one matched 16-point FE-nodal-q interaction EDI with
%       r_inner/a0 = 0.10
%       r_outer/a0 = 0.65
%
% No physical K target is imposed. In particular, KII is measured, not
% compared against an expected sign or magnitude.
%
% Recommended call while R0 and F from Stage II-C remain in memory:
%
%   Rphys = main_stage2_theta0_physical_solve( ...
%       'FrozenState',R0, ...
%       'Candidate',F.candidate, ...
%       'AllowSolve',true);
%
% A valid existing checkpoint is reused without a second solve.

    ip=inputParser;
    addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FrozenState',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'Candidate',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'CandidateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'CheckpointFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'SaveFile','',@(x)ischar(x)||isstring(x));
    % Explicit escape hatch for a separately qualified nonreference exterior
    % mesh. Default behavior retains the historical branch/count fingerprints.
    addParameter(ip,'AlternativeQualifiedCandidate',false, ...
        @(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(mfilename('fullpath'));
    if ~opt.AlternativeQualifiedCandidate
        local_assert_branch(root);
    end
    outDir=fullfile(root,'verification','crack_path');
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    cp=char(opt.CheckpointFile);
    if isempty(cp)
        cp=fullfile(outDir,'stage2_theta0_physical_solved.mat');
    end
    saveFile=char(opt.SaveFile);
    if isempty(saveFile)
        saveFile=fullfile(outDir,'stage2_theta0_physical_small_data.mat');
    end

    candidateFile=char(opt.CandidateFile);
    if isempty(candidateFile)
        candidateFile=fullfile(outDir, ...
            'stage2_full_domain_scaled_core_theta0_candidate_T3.mat');
    end

    [R0,frozenSource]=local_load_frozen(root,opt.FrozenState);
    [candidate,candidateSource]=local_load_candidate( ...
        opt.Candidate,candidateFile);

    for name={'pcg','symamd'}
        if ~local_matlab_callable(name{1})
            error('stage2phys:MissingIterativeTool', ...
                'Physical solve requires callable MATLAB %s.',name{1});
        end
    end

    % ------------------------------------------------------------------
    % Frozen physical configuration and exact candidate provenance.
    % ------------------------------------------------------------------
    assert(isfield(R0,'summary')&&istable(R0.summary)&&height(R0.summary)==1, ...
        'stage2phys:FrozenSummary','Frozen R0.summary is required.');
    assert(isfield(R0,'C')&&isstruct(R0.C), ...
        'stage2phys:FrozenConfig','Frozen R0.C is required.');

    T0=R0.summary;
    req={'a0_reserved_m','x_star_m','y_star_m', ...
        'nmat_x','nmat_y','stage1_pass'};
    local_require_table_variables(T0,req);
    row=T0(1,:);
    assert(logical(row.stage1_pass),'stage2phys:Stage1NotPassed', ...
        'Frozen Stage-I state did not pass.');

    C=R0.C;
    a0=row.a0_reserved_m;
    mouth=[row.x_star_m,row.y_star_m];
    e1=[row.nmat_x,row.nmat_y];e1=e1/norm(e1);
    tip=mouth+a0*e1;

    local_require_candidate(candidate);

    if ~candidate.scientificallyReadyForOnePhysicalSolve
        error('stage2phys:CandidateNotQualified', ...
            'Full-domain candidate is not marked ready for one physical solve.');
    end
    if ~all(structfun(@logical,candidate.gates)) || ...
            ~all(structfun(@logical,candidate.syntheticGates))
        error('stage2phys:CandidateGateFailure', ...
            'Full-domain candidate does not retain all structural/synthetic gates.');
    end

    P=candidate.p;
    T=candidate.t;
    crack=candidate.crack;

    if ~opt.AlternativeQualifiedCandidate
        assert(size(P,1)==24389 && size(T,1)==47828, ...
            'stage2phys:CandidateT3Fingerprint', ...
            'Qualified theta0 T3 count fingerprint changed.');
    else
        assert(isfield(candidate,'exteriorMeshControls') && ...
            isstruct(candidate.exteriorMeshControls) && ...
            isfield(candidate.exteriorMeshControls,'isReferenceProductionExterior') && ...
            ~logical(candidate.exteriorMeshControls.isReferenceProductionExterior), ...
            'stage2phys:AlternativeCandidateProvenance', ...
            'AlternativeQualifiedCandidate requires explicit nonreference exterior provenance.');
    end
    assert(abs(candidate.a0-a0)<=1e-14 && candidate.theta1==0, ...
        'stage2phys:CandidateCrackParameter', ...
        'Candidate a0/theta1 differs from qualified theta0 problem.');
    assert(norm(crack.Pmid(1,:)-mouth)<=2e-12 && ...
           norm(crack.Pmid(end,:)-tip)<=2e-12, ...
        'stage2phys:CandidateCrackGeometry', ...
        'Candidate mouth/tip differs from frozen theta0 geometry.');
    assert(numel(candidate.pairedElementIDs)==12678, ...
        'stage2phys:CoreFingerprint','Qualified core element count changed.');
    assert(numel(candidate.primarySupportElementIDs)==11316, ...
        'stage2phys:SupportFingerprint','Qualified EDI support count changed.');

    [P6,T6]=T3toT6_fast(P,T);
    if ~opt.AlternativeQualifiedCandidate
        assert(size(P6,1)==96606, ...
            'stage2phys:T6Fingerprint','Qualified theta0 T6 node count changed.');
    end

    mesh=struct('coord3',P,'connect3',T,'coord',P6,'connect',T6);

    assert(isfield(candidate.mat,'E')&&isfield(candidate.mat,'nu')&& ...
           isfield(candidate.mat,'ps') && ...
           abs(candidate.mat.E-C.E)<=1e-12*max(1,abs(C.E)) && ...
           abs(candidate.mat.nu-C.nu)<=1e-14 && candidate.mat.ps==C.ps, ...
        'stage2phys:CandidateMaterialMismatch', ...
        'Qualified candidate material differs from frozen Stage-I material.');

    mat=local_material_with_D(candidate.mat,C);

    % Physics must be the exact frozen Stage-I material/loading/anchoring.
    assert(abs(mat.E-C.E)<=1e-12*max(1,abs(C.E)) && ...
           abs(mat.nu-C.nu)<=1e-14 && mat.ps==C.ps, ...
        'stage2phys:MaterialMismatch','Candidate/frozen material mismatch.');
    assert(strcmp(C.load.type,'remote_tension_y') && ...
           abs(C.load.sig0-1)<=1e-14 && ...
           strcmp(C.bc.anchor_mode,'minimal'), ...
        'stage2phys:PhysicalConfig', ...
        'Expected frozen unit remote-y loading with minimal anchoring.');

    % Pre-solve native sampling fingerprint.
    zeroU=zeros(2*size(P6,1),1);
    [rZero,~,face0]=native_COD_audit(mesh,zeroU,mat,crack,8);
    rrZero=rZero/a0;
    windows=[.04 .20;.04 .30;.08 .30;.12 .30];
    sampleN=zeros(4,1);
    for k=1:4
        sampleN(k)=nnz(rrZero>=windows(k,1)&rrZero<=windows(k,2));
    end
    if face0.nUpper~=138 || face0.nLower~=138 || ...
            face0.gridMismatch>1e-12 || any(sampleN~=[38;55;44;34])
        error('stage2phys:NativeSampling', ...
            'Expected exact qualified native sampling 38/55/44/34.');
    end

    % Plate geometry and minimal anchors.
    xmin=min(P(:,1));xmax=max(P(:,1));
    ymin=min(P(:,2));ymax=max(P(:,2));
    A=C.A;B=C.B;
    if abs(xmin)>1e-12 || abs(xmax-A)>1e-12 || ...
            abs(ymin+B)>1e-12 || abs(ymax-B)>1e-12
        error('stage2phys:PlateBoundary','Plate boundary fingerprint changed.');
    end

    [~,iLB]=min(sum((P-[0,-B]).^2,2));
    [~,iRB]=min(sum((P-[A,-B]).^2,2));
    if norm(P(iLB,:)-[0,-B])>1e-12 || norm(P(iRB,:)-[A,-B])>1e-12
        error('stage2phys:Corners','Could not identify exact plate bottom corners.');
    end

    fixvar=unique([2*iLB-1;2*iLB;2*iRB]);
    ndof=2*size(P6,1);
    freeMask=true(ndof,1);
    freeMask(fixvar)=false;
    free=find(freeMask);

    [nip2,xip2,w2,Nextr]=integr();
    quad=struct('nip2',nip2,'xip2',xip2,'w2',w2,'Nextr',Nextr);

    pcgTol=1e-10;
    pcgMaxIt=5000;
    preconditionerName='symmetric_gauss_seidel';
    loadFactor=1.0;

    fprintf('\n============================================================\n');
    fprintf('STAGE II: ONE PHYSICAL theta_1=0 SOLVE\n');
    fprintf('============================================================\n');
    fprintf('  Frozen source    : %s\n',frozenSource);
    fprintf('  Candidate source : %s\n',candidateSource);
    fprintf('  a0               : %.9f mm\n',1e3*a0);
    fprintf('  theta_1          : 0 deg\n');
    fprintf('  T3               : %d nodes / %d elements\n',size(P,1),size(T,1));
    fprintf('  T6               : %d nodes / %d elements\n',size(P6,1),size(T6,1));
    fprintf('  DOF/free DOF     : %d / %d\n',ndof,numel(free));
    fprintf('  load             : unit remote-y traction\n');
    fprintf('  solver           : PCG + parameter-free SGS\n');
    fprintf('  PCG tol/maxit    : %.1e / %d\n',pcgTol,pcgMaxIt);
    fprintf('  EDI annulus      : [%.6f, %.6f] mm\n', ...
        1e3*.10*a0,1e3*.65*a0);

    % ------------------------------------------------------------------
    % Phase 1. Reuse exact checkpoint or execute ONE authorized solve.
    % ------------------------------------------------------------------
    newSolve=false;

    if exist(cp,'file')==2
        s0=load(cp,'meta','mesh','U','mat','crack','a0','solverInfo');
        local_validate_checkpoint(s0,mesh,crack,a0,ndof,opt.AlternativeQualifiedCandidate);
        fprintf('\nPHASE 1: valid physical checkpoint found; NO new solve.\n');
    else
        if ~opt.AllowSolve
            error('stage2phys:ExplicitSolveApprovalRequired', ...
                ['Physical driver is prepared but guarded. Rerun with ', ...
                 '''AllowSolve'',true to authorize exactly one theta_1=0 solve.']);
        end

        fprintf('\nPHASE 1: ASSEMBLE UNCLAMPED SYMMETRIC SYSTEM\n');
        fprintf('  stif_assem(...,fixvar=[]) -- no row clamping.\n');

        K=stif_assem(mesh,mat,quad,[]);
        if size(K,1)~=ndof || size(K,2)~=ndof
            error('stage2phys:StiffnessSize','Unexpected stiffness dimensions.');
        end

        symErr=norm(K-K.',1)/max(1,norm(K,1));
        if ~isfinite(symErr) || symErr>5e-13
            error('stage2phys:StiffnessSymmetry', ...
                'Unclamped stiffness symmetry error %.3e exceeds gate.',symErr);
        end

        % Exact unit remote-y traction vector.
        Fload=zeros(ndof,1);
        eps1=1e-8*max(A,2*B);
        elod=edge_loads_T6(mesh.coord,B,eps1);
        qy=loadFactor*C.load.sig0;
        for k=1:size(elod,1)
            node=elod(k,1);
            w=elod(k,2);
            Fload(2*node)=Fload(2*node)+qy*w;
        end

        topMask=abs(mesh.coord(elod(:,1),2)-B)<eps1;
        botMask=abs(mesh.coord(elod(:,1),2)+B)<eps1;
        topResultant=qy*sum(elod(topMask,2));
        bottomResultant=qy*sum(elod(botMask,2));
        rawNetResultant=topResultant+bottomResultant;

        if abs(topResultant-A)>1e-12 || ...
           abs(bottomResultant+A)>1e-12 || ...
           abs(rawNetResultant)>1e-12
            error('stage2phys:LoadResultant', ...
                'Remote traction integration failed: top=%g bottom=%g raw net=%g.', ...
                topResultant,bottomResultant,rawNetResultant);
        end

        % Homogeneous essential values: discard force entries at constrained
        % DOFs exactly as in the audited Step67A free-DOF formulation.
        Fload(fixvar)=0;

        Kff=K(free,free);
        Ff=Fload(free);
        Kff=(Kff+Kff.')/2;

        if any(~isfinite(diag(Kff))) || any(diag(Kff)<=0)
            error('stage2phys:NonpositiveDiagonal', ...
                'Free-DOF stiffness has a nonpositive/nonfinite diagonal.');
        end

        fprintf('  stiffness symmetry error = %.3e\n',symErr);
        fprintf('  top/bottom resultants    = %+g / %+g\n', ...
            topResultant,bottomResultant);
        fprintf('  Building symamd + SGS preconditioner.\n');

        p=symamd(Kff);
        Ap=Kff(p,p);
        bp=Ff(p);
        clear Kff Ff

        d=diag(Ap);
        if any(~isfinite(d)) || any(d<=0)
            error('stage2phys:NonpositiveSGSDiagonal', ...
                'SGS requires a finite positive diagonal.');
        end

        tPre=tic;
        M1=tril(Ap);              % D + L
        Dinv=spdiags(1./d,0,numel(d),numel(d));
        M2=Dinv*M1.';             % D^{-1}(D + U)
        preconditionerSeconds=toc(tPre);

        nnzM1=nnz(M1);
        nnzM2=nnz(M2);
        preconditionerGiB=(16*(nnzM1+nnzM2)+ ...
            8*(size(M1,2)+size(M2,2)+2))/1024^3;

        fprintf('  SGS storage estimate      = %.4f GiB\n',preconditionerGiB);
        fprintf('  Starting exactly ONE authorized physical PCG solve.\n');

        tSolve=tic;
        [xp,flag,relres,iter,resvec]=pcg( ...
            Ap,bp,pcgTol,pcgMaxIt,M1,M2);
        solveSeconds=toc(tSolve);

        if flag~=0 || ~isfinite(relres) || relres>pcgTol || ...
                any(~isfinite(xp))
            error('stage2phys:PCGFailed', ...
                'PCG failed: flag=%d relres=%.3e iter=%d.', ...
                flag,relres,iter);
        end

        uf=zeros(numel(free),1);
        uf(p)=xp;
        U=zeros(ndof,1);
        U(free)=uf;
        U(fixvar)=0;

        trueRelResidual=norm(K(free,:)*U-Fload(free))/ ...
            max(norm(Fload(free)),eps);
        constraintInf=max(abs(U(fixvar)));

        if trueRelResidual>5e-10 || constraintInf>1e-14
            error('stage2phys:ResidualGate', ...
                'True residual/constraint gate failed: %.3e / %.3e.', ...
                trueRelResidual,constraintInf);
        end

        solverInfo=struct( ...
            'method','pcg_free_dof_spd_sgs', ...
            'pcgTol',pcgTol,'pcgMaxIt',pcgMaxIt, ...
            'preconditioner',preconditionerName, ...
            'flag',flag,'relres',relres,'iter',iter, ...
            'trueRelResidual',trueRelResidual, ...
            'constraintInf',constraintInf, ...
            'resvecFinal',resvec(end), ...
            'resvecLength',numel(resvec), ...
            'symmetryError',symErr, ...
            'nnzK',nnz(K),'nnzAp',nnz(Ap), ...
            'nnzM1',nnzM1,'nnzM2',nnzM2, ...
            'preconditionerGiBApprox',preconditionerGiB, ...
            'preconditionerSeconds',preconditionerSeconds, ...
            'solveSeconds',solveSeconds, ...
            'topResultant',topResultant, ...
            'bottomResultant',bottomResultant, ...
            'rawNetVerticalResultant',rawNetResultant, ...
            'onePhysicalLinearSolve',true, ...
            'noDirectBackslash',true);

        meta=struct( ...
            'stage','stage2_theta0_physical_solved', ...
            'source','one guarded full-domain theta_1=0 physical solve', ...
            'nT3Nodes',size(P,1),'nT3',size(T,1), ...
            'nT6Nodes',size(P6,1),'ndof',ndof, ...
            'a0',a0,'theta1',0,'theta1Deg',0, ...
            'mouth',mouth,'tip',tip, ...
            'loadType',C.load.type,'sig0',C.load.sig0, ...
            'loadFactor',loadFactor, ...
            'anchorMode',C.bc.anchor_mode, ...
            'cornerNodes',[iLB iRB], ...
            'solver','pcg_free_dof_spd_sgs', ...
            'pcgTol',pcgTol,'preconditioner',preconditionerName, ...
            'postprocessingPerformedBeforeCheckpoint',false, ...
            'singlePhysicalAngleOnly',true, ...
            'noAngleSweep',true, ...
            'alternativeQualifiedCandidate',logical(opt.AlternativeQualifiedCandidate));

        tmp=[cp '.incomplete.mat'];
        if exist(tmp,'file')==2
            error('stage2phys:InterruptedSave', ...
                'Inspect/remove prior incomplete checkpoint manually: %s',tmp);
        end

        [folder,~,~]=fileparts(cp);
        if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end

        save(tmp,'mesh','U','mat','crack','a0','meta', ...
            'solverInfo','C','-v7.3');
        [ok,msg]=movefile(tmp,cp);
        if ~ok,error('stage2phys:CheckpointSave','%s',msg);end

        fprintf('  Physical field checkpointed BEFORE postprocessing:\n    %s\n',cp);
        newSolve=true;

        clear K Fload Ap bp M1 M2 Dinv xp uf U
    end

    % ------------------------------------------------------------------
    % Phase 2. Postprocess ONLY from the checkpoint.
    % ------------------------------------------------------------------
    fprintf('\nPHASE 2: POSTPROCESS SAVED PHYSICAL FIELD\n');
    s=load(cp,'mesh','U','mat','crack','a0','meta','solverInfo');
    local_validate_checkpoint(s,mesh,crack,a0,ndof);

    [r,app,face]=native_COD_audit( ...
        s.mesh,s.U,s.mat,s.crack,8);
    rr=r/s.a0;

    fitDegrees=[1 2];
    fitRows=nan(8,9);
    qRaw=app(:,2)./app(:,1);
    n=0;
    for iw=1:size(windows,1)
        ids=find(rr>=windows(iw,1)&rr<=windows(iw,2));
        for deg=fitDegrees
            pI=polyfit(rr(ids),app(ids,1),deg);
            pII=polyfit(rr(ids),app(ids,2),deg);
            predII=polyval(pII,rr(ids));
            n=n+1;
            fitRows(n,:)=[windows(iw,:),deg,numel(ids), ...
                pI(end),pII(end),pII(end)/pI(end), ...
                sqrt(mean((app(ids,2)-predII).^2)), ...
                median(qRaw(ids))];
        end
    end

    fitTable=array2table(fitRows, ...
        'VariableNames',{'lower_r_over_a0','upper_r_over_a0','degree', ...
        'n_native','KI_COD','KII_COD','ratio_COD','RMSE_KII', ...
        'median_pointwise_ratio'});

    % Exactly one matched physical EDI.
    ri=.10*s.a0;
    ro=.65*s.a0;
    fprintf('  Running ONE matched physical EDI: r=[%.6f, %.6f] mm.\n', ...
        1e3*ri,1e3*ro);

    [KI,KII,Aux]=SIF_LEFM_interaction_EDI( ...
        s.mesh,s.U,s.crack.Pmid,s.mat, ...
        struct('r_inner',ri,'r_outer',ro), ...
        'UsePlaneStrain',s.mat.ps==1, ...
        'Verbose',false, ...
        'WeightFunction','fe_nodal', ...
        'QuadratureRule',16, ...
        'StoreGPDiagnostics',false);

    qEDI=KII/KI;

    % Linear rescaling to frozen initiation load, if recorded. This is NOT
    % another solve and does not affect the local-symmetry ratio.
    lambdaIni=NaN;
    if ismember('lambda_ini',T0.Properties.VariableNames)
        lambdaIni=T0.lambda_ini;
    end
    KIatIni=lambdaIni*KI;
    KIIatIni=lambdaIni*KII;

    gates=struct();
    gates.candidateQualified=true;
    gates.pcgConverged=s.solverInfo.flag==0;
    gates.pcgReportedResidual=s.solverInfo.relres<=pcgTol;
    gates.trueResidual=s.solverInfo.trueRelResidual<=5e-10;
    gates.constraintsExact=s.solverInfo.constraintInf<=1e-14;
    gates.unitLoadResultants= ...
        abs(s.solverInfo.topResultant-A)<=1e-12 && ...
        abs(s.solverInfo.bottomResultant+A)<=1e-12 && ...
        abs(s.solverInfo.rawNetVerticalResultant)<=1e-12;
    gates.nativeFacePairing= ...
        face.nUpper==138 && face.nLower==138 && face.gridMismatch<=1e-12;
    gates.nativeSamplingExact=isequal( ...
        fitTable.n_native,[38;38;55;55;44;44;34;34]);
    gates.physicalEDIFinite=isfinite(KI)&&isfinite(KII)&&isfinite(qEDI);
    gates.modeIPositive=KI>0;
    gates.matchedEDISupport=Aux.nElem_used==11316;
    gates.singlePhysicalAngle=logical(s.meta.singlePhysicalAngleOnly);
    gates.noAngleSweep=logical(s.meta.noAngleSweep);
    gates.noDirectBackslash=logical(s.solverInfo.noDirectBackslash);

    pass=all(structfun(@logical,gates));

    EDI=table(KI,KII,qEDI,Aux.nElem_used,Aux.nGP_used, ...
        lambdaIni,KIatIni,KIIatIni, ...
        'VariableNames',{'KI_unit','KII_unit','KII_over_KI', ...
        'EDI_elements','EDI_Gauss_points', ...
        'lambda_ini','KI_at_lambda_ini','KII_at_lambda_ini'});

    Solver=struct2table(s.solverInfo,'AsArray',true);

    fprintf('\nPHYSICAL COD FITS\n');
    disp(fitTable);
    fprintf('\nPHYSICAL EDI\n');
    disp(EDI);
    fprintf('\nSOLVER INFO\n');
    disp(Solver);
    fprintf('\nPHYSICAL-SOLVE GATES\n');
    disp(gates);

    fprintf('\nPRIMARY LOCAL-SYMMETRY OBSERVABLE\n');
    fprintf('  KI(0)       = %.12g MPa*sqrt(m) at unit traction\n',KI);
    fprintf('  KII(0)      = %+.12g MPa*sqrt(m) at unit traction\n',KII);
    fprintf('  KII/KI      = %+.12g\n',qEDI);
    if isfinite(lambdaIni)
        fprintf('  lambda_ini  = %.12g\n',lambdaIni);
        fprintf('  scaled KI   = %.12g MPa*sqrt(m)\n',KIatIni);
        fprintf('  scaled KII  = %+.12g MPa*sqrt(m)\n',KIIatIni);
    end

    if pass
        fprintf('\nSTAGE-II theta_1=0 PHYSICAL SOLVE PASS.\n');
        fprintf('  One physical angle only; no angle sweep was performed.\n');
        fprintf('  KII(0) is now available to choose the next local probe.\n');
    else
        fprintf('\nSTAGE-II theta_1=0 PHYSICAL SOLVE STOP.\n');
        fprintf('  Do not use KII(0) for angle selection until failed gates are resolved.\n');
    end

    Summary=table( ...
        size(P,1),size(T,1),size(P6,1),ndof,numel(free), ...
        s.solverInfo.iter,s.solverInfo.relres,s.solverInfo.trueRelResidual, ...
        KI,KII,qEDI,Aux.nElem_used,Aux.nGP_used,newSolve,pass, ...
        'VariableNames',{ ...
        'T3_nodes','T3_elements','T6_nodes','DOF','free_DOF', ...
        'PCG_iterations','PCG_relres','true_rel_residual', ...
        'KI_unit','KII_unit','KII_over_KI', ...
        'EDI_elements','EDI_Gauss_points','newSolve','pass'});

    R=struct();
    R.summary=Summary;
    R.fitTable=fitTable;
    R.EDI=EDI;
    R.solverInfo=s.solverInfo;
    R.gates=gates;
    R.pass=pass;
    R.newSolve=newSolve;
    R.checkpointPath=cp;
    R.candidateSource=candidateSource;
    R.frozenSource=frozenSource;
    R.alternativeQualifiedCandidate=logical(opt.AlternativeQualifiedCandidate);
    if isfield(candidate,'exteriorMeshControls')
        R.exteriorMeshControls=candidate.exteriorMeshControls;
    end
    R.interpretation=[ ...
        'One qualified theta_1=0 physical LEFM solve. ', ...
        'The primary first-segment direction observable is signed KII(0). ', ...
        'No angular derivative or root bracket has yet been measured.'];

    save(saveFile,'R','-v7');
    fprintf('  Compact physical result saved: %s\n',saveFile);
end


% =========================================================================
function [R0,label]=local_load_frozen(root,Rin)
    if ~isempty(Rin)
        R0=Rin;
        label='<in-memory FrozenState>';
        return
    end

    f=fullfile(root,'verification','crack_path','stage1_starting_state.mat');
    if exist(f,'file')~=2
        error('stage2phys:MissingFrozenState', ...
            ['Frozen Stage-I MAT not found. Pass ''FrozenState'',R0 ', ...
             'if the accepted struct is still in memory.']);
    end
    d=load(f);
    assert(isfield(d,'R0')&&isstruct(d.R0), ...
        'stage2phys:BadFrozenState','Frozen MAT must contain R0.');
    R0=d.R0;
    label=f;
end

function [candidate,label]=local_load_candidate(cin,file)
    if ~isempty(cin)
        candidate=cin;
        label='<in-memory full-domain candidate>';
        return
    end
    if exist(file,'file')~=2
        error('stage2phys:MissingCandidate', ...
            ['Qualified full-domain candidate not found. Pass ', ...
             '''Candidate'',F.candidate if F is still in memory.']);
    end
    d=load(file);
    assert(isfield(d,'candidate')&&isstruct(d.candidate), ...
        'stage2phys:BadCandidateMAT','Candidate MAT must contain candidate.');
    candidate=d.candidate;
    label=file;
end

function local_require_candidate(c)
    req={'p','t','crack','mat','a0','theta1', ...
        'pairedElementIDs','primarySupportElementIDs', ...
        'gates','syntheticGates','scientificallyReadyForOnePhysicalSolve'};
    for k=1:numel(req)
        if ~isfield(c,req{k}) || isempty(c.(req{k}))
            error('stage2phys:CandidateField', ...
                'Candidate missing required field %s.',req{k});
        end
    end
end

function local_require_table_variables(T,names)
    miss=names(~ismember(names,T.Properties.VariableNames));
    if ~isempty(miss)
        error('stage2phys:FrozenFields', ...
            'Missing frozen summary fields: %s',strjoin(miss,', '));
    end
end

function mat=local_material_with_D(mat0,C)
    mat=mat0;
    mat.E=C.E;
    mat.nu=C.nu;
    mat.ps=C.ps;
    E=mat.E;nu=mat.nu;
    if mat.ps==1
        coef=E/((1+nu)*(1-2*nu));
        D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
    else
        coef=E/(1-nu^2);
        D=coef*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
    end
    mat.D=D;
    mat.Dmat=D;
end

function local_validate_checkpoint(s,mesh,crack,a0,ndof,alternativeQualifiedCandidate)
    req={'meta','mesh','U','mat','crack','a0','solverInfo'};
    for k=1:numel(req)
        if ~isfield(s,req{k})
            error('stage2phys:CheckpointField', ...
                'Physical checkpoint missing %s.',req{k});
        end
    end

    commonMismatch = ...
            ~strcmp(s.meta.stage,'stage2_theta0_physical_solved') || ...
            s.meta.ndof~=ndof || ...
            abs(s.a0-a0)>1e-14 || s.meta.theta1~=0 || ...
            ~strcmp(s.meta.solver,'pcg_free_dof_spd_sgs') || ...
            ~isequal(s.mesh.connect3,mesh.connect3) || ...
            max(abs(s.mesh.coord3(:)-mesh.coord3(:)))>1e-12 || ...
            numel(s.U)~=ndof || any(~isfinite(s.U)) || ...
            norm(s.crack.Pmid-crack.Pmid,'fro')>2e-12;

    if alternativeQualifiedCandidate
        countMismatch = s.meta.nT3Nodes~=size(mesh.coord3,1) || ...
            s.meta.nT3~=size(mesh.connect3,1) || ...
            s.meta.nT6Nodes~=size(mesh.coord,1);
    else
        countMismatch = s.meta.nT3Nodes~=24389 || s.meta.nT3~=47828 || ...
            s.meta.nT6Nodes~=96606;
    end

    if commonMismatch || countMismatch
        error('stage2phys:CheckpointMismatch', ...
            'Existing physical checkpoint is not the exact current theta0 candidate.');
    end
end

function tf=local_matlab_callable(name)
    tf=exist(name,'file')~=0 || exist(name,'builtin')~=0;
end


function local_assert_branch(root)
    [status,b]=system(sprintf('git -C "%s" branch --show-current',root));
    assert(status==0 && strcmp(strtrim(b),'stage2-theta0-physical-solve'), ...
        'stage2phys:Branch', ...
        ['Run the guarded theta0 physical qualification only on ', ...
         'stage2-theta0-physical-solve.']);
end
