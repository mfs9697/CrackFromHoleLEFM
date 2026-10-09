function R = solve_incremental_crack_tip(candidate,varargin)
%SOLVE_INCREMENTAL_CRACK_TIP
% Perform exactly one guarded physical LEFM solve at the current tip of an
% already-qualified arbitrary polyline crack candidate.
%
% The crack path is fixed geometry. This function does not append a segment,
% vary an angle, or search over directions. It measures KI,KII at the current
% tip, then predicts the NEXT incremental turn by MTS.
%
% DEFAULT IS SAFE: AllowSolve=false.
%
% Required candidate source: qualify_incremental_crack_candidate.

    assert(isstruct(candidate)&&isscalar(candidate), ...
        'pathsolve:Candidate','Qualified candidate struct is required.');

    ip=inputParser;
    addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FastEDI',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FrozenState',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'CheckpointFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'SaveFile','',@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    incremental_profile_clock('begin','physical',candidate.nSegments);
    profileCleanup=onCleanup(@()incremental_profile_clock('end','physical'));
    root=fileparts(mfilename('fullpath'));
    outDir=fullfile(root,'verification','crack_path');
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    [R0,frozenSource]=local_load_frozen(root,opt.FrozenState);
    local_require_candidate(candidate);
    candidateSource='<in-memory qualified incremental candidate>';

    nSegments=candidate.nSegments;
    cp=char(opt.CheckpointFile);
    if isempty(cp)
        cp=fullfile(outDir,sprintf('incremental_step_%03d_physical_solved.mat',nSegments));
    end
    saveFile=char(opt.SaveFile);
    if isempty(saveFile)
        saveFile=fullfile(outDir,sprintf('incremental_step_%03d_physical_small.mat',nSegments));
    end

    for name={'pcg','symamd'}
        if ~local_matlab_callable(name{1})
            error('pathsolve:MissingIterativeTool', ...
                'Physical solve requires callable MATLAB %s.',name{1});
        end
    end

    % ------------------------------------------------------------------
    % Frozen physical configuration and exact candidate provenance.
    % ------------------------------------------------------------------
    assert(isfield(R0,'summary')&&istable(R0.summary)&&height(R0.summary)==1, ...
        'pathsolve:FrozenSummary','Frozen R0.summary is required.');
    assert(isfield(R0,'C')&&isstruct(R0.C), ...
        'pathsolve:FrozenConfig','Frozen R0.C is required.');

    T0=R0.summary;
    req={'a0_reserved_m','x_star_m','y_star_m', ...
        'nmat_x','nmat_y','that_x','that_y','stage1_pass'};
    local_require_table_variables(T0,req);
    row=T0(1,:);
    assert(logical(row.stage1_pass),'pathsolve:Stage1NotPassed', ...
        'Frozen Stage-I state did not pass.');

    C=R0.C;
    increment=row.a0_reserved_m;
    mouth=[row.x_star_m,row.y_star_m];
    nMat=[row.nmat_x,row.nmat_y];nMat=nMat/norm(nMat);
    tHat=[row.that_x,row.that_y];tHat=tHat/norm(tHat);

    path=candidate.path;
    assert(size(path,1)==nSegments+1&&nSegments>=2, ...
        'pathsolve:PathSize','Candidate path/segment count is inconsistent.');
    seg=diff(path,1,1);
    segLength=vecnorm(seg,2,2);
    assert(max(abs(segLength-increment))<=2e-12, ...
        'pathsolve:IncrementLength','Candidate does not retain the frozen increment.');

    directions=seg./segLength;
    thetaSegments=atan2(directions*tHat(:),directions*nMat(:));
    thetaSegmentsDeg=rad2deg(thetaSegments);
    thetaCurrent=thetaSegments(end);
    thetaCurrentDeg=thetaSegmentsDeg(end);

    p0=path(1,:);
    priorTip=path(end-1,:);
    tip=path(end,:);
    currentIncrement=segLength(end);
    totalCrack=sum(segLength);

    if ~candidate.scientificallyReadyForIncrementalPhysicalSolve
        error('pathsolve:CandidateNotQualified', ...
            'Candidate is not marked ready for an incremental physical solve.');
    end
    if ~all(structfun(@logical,candidate.gates)) || ...
            ~all(structfun(@logical,candidate.syntheticGates))
        error('pathsolve:CandidateGateFailure', ...
            'Full-domain candidate does not retain all structural/synthetic gates.');
    end

    P=candidate.p;
    T=candidate.t;
    crack=candidate.crack;

    assert(candidate.nSegments==nSegments && ...
           abs(candidate.currentIncrementLength-currentIncrement)<=1e-14 && ...
           abs(candidate.totalCrackLength-totalCrack)<=2e-12 && ...
           norm(candidate.path-path,'fro')<=2e-12 && ...
           max(abs(candidate.thetaSegmentsDeg-thetaSegmentsDeg))<=1e-10, ...
        'pathsolve:CandidateCrackParameter', ...
        'Candidate path metadata does not match its qualified geometry.');
    assert(norm(path(1,:)-mouth)<=2e-12 && ...
           norm(candidate.priorTip-priorTip)<=2e-12 && ...
           norm(crack.Pmid-path,'fro')<=2e-12, ...
        'pathsolve:CandidateCrackGeometry', ...
        'Candidate mouth/prior-tip/current-tip geometry changed.');
    coreFP=local_core_fingerprint(candidate,currentIncrement);
    assert(numel(candidate.pairedElementIDs)==coreFP.expectedCoreT3, ...
        'pathsolve:CoreFingerprint','Qualified core element count changed.');
    assert(numel(candidate.primarySupportElementIDs)==coreFP.expectedEDIElements, ...
        'pathsolve:SupportFingerprint','Qualified production EDI support count changed.');
    if isfield(candidate,'skipConstantSupportElementIDs')
        assert(numel(candidate.skipConstantSupportElementIDs)== ...
                coreFP.expectedOptimizedSupport, ...
            'pathsolve:OptimizedSupportFingerprint', ...
            'Qualified skip-constant support count changed.');
        assert(all(ismember(candidate.skipConstantSupportElementIDs, ...
                candidate.primarySupportElementIDs)), ...
            'pathsolve:OptimizedSupportContainment', ...
            'Skip-constant diagnostic support is not contained in production EDI support.');
    end
    incremental_profile_clock('phase','physical','t3_t6');
    [P6,T6]=T3toT6_fast(P,T);
    mesh=struct('coord3',P,'connect3',T,'coord',P6,'connect',T6);

    incremental_profile_clock('phase','physical','preflight');
    assert(isfield(candidate.mat,'E')&&isfield(candidate.mat,'nu')&& ...
           isfield(candidate.mat,'ps') && ...
           abs(candidate.mat.E-C.E)<=1e-12*max(1,abs(C.E)) && ...
           abs(candidate.mat.nu-C.nu)<=1e-14 && candidate.mat.ps==C.ps, ...
        'pathsolve:CandidateMaterialMismatch', ...
        'Qualified candidate material differs from frozen Stage-I material.');

    mat=local_material_with_D(candidate.mat,C);

    % Physics must be the exact frozen Stage-I material/loading/anchoring.
    assert(abs(mat.E-C.E)<=1e-12*max(1,abs(C.E)) && ...
           abs(mat.nu-C.nu)<=1e-14 && mat.ps==C.ps, ...
        'pathsolve:MaterialMismatch','Candidate/frozen material mismatch.');
    assert(strcmp(C.load.type,'remote_tension_y') && ...
           abs(C.load.sig0-1)<=1e-14 && ...
           strcmp(C.bc.anchor_mode,'minimal'), ...
        'pathsolve:PhysicalConfig', ...
        'Expected frozen unit remote-y loading with minimal anchoring.');

    % Pre-solve native sampling fingerprint.
    zeroU=zeros(2*size(P6,1),1);
    [rZero,~,face0]=native_COD_polyline_audit(mesh,zeroU,mat,crack,8);
    rrZero=rZero/currentIncrement;
    windows=[.04 .20;.04 .30;.08 .30;.12 .30];
    sampleN=zeros(4,1);
    for k=1:4
        sampleN(k)=nnz(rrZero>=windows(k,1)&rrZero<=windows(k,2));
    end
    if face0.nUpper~=face0.nLower || face0.nUpper<34 || ...
            ~logical(face0.usesLastSegmentFrame) || ...
            face0.gridMismatch>1e-12 || any(sampleN~=coreFP.expectedNativeSamples)
        error('pathsolve:NativeSampling', ...
            'Native sampling does not match the qualified core-family fingerprint.');
    end

    % Plate geometry and minimal anchors.
    xmin=min(P(:,1));xmax=max(P(:,1));
    ymin=min(P(:,2));ymax=max(P(:,2));
    A=C.A;B=C.B;
    if abs(xmin)>1e-12 || abs(xmax-A)>1e-12 || ...
            abs(ymin+B)>1e-12 || abs(ymax-B)>1e-12
        error('pathsolve:PlateBoundary','Plate boundary fingerprint changed.');
    end

    [~,iLB]=min(sum((P-[0,-B]).^2,2));
    [~,iRB]=min(sum((P-[A,-B]).^2,2));
    if norm(P(iLB,:)-[0,-B])>1e-12 || norm(P(iRB,:)-[A,-B])>1e-12
        error('pathsolve:Corners','Could not identify exact plate bottom corners.');
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
    fprintf('INCREMENTAL PATH: PHYSICAL SOLVE AT CURRENT TIP\n');
    fprintf('============================================================\n');
    fprintf('  Frozen source    : %s\n',frozenSource);
    fprintf('  Candidate source : %s\n',candidateSource);
    fprintf('  first leg        : %.9f mm, theta_1=%+.12g deg\n',1e3*increment,0);
    fprintf('  segments         : %d\n',nSegments);
    fprintf('  current angle    : %+.12g deg (PRESCRIBED GEOMETRY)\n',thetaCurrentDeg);
    fprintf('  total crack      : %.9f mm\n',1e3*totalCrack);
    fprintf('  current increment: %.9f mm\n',1e3*currentIncrement);
    fprintf('  T3               : %d nodes / %d elements\n',size(P,1),size(T,1));
    fprintf('  T6               : %d nodes / %d elements\n',size(P6,1),size(T6,1));
    fprintf('  DOF/free DOF     : %d / %d\n',ndof,numel(free));
    fprintf('  load             : unit remote-y traction\n');
    fprintf('  solver           : PCG + parameter-free SGS\n');
    fprintf('  PCG tol/maxit    : %.1e / %d\n',pcgTol,pcgMaxIt);
    fprintf('  EDI annulus      : [%.6f, %.6f] mm\n', ...
        1e3*.10*currentIncrement,1e3*.65*currentIncrement);

    % ------------------------------------------------------------------
    % Phase 1. Reuse exact checkpoint or execute ONE authorized solve.
    % ------------------------------------------------------------------
    newSolve=false;

    incremental_profile_clock('phase','physical','checkpoint_validation');
    if exist(cp,'file')==2
        s0=load(cp,'meta','mesh','U','mat','crack','currentIncrement','solverInfo','C');
        local_validate_checkpoint(s0,mesh,crack,currentIncrement,ndof,path,thetaSegmentsDeg,mat,C);
        fprintf('\nPHASE 1: valid physical checkpoint found; NO new solve.\n');
    else
        if ~opt.AllowSolve
            error('pathsolve:ExplicitSolveApprovalRequired', ...
                ['Physical driver is prepared but guarded. Rerun with ', ...
                 '''AllowSolve'',true to authorize exactly one two-leg new-tip solve.']);
        end

        fprintf('\nPHASE 1: ASSEMBLE UNCLAMPED SYMMETRIC SYSTEM\n');
        fprintf('  stif_assem(...,fixvar=[]) -- no row clamping.\n');

        incremental_profile_clock('phase','physical','stiffness_assembly');
        K=stif_assem(mesh,mat,quad,[]);
        incremental_profile_clock('phase','physical','stiffness_checks');
        if size(K,1)~=ndof || size(K,2)~=ndof
            error('pathsolve:StiffnessSize','Unexpected stiffness dimensions.');
        end

        symErr=norm(K-K.',1)/max(1,norm(K,1));
        if ~isfinite(symErr) || symErr>5e-13
            error('pathsolve:StiffnessSymmetry', ...
                'Unclamped stiffness symmetry error %.3e exceeds gate.',symErr);
        end

        % Exact unit remote-y traction vector.
        incremental_profile_clock('phase','physical','loads');
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
            error('pathsolve:LoadResultant', ...
                'Remote traction integration failed: top=%g bottom=%g raw net=%g.', ...
                topResultant,bottomResultant,rawNetResultant);
        end

        % Homogeneous essential values: discard force entries at constrained
        % DOFs exactly as in the audited Step67A free-DOF formulation.
        Fload(fixvar)=0;

        incremental_profile_clock('phase','physical','free_dof_extraction');
        Kff=K(free,free);
        Ff=Fload(free);
        Kff=(Kff+Kff.')/2;

        if any(~isfinite(diag(Kff))) || any(diag(Kff)<=0)
            error('pathsolve:NonpositiveDiagonal', ...
                'Free-DOF stiffness has a nonpositive/nonfinite diagonal.');
        end

        fprintf('  stiffness symmetry error = %.3e\n',symErr);
        fprintf('  top/bottom resultants    = %+g / %+g\n', ...
            topResultant,bottomResultant);
        fprintf('  Building symamd + SGS preconditioner.\n');

        incremental_profile_clock('phase','physical','symamd_and_permutation');
        p=symamd(Kff);
        Ap=Kff(p,p);
        bp=Ff(p);
        clear Kff Ff

        d=diag(Ap);
        if any(~isfinite(d)) || any(d<=0)
            error('pathsolve:NonpositiveSGSDiagonal', ...
                'SGS requires a finite positive diagonal.');
        end

        incremental_profile_clock('phase','physical','sgs_construction');
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
        fprintf('  Starting one authorized physical linear solve (one same-settings PCG restart permitted only on flag=3 stagnation).\n');

        incremental_profile_clock('phase','physical','pcg');
        tSolve=tic;
        [xp,flag,relres,iter,resvec]=pcg( ...
            Ap,bp,pcgTol,pcgMaxIt,M1,M2);
        firstSolveSeconds=toc(tSolve);
        solveSeconds=firstSolveSeconds;

        firstFlag=flag;
        firstRelres=relres;
        firstIter=iter;
        firstResvec=resvec;
        firstTrueRelResidual=norm(Ap*xp-bp)/max(norm(bp),eps);
        restarted=false;
        restartIter=0;
        restartFlag=NaN;
        restartRelres=NaN;
        restartTrueRelResidual=NaN;
        restartResvec=[];

        % A single Krylov restart is permitted only for MATLAB flag=3
        % (stagnation) with a finite returned iterate. It uses the SAME
        % matrix, RHS, SGS factors, tolerance and maxit, with the stalled
        % iterate as x0. No physical/numerical acceptance gate is relaxed.
        if flag==3 && isfinite(relres) && all(isfinite(xp))
            fprintf('  PCG stagnated; performing ONE same-settings restart from stalled iterate.\n');
            fprintf('    first relres / true rel = %.16e / %.16e\n', ...
                firstRelres,firstTrueRelResidual);
            tRestart=tic;
            [xp,flag,relres,restartIter,restartResvec]=pcg( ...
                Ap,bp,pcgTol,pcgMaxIt,M1,M2,xp);
            restartSeconds=toc(tRestart);
            solveSeconds=solveSeconds+restartSeconds;
            restartFlag=flag;
            restartRelres=relres;
            restartTrueRelResidual=norm(Ap*xp-bp)/max(norm(bp),eps);
            restarted=true;
            iter=firstIter+restartIter;
            resvec=restartResvec;
            fprintf('    restart flag / iter     = %d / %d\n',restartFlag,restartIter);
            fprintf('    restart relres / true   = %.16e / %.16e\n', ...
                restartRelres,restartTrueRelResidual);
        else
            restartSeconds=0;
        end

        incremental_profile_clock('phase','physical','residual_gates');

        % Always reconstruct the returned iterate and evaluate the true
        % free-system residual before accepting OR rejecting PCG.
        uf=zeros(numel(free),1);
        uf(p)=xp;
        U=zeros(ndof,1);
        U(free)=uf;
        U(fixvar)=0;

        trueRelResidual=norm(K(free,:)*U-Fload(free))/ ...
            max(norm(Fload(free)),eps);
        constraintInf=max(abs(U(fixvar)));

        if flag~=0 || ~isfinite(relres) || relres>pcgTol || ...
                any(~isfinite(xp))
            tailCount=min(8,numel(resvec));
            tail=resvec(end-tailCount+1:end)/max(norm(bp),eps);
            fprintf('  PCG REJECTED before postprocessing.\n');
            fprintf('    final flag / total iter = %d / %d\n',flag,iter);
            fprintf('    reported relres         = %.16e\n',relres);
            fprintf('    recomputed true rel     = %.16e\n',trueRelResidual);
            fprintf('    final residual tail     = %s\n',mat2str(tail(:).',8));
            error('pathsolve:PCGFailed', ...
                ['PCG failed unchanged acceptance gates: flag=%d, ', ...
                 'relres=%.3e, trueRel=%.3e, totalIter=%d.'], ...
                flag,relres,trueRelResidual,iter);
        end

        if trueRelResidual>5e-10 || constraintInf>1e-14
            error('pathsolve:ResidualGate', ...
                'True residual/constraint gate failed: %.3e / %.3e.', ...
                trueRelResidual,constraintInf);
        end

        solverInfo=struct( ...
            'method','pcg_free_dof_spd_sgs', ...
            'pcgTol',pcgTol,'pcgMaxIt',pcgMaxIt, ...
            'preconditioner',preconditionerName, ...
            'flag',flag,'relres',relres,'iter',iter, ...
            'trueRelResidual',trueRelResidual, ...
            'pcgRestarted',restarted, ...
            'pcgCalls',1+double(restarted), ...
            'firstFlag',firstFlag,'firstRelres',firstRelres, ...
            'firstIter',firstIter,'firstTrueRelResidual',firstTrueRelResidual, ...
            'restartFlag',restartFlag,'restartRelres',restartRelres, ...
            'restartIter',restartIter, ...
            'restartTrueRelResidual',restartTrueRelResidual, ...
            'firstSolveSeconds',firstSolveSeconds, ...
            'restartSeconds',restartSeconds, ...
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
            'stage','incremental_crack_tip_physical_solved', ...
            'source','one guarded physical solve of a qualified incremental polyline state', ...
            'nT3Nodes',size(P,1),'nT3',size(T,1), ...
            'nT6Nodes',size(P6,1),'ndof',ndof, ...
            'nSegments',nSegments, ...
            'currentIncrement',currentIncrement,'totalCrackLength',totalCrack, ...
            'thetaSegments',thetaSegments,'thetaSegmentsDeg',thetaSegmentsDeg, ...
            'thetaCurrent',thetaCurrent,'thetaCurrentDeg',thetaCurrentDeg, ...
            'path',path,'mouth',p0,'priorTip',priorTip,'tip',tip, ...
            'loadType',C.load.type,'sig0',C.load.sig0, ...
            'loadFactor',loadFactor, ...
            'anchorMode',C.bc.anchor_mode, ...
            'cornerNodes',[iLB iRB], ...
            'solver','pcg_free_dof_spd_sgs', ...
            'pcgTol',pcgTol,'preconditioner',preconditionerName, ...
            'postprocessingPerformedBeforeCheckpoint',false, ...
            'singlePhysicalStateOnly',true, ...
            'currentPathPrescribed',true, ...
            'noThirdLegGenerated',true, ...
            'noAngleSweep',true);

        incremental_profile_clock('phase','physical','checkpoint_save');
        tmp=[cp '.incomplete.mat'];
        if exist(tmp,'file')==2
            error('pathsolve:InterruptedSave', ...
                'Inspect/remove prior incomplete checkpoint manually: %s',tmp);
        end

        [folder,~,~]=fileparts(cp);
        if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end

        save(tmp,'mesh','U','mat','crack','currentIncrement','meta', ...
            'solverInfo','C','-v7.3');
        [ok,msg]=movefile(tmp,cp);
        if ~ok,error('pathsolve:CheckpointSave','%s',msg);end

        fprintf('  Physical field checkpointed BEFORE postprocessing:\n    %s\n',cp);
        newSolve=true;

        clear K Fload Ap bp M1 M2 Dinv xp uf U
    end

    % ------------------------------------------------------------------
    % Phase 2. Postprocess ONLY from the checkpoint.
    % ------------------------------------------------------------------
    fprintf('\nPHASE 2: POSTPROCESS SAVED PHYSICAL FIELD\n');
    incremental_profile_clock('phase','physical','checkpoint_reload');
    s=load(cp,'mesh','U','mat','crack','currentIncrement','meta','solverInfo','C');
    local_validate_checkpoint(s,mesh,crack,currentIncrement,ndof,path,thetaSegmentsDeg,mat,C);

    incremental_profile_clock('phase','physical','native_cod');
    [r,app,face]=native_COD_polyline_audit( ...
        s.mesh,s.U,s.mat,s.crack,8);
    incremental_profile_clock('phase','physical','cod_fits');
    rr=r/currentIncrement;

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
        'VariableNames',{'lower_r_over_DeltaA','upper_r_over_DeltaA','degree', ...
        'n_native','KI_COD','KII_COD','ratio_COD','RMSE_KII', ...
        'median_pointwise_ratio'});
    deltaCOD=nan(height(fitTable),1);
    for k=1:height(fitTable)
        [~,deltaCOD(k)]=kink_angle_LEFM_MTS( ...
            fitTable.KI_COD(k),fitTable.KII_COD(k));
    end
    fitTable.delta_theta_next_MTS_deg=deltaCOD;

    % Exactly one matched physical EDI.
    ri=.10*currentIncrement;
    ro=.65*currentIncrement;
    fprintf('  Running ONE matched physical EDI: r=[%.6f, %.6f] mm.\n', ...
        1e3*ri,1e3*ro);

    incremental_profile_clock('phase','physical','physical_edi');
    [KI,KII,Aux]=SIF_LEFM_interaction_EDI( ...
        s.mesh,s.U,s.crack.Pmid,s.mat, ...
        struct('r_inner',ri,'r_outer',ro), ...
        'UsePlaneStrain',s.mat.ps==1, ...
        'Verbose',false, ...
        'WeightFunction','fe_nodal', ...
        'QuadratureRule',16, ...
        'StoreGPDiagnostics',false,'SkipUnusedAuxWork',opt.FastEDI);

    incremental_profile_clock('phase','physical','postprocessing_gates');
    qEDI=KII/KI;

    [deltaThetaNext,deltaThetaNextDeg,sigmaTTMTS]=kink_angle_LEFM_MTS(KI,KII);
    thetaNext=thetaCurrent+deltaThetaNext;
    thetaNextDeg=rad2deg(thetaNext);
    globalThetaCurrentDeg=rad2deg(atan2(nMat(2),nMat(1)))+thetaCurrentDeg;
    globalThetaNextDeg=globalThetaCurrentDeg+deltaThetaNextDeg;

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
        face.nUpper==face.nLower && face.nUpper>=34 && face.gridMismatch<=1e-12;
    gates.lastSegmentCODFrame=logical(face.usesLastSegmentFrame);
    expectedFitNative=repelem(coreFP.expectedNativeSamples,2);
    gates.nativeSamplingExact=isequal(fitTable.n_native,expectedFitNative);
    gates.physicalEDIFinite=isfinite(KI)&&isfinite(KII)&&isfinite(qEDI);
    gates.modeIPositive=KI>0;
    gates.matchedEDISupport=Aux.nElem_used==coreFP.expectedEDIElements;
    gates.mtsFinite=isfinite(deltaThetaNext)&&isfinite(deltaThetaNextDeg)&& ...
        isfinite(thetaNextDeg)&&isfinite(globalThetaNextDeg)&&isfinite(sigmaTTMTS);
    gates.codMTSFinite=all(isfinite(fitTable.delta_theta_next_MTS_deg));
    gates.checkpointBeforePostprocessing= ...
        ~logical(s.meta.postprocessingPerformedBeforeCheckpoint);
    gates.singlePhysicalState=logical(s.meta.singlePhysicalStateOnly);
    gates.currentPathPrescribed=logical(s.meta.currentPathPrescribed) && ...
        norm(s.meta.path-path,'fro')<=2e-12;
    gates.currentAngleMatchesCandidate= ...
        abs(candidate.thetaCurrentDeg-thetaCurrentDeg)<=1e-10;
    gates.pathHasAtLeastTwoSegments=nSegments>=2;
    gates.noThirdLegGenerated=logical(s.meta.noThirdLegGenerated);
    gates.noAngleSweep=logical(s.meta.noAngleSweep);
    gates.noDirectBackslash=logical(s.solverInfo.noDirectBackslash);

    pass=all(structfun(@logical,gates));

    EDI=table(KI,KII,qEDI,Aux.nElem_used,Aux.nGP_used, ...
        deltaThetaNextDeg,thetaNextDeg,globalThetaNextDeg,sigmaTTMTS, ...
        lambdaIni,KIatIni,KIIatIni, ...
        'VariableNames',{'KI_unit','KII_unit','KII_over_KI', ...
        'EDI_elements','EDI_Gauss_points', ...
        'delta_theta_next_MTS_deg','theta_next_local_deg','theta_next_global_deg', ...
        'sigmaTT_MTS','lambda_ini','KI_at_lambda_ini','KII_at_lambda_ini'});

    Prediction=table(thetaCurrentDeg,globalThetaCurrentDeg,deltaThetaNextDeg, ...
        thetaNextDeg,globalThetaNextDeg,sigmaTTMTS, ...
        'VariableNames',{'theta_current_local_deg','theta_current_global_deg', ...
        'delta_theta_next_MTS_deg','theta_next_local_deg','theta_next_global_deg','sigmaTT_MTS'});

    Solver=struct2table(s.solverInfo,'AsArray',true);

    fprintf('\nPHYSICAL COD FITS\n');
    disp(fitTable);
    fprintf('\nPHYSICAL EDI\n');
    disp(EDI);
    fprintf('\nMTS PREDICTION FOR NEXT SEGMENT\n');
    disp(Prediction);
    fprintf('\nSOLVER INFO\n');
    disp(Solver);
    fprintf('\nPHYSICAL-SOLVE GATES\n');
    disp(gates);

    fprintf('\nNEW-TIP PHYSICAL OBSERVABLES AND MTS PREDICTION\n');
    fprintf('  theta_k (fixed) = %+.12g deg local\n',thetaCurrentDeg);
    fprintf('  KI(P_k)          = %.12g MPa*sqrt(m) at unit traction\n',KI);
    fprintf('  KII(P_k)         = %+.12g MPa*sqrt(m) at unit traction\n',KII);
    fprintf('  KII/KI          = %+.12g\n',qEDI);
    fprintf('  Delta theta_{k+1}   = %+.12g deg (MTS, relative to current leg)\n',deltaThetaNextDeg);
    fprintf('  theta_{k+1} local   = %+.12g deg\n',thetaNextDeg);
    fprintf('  theta_{k+1} global  = %+.12g deg\n',globalThetaNextDeg);
    if isfinite(lambdaIni)
        fprintf('  lambda_ini      = %.12g\n',lambdaIni);
        fprintf('  scaled KI       = %.12g MPa*sqrt(m)\n',KIatIni);
        fprintf('  scaled KII      = %+.12g MPa*sqrt(m)\n',KIIatIni);
    end

    if pass
        fprintf('\nINCREMENTAL TIP PHYSICAL SOLVE PASS.\n');
        fprintf('  The current crack path remained prescribed; no angle search was performed.\n');
        fprintf('  MTS predicts Delta theta_{k+1}, but no next leg was generated.\n');
    else
        fprintf('\nINCREMENTAL TIP PHYSICAL SOLVE STOP.\n');
        fprintf('  Do not use the MTS prediction until failed gates are resolved.\n');
    end

    Summary=table( ...
        nSegments,currentIncrement,totalCrack,thetaCurrentDeg,coreFP.scale,coreFP.hTip_m, ...
        size(P,1),size(T,1),size(P6,1),ndof,numel(free), ...
        s.solverInfo.iter,s.solverInfo.relres,s.solverInfo.trueRelResidual, ...
        KI,KII,qEDI,deltaThetaNextDeg,thetaNextDeg,globalThetaNextDeg, ...
        Aux.nElem_used,Aux.nGP_used,newSolve,pass, ...
        'VariableNames',{ ...
        'n_segments','increment_m','total_crack_m','theta_current_deg','core_scale','hTip_m', ...
        'T3_nodes','T3_elements','T6_nodes','DOF','free_DOF', ...
        'PCG_iterations','PCG_relres','true_rel_residual', ...
        'KI_unit','KII_unit','KII_over_KI','delta_theta_next_MTS_deg', ...
        'theta_next_local_deg','theta_next_global_deg', ...
        'EDI_elements','EDI_Gauss_points','newSolve','pass'});

    R=struct();
    R.summary=Summary;
    R.coreMeshControls=coreFP;
    if isfield(candidate,'exteriorMeshControls')
        R.exteriorMeshControls=candidate.exteriorMeshControls;
    end
    R.frozenPhysics=crack_physics_signature(C);
    R.fitTable=fitTable;
    R.EDI=EDI;
    R.prediction=Prediction;
    R.solverInfo=s.solverInfo;
    R.gates=gates;
    R.pass=pass;
    R.newSolve=newSolve;
    R.fastEDI=opt.FastEDI;
    R.checkpointPath=cp;
    R.candidateSource=candidateSource;
    R.frozenSource=frozenSource;
    R.nSegments=nSegments;
    R.thetaSegments=thetaSegments;
    R.thetaSegmentsDeg=thetaSegmentsDeg;
    R.thetaCurrentDeg=thetaCurrentDeg;
    R.thetaCurrent=thetaCurrent;
    R.deltaThetaNextDeg=deltaThetaNextDeg;
    R.deltaThetaNext=deltaThetaNext;
    R.thetaNextDeg=thetaNextDeg;
    R.thetaNext=thetaNext;
    R.thetaNextGlobalDeg=globalThetaNextDeg;
    R.pathFixed=path;
    R.currentIncrementLength=currentIncrement;
    R.nextLegGenerated=false;
    R.interpretation=[ ...
        'One qualified physical LEFM solve of a fixed incremental polyline state. ', ...
        'The measured current-tip KI,KII predict the next MTS turn. ', ...
        'No next segment is generated by this function.'];

    incremental_profile_clock('phase','physical','compact_save');
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
        error('pathsolve:MissingFrozenState', ...
            ['Frozen Stage-I MAT not found. Pass ''FrozenState'',R0 ', ...
             'if the accepted struct is still in memory.']);
    end
    d=load(f);
    assert(isfield(d,'R0')&&isstruct(d.R0), ...
        'pathsolve:BadFrozenState','Frozen MAT must contain R0.');
    R0=d.R0;
    label=f;
end

function fp=local_core_fingerprint(c,increment)
    if isfield(c,'coreMeshControls') && isstruct(c.coreMeshControls)
        fp=c.coreMeshControls;

        % Transitional compatibility with the first diagnostic draft.
        if ~isfield(fp,'expectedEDIElements') && isfield(fp,'expectedEDISupport')
            if abs(fp.scale-1)<=1e-14
                fp.expectedEDIElements=11316;
                fp.expectedOptimizedSupport=10278;
            elseif abs(fp.scale-.5)<=1e-14
                fp.expectedEDIElements=44130;
                fp.expectedOptimizedSupport=40146;
            elseif abs(fp.scale-2)<=1e-14
                fp.expectedEDIElements=2976;
                fp.expectedOptimizedSupport=2700;
            end
        end

        req={'scale','hTip_m','rInner_m','rOuter_m','rCore_m', ...
            'expectedCoreT3','expectedEDIElements', ...
            'expectedOptimizedSupport','expectedNativeSamples'};
        for j=1:numel(req)
            if ~isfield(fp,req{j})
                error('pathsolve:CoreFingerprintField', ...
                    'Candidate coreMeshControls missing %s.',req{j});
            end
        end
    else
        % Backward compatibility for already-saved reference candidates
        % created before coreMeshControls metadata was introduced.
        fp=struct('family','L0_legacy','scale',1, ...
            'hTip_m',0.00675308135*increment, ...
            'rInner_m',.10*increment,'rOuter_m',.65*increment, ...
            'rCore_m',.75*increment,'expectedCoreT3',12678, ...
            'expectedEDIElements',11316, ...
            'expectedOptimizedSupport',10278, ...
            'expectedNativeSamples',[38;55;44;34]);
    end
    fp.expectedNativeSamples=fp.expectedNativeSamples(:);
    if abs(fp.scale-1)<=1e-14
        expected=[12678,11316,10278];
        native=[38;55;44;34];
    elseif abs(fp.scale-.5)<=1e-14
        expected=[49518,44130,40146];
        native=[74;108;86;67];
    elseif abs(fp.scale-2)<=1e-14
        expected=[3318,2976,2700];
        native=[19;28;23;18];
    else
        error('pathsolve:CoreScale','Unsupported structured-core scale %.16g.',fp.scale);
    end
    if fp.expectedCoreT3~=expected(1) || ...
            fp.expectedEDIElements~=expected(2) || ...
            fp.expectedOptimizedSupport~=expected(3) || ...
            ~isequal(fp.expectedNativeSamples,native) || ...
            abs(fp.hTip_m-fp.scale*0.00675308135*increment)>1e-14 || ...
            abs(fp.rInner_m-.10*increment)>1e-14 || ...
            abs(fp.rOuter_m-.65*increment)>1e-14 || ...
            abs(fp.rCore_m-.75*increment)>1e-14
        error('pathsolve:CoreFingerprintMetadata', ...
            'Candidate core-family metadata does not match the audited fingerprint.');
    end
end

function local_require_candidate(c)
    req={'p','t','crack','mat','nSegments','segmentLengths', ...
        'currentIncrementLength','totalCrackLength', ...
        'thetaSegments','thetaSegmentsDeg','thetaCurrent','thetaCurrentDeg', ...
        'path','priorTip','pairedElementIDs','primarySupportElementIDs', ...
        'gates','syntheticGates','scientificallyReadyForIncrementalPhysicalSolve'};
    for k=1:numel(req)
        if ~isfield(c,req{k}) || isempty(c.(req{k}))
            error('pathsolve:CandidateField', ...
                'Candidate missing required field %s.',req{k});
        end
    end
end

function local_require_table_variables(T,names)
    miss=names(~ismember(names,T.Properties.VariableNames));
    if ~isempty(miss)
        error('pathsolve:FrozenFields', ...
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

function local_validate_checkpoint(s,mesh,crack,currentIncrement,ndof,path,thetaSegmentsDeg,mat,C)
    req={'meta','mesh','U','mat','crack','currentIncrement','solverInfo'};
    for k=1:numel(req)
        if ~isfield(s,req{k})
            error('pathsolve:CheckpointField', ...
                'Physical checkpoint missing %s.',req{k});
        end
    end

    assert_crack_checkpoint_physics(s,mat,C,'pathsolve:CheckpointMismatch');
    if ~strcmp(s.meta.stage,'incremental_crack_tip_physical_solved') || ...
            s.meta.nT3Nodes~=size(mesh.coord3,1) || ...
            s.meta.nT3~=size(mesh.connect3,1) || ...
            s.meta.nT6Nodes~=size(mesh.coord,1) || ...
            s.meta.ndof~=ndof || ...
            abs(s.currentIncrement-currentIncrement)>1e-14 || ...
            abs(s.meta.currentIncrement-currentIncrement)>1e-14 || ...
            ~isfield(s.meta,'currentPathPrescribed') || ...
            ~logical(s.meta.currentPathPrescribed) || ...
            ~logical(s.meta.singlePhysicalStateOnly) || ...
            ~logical(s.meta.noThirdLegGenerated) || ...
            ~logical(s.meta.noAngleSweep) || ...
            ~strcmp(s.meta.solver,'pcg_free_dof_spd_sgs') || ...
            ~isequal(s.mesh.connect3,mesh.connect3) || ...
            ~isequal(s.mesh.connect,mesh.connect) || ...
            ~isequal(size(s.mesh.coord),size(mesh.coord)) || ...
            max(abs(s.mesh.coord(:)-mesh.coord(:)))>1e-12 || ...
            max(abs(s.mesh.coord3(:)-mesh.coord3(:)))>1e-12 || ...
            numel(s.U)~=ndof || any(~isfinite(s.U)) || ...
            norm(s.crack.Pmid-crack.Pmid,'fro')>2e-12 || ...
            norm(s.meta.path-path,'fro')>2e-12 || ...
            max(abs(s.meta.thetaSegmentsDeg-thetaSegmentsDeg))>1e-10
        error('pathsolve:CheckpointMismatch', ...
            'Existing checkpoint is not this exact qualified incremental state.');
    end
end

function tf=local_matlab_callable(name)
    tf=exist(name,'file')~=0 || exist(name,'builtin')~=0;
end
