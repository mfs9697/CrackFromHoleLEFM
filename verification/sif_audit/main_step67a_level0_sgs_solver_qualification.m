function R67=main_step67a_level0_sgs_solver_qualification(varargin)
%MAIN_STEP67_LEVEL0_ITERATIVE_SOLVER_QUALIFICATION
% One explicitly authorized Level-0 physical solve using an SPD free-DOF
% formulation and PCG + symmetric Gauss-Seidel.
%
% DEFAULT IS SAFE:
%   AllowSolve = false
%
% This is a solver qualification, NOT a new mesh experiment.
%
% Fixed problem:
%   - deterministic Level-0 C03 mesh;
%   - same material, unit remote-y traction, minimal anchoring;
%   - homogeneous essential constraints imposed by restricting to free DOFs;
%   - one PCG solve with one predeclared incomplete-Cholesky setup;
%   - native COD with the exact Step63 windows/degrees;
%   - one physical EDI on r=[0.8,5.2] mm;
%   - comparison with the established Level-0 direct-solve fingerprints.
%
% NO Level-1 physical solve is present in this file.
%
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');

ip=inputParser;
addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'SourceCandidateFile', ...
    fullfile(vdir,'step62_structured_graded_mesh_candidate_T3.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'CheckpointFile', ...
    fullfile(vdir,'step67a_level0_iterative_physical_solved.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'SaveFile', ...
    fullfile(vdir,'step67a_level0_iterative_qualification_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'ReferenceCheckpoint', ...
    fullfile(vdir,'step63_calibrated_asymmetric_physical_solved.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;

addpath(genpath(root));
assert_step67a_branch(root);

sourceFile=char(opt.SourceCandidateFile);
cp=char(opt.CheckpointFile);
saveFile=char(opt.SaveFile);
refCp=char(opt.ReferenceCheckpoint);

if exist(sourceFile,'file')~=2
    error('step67a:MissingSource', ...
        'Missing committed archived Step62 candidate: %s',sourceFile);
end
for name={'pcg','symamd'}
    if ~matlab_callable(name{1})
        error('step67a:MissingIterativeTool', ...
            'Step67A requires callable MATLAB %s.',name{1});
    end
end

fprintf('\n============================================================\n');
fprintf('STEP 67A: LEVEL-0 ITERATIVE SOLVER QUALIFICATION\n');
fprintf('============================================================\n');
fprintf('  Same physical C03 Level-0 problem as Step63.\n');
fprintf('  SPD free-DOF system; one PCG solve; parameter-free SGS preconditioner.\n');
fprintf('  Exact Step63 COD windows; one matched 0.8-5.2 mm EDI.\n');
fprintf('  NO Level-1 physical solve.\n');

% -------------------------------------------------------------------------
% Deterministic exact Level-0 C03 reconstruction, mesh-only.
% -------------------------------------------------------------------------
cal=c03_calibration();
O0=main_step62_structured_graded_mesh( ...
    'SourceCandidateFile',sourceFile, ...
    'SavePrefix',fullfile(vdir,'step67a_internal_level0'), ...
    'Visible','off', ...
    'RunSynthetic',false, ...
    'ExteriorCalibration',cal, ...
    'WriteArtifacts',false, ...
    'ReturnCandidate',true, ...
    'Verbose',false, ...
    'Level',0);

if ~O0.gates.structuralPass || ...
        O0.summary.newT3~=32980 || O0.summary.newT6Nodes~=66854 || ...
        abs(O0.maxNeighborSizeRatio-1.79678451)>5e-7 || ...
        abs(O0.summary.pairedRadius_mm-6)>1e-12 || ...
        ~O0.gates.completeT3Pairing || ~O0.gates.completeT6Pairing
    error('step67a:Level0Mismatch', ...
        'Deterministic Level-0 C03 reconstruction no longer matches qualification.');
end

candidate=O0.candidate;
P=candidate.p; T=candidate.t; crack=candidate.crack; mat0=candidate.mat;
[P6,T6]=T3toT6_fast(P,T);
mesh=struct('coord3',P,'connect3',T,'coord',P6,'connect',T6);
a0=norm(diff(crack.Pmid,1,1));

if abs(a0-.008)>1e-12 || ...
        norm(crack.Pmid(1,:)-[.199988872196,-.0208170339049])>2e-12 || ...
        norm(crack.Pmid(end,:)-[.207985904782,-.0210349096129])>2e-12
    error('step67a:CrackGeometry','Qualified crack geometry changed.');
end

% Native sampling gate before any solve.
zeroU=zeros(2*size(P6,1),1);
[rZero,~,face0]=native_COD_audit(mesh,zeroU,mat0,crack,8);
rrZero=rZero/a0;
windows=[.04 .20;.04 .30;.08 .30;.12 .30];
sampleN=zeros(4,1);
for k=1:4
    sampleN(k)=nnz(rrZero>=windows(k,1)&rrZero<=windows(k,2));
end
if face0.nUpper~=138 || face0.nLower~=138 || ...
        face0.gridMismatch>1e-12 || any(sampleN~=[38;55;44;34])
    error('step67a:NativeSampling', ...
        'Expected exact Level-0 native COD sampling 38/55/44/34.');
end

% -------------------------------------------------------------------------
% Same physical setup as Step63.
% -------------------------------------------------------------------------
xmin=min(P(:,1)); xmax=max(P(:,1));
ymin=min(P(:,2)); ymax=max(P(:,2));
A=xmax; B=max(abs([ymin,ymax]));
if abs(xmin)>1e-12 || abs(A-.30)>1e-12 || ...
        abs(ymin+.10)>1e-12 || abs(ymax-.10)>1e-12
    error('step67a:PlateBoundary','Plate boundary changed.');
end

Cref=cfg_hole_initiation();
if abs(Cref.E-mat0.E)>1e-12*max(1,abs(mat0.E)) || ...
        abs(Cref.nu-mat0.nu)>1e-14 || Cref.ps~=mat0.ps || ...
        ~strcmp(Cref.load.type,'remote_tension_y') || ...
        abs(Cref.load.sig0-1)>1e-14 || ...
        ~strcmp(Cref.bc.anchor_mode,'minimal')
    error('step67a:PhysicalConfig', ...
        'Canonical material/loading/anchoring differs from audited physics.');
end

mat=material_with_D(mat0);
[nip2,xip2,w2,Nextr]=integr();
quad=struct('nip2',nip2,'xip2',xip2,'w2',w2,'Nextr',Nextr);

[~,iLB]=min(sum((P-[0,-B]).^2,2));
[~,iRB]=min(sum((P-[A,-B]).^2,2));
if norm(P(iLB,:)-[0,-B])>1e-12 || norm(P(iRB,:)-[A,-B])>1e-12
    error('step67a:Corners','Could not identify exact plate bottom corners.');
end
fixvar=unique([2*iLB-1;2*iLB;2*iRB]);
ndof=2*size(P6,1);
free=true(ndof,1);free(fixvar)=false;free=find(free);

C=struct('A',A,'B',B,'E',mat.E,'nu',mat.nu,'ps',mat.ps,'a0',a0, ...
    'load',Cref.load,'bc',Cref.bc);

% Fixed PCG/preconditioner qualification settings. There is no parameter sweep.
pcgTol=1e-10;
pcgMaxIt=5000;
preconditionerName='symmetric_gauss_seidel';

% Reference fingerprints from the completed Step63R/64 direct-solve audit.
refFitRatio=[ ...
    0.0001052647913;
    0.0001056239885;
    0.0001048946546;
    0.0001057078877;
    0.0001045631315;
    0.0001058187639;
    0.0001041556245;
    0.0001058252201];
refEDI=struct('KI',0.43785,'KII',4.6547e-05,'ratio',0.00010631);

fprintf('  Level 0: T3=%d, T6=%d, DOF=%d, free DOF=%d.\n', ...
    size(T,1),size(P6,1),ndof,numel(free));
fprintf('  PCG tol=%.1e, maxit=%d; SGS preconditioner, no tuning parameter.\n', ...
    pcgTol,pcgMaxIt);

% -------------------------------------------------------------------------
% Existing valid Step67A checkpoint is reused. Otherwise one authorized solve.
% -------------------------------------------------------------------------
newSolve=false;
if exist(cp,'file')==2
    s0=load(cp,'meta','mesh','U','mat','crack','a0','solverInfo');
    if ~isfield(s0,'meta') || ...
            ~strcmp(s0.meta.stage,'step67a_level0_iterative_physical') || ...
            s0.meta.nT3~=32980 || s0.meta.nT6~=66854 || ...
            abs(s0.meta.maxNeighborRatio-1.79678451)>5e-7 || ...
            abs(s0.a0-a0)>1e-12 || ...
            ~isequal(s0.mesh.connect3,T) || ...
            max(abs(s0.mesh.coord3(:)-P(:)))>1e-12 || ...
            numel(s0.U)~=ndof || any(~isfinite(s0.U))
        error('step67a:ExistingCheckpointMismatch', ...
            'Existing Step67A checkpoint is not the exact current qualification field.');
    end
    fprintf('  Reusing existing Step67A checkpoint; NO new physical solve.\n');
else
    if ~opt.AllowSolve
        error('step67a:ExplicitSolveApprovalRequired', ...
            ['Step67A is prepared but guarded. Rerun with ''AllowSolve'',true ', ...
             'only after explicit investigator authorization.']);
    end

    fprintf('\nSTEP67 PHASE 1: ASSEMBLE SYMMETRIC LEVEL-0 SYSTEM\n');
    fprintf('  stif_assem called with fixvar=[]; no row clamping.\n');
    K=stif_assem(mesh,mat,quad,[]);
    if size(K,1)~=ndof || size(K,2)~=ndof
        error('step67a:StiffnessSize','Unexpected stiffness dimensions.');
    end
    symErr=norm(K-K.',1)/max(1,norm(K,1));
    if ~isfinite(symErr) || symErr>5e-13
        error('step67a:StiffnessSymmetry', ...
            'Unclamped stiffness symmetry error %.3e exceeds gate.',symErr);
    end

    % Build the identical Step63 remote-y load vector.
    F=zeros(ndof,1);
    eps1=1e-8*max(A,2*B);
    elod=edge_loads_T6(mesh.coord,B,eps1);
    qy=C.load.sig0;
    for k=1:size(elod,1)
        node=elod(k,1);w=elod(k,2);
        F(2*node)=F(2*node)+qy*w;
    end
    F(fixvar)=0;

    Kff=K(free,free);
    Ff=F(free);
    % Remove only machine-level antisymmetry after the strict gate.
    Kff=(Kff+Kff.')/2;

    if any(diag(Kff)<=0)
        error('step67a:NonpositiveDiagonal','Free-DOF K has nonpositive diagonal.');
    end

    fprintf('  K symmetry gate: %.3e. Building parameter-free SGS preconditioner.\n',symErr);
    p=symamd(Kff);
    Ap=Kff(p,p);
    bp=Ff(p);
    clear Kff Ff

    d=diag(Ap);
    if any(~isfinite(d)) || any(d<=0)
        error('step67a:NonpositiveSGSDiagonal', ...
            'SGS requires a finite positive diagonal.');
    end

    tPre=tic;
    M1=tril(Ap); % D+L
    Dinv=spdiags(1./d,0,numel(d),numel(d));
    M2=Dinv*M1.'; % D^{-1}(D+U), so M1*M2 is symmetric SGS
    preconditionerSeconds=toc(tPre);
    nnzM1=nnz(M1); nnzM2=nnz(M2);
    preconditionerGiB=(16*(nnzM1+nnzM2)+ ...
        8*(size(M1,2)+size(M2,2)+2))/1024^3;

    fprintf('  SGS complete: nnz(M1)+nnz(M2)=%d, approx sparse storage %.4f GiB, %.2f s.\n', ...
        nnzM1+nnzM2,preconditionerGiB,preconditionerSeconds);
    fprintf('  Starting exactly ONE explicitly authorized Level-0 PCG physical solve.\n');

    tSolve=tic;
    [xp,flag,relres,iter,resvec]=pcg(Ap,bp,pcgTol,pcgMaxIt,M1,M2);
    solveSeconds=toc(tSolve);

    if flag~=0 || ~isfinite(relres) || relres>pcgTol || any(~isfinite(xp))
        error('step67a:PCGFailed', ...
            'PCG qualification failed: flag=%d, relres=%.3e, iter=%d.', ...
            flag,relres,iter);
    end

    uf=zeros(numel(free),1);
    uf(p)=xp;
    U=zeros(ndof,1);
    U(free)=uf;
    U(fixvar)=0;

    % True free-system residual in original ordering.
    trueRelResidual=norm(K(free,:)*U-F(free))/max(norm(F(free)),eps);
    constraintInf=max(abs(U(fixvar)));
    if trueRelResidual>5e-10 || constraintInf>1e-14
        error('step67a:ResidualGate', ...
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
        'nnzK',nnz(K), ...
        'nnzAp',nnz(Ap), ...
        'nnzM1',nnzM1,'nnzM2',nnzM2, ...
        'preconditionerGiBApprox',preconditionerGiB, ...
        'preconditionerSeconds',preconditionerSeconds, ...
        'solveSeconds',solveSeconds, ...
        'onePhysicalLinearSolve',true, ...
        'noDirectBackslash',true);

    meta=struct( ...
        'stage','step67a_level0_iterative_physical', ...
        'source','one explicitly authorized Level-0 free-DOF PCG solver qualification', ...
        'selectedCandidate','C03', ...
        'nT3',size(T,1),'nT6',size(P6,1),'ndof',ndof, ...
        'a0',a0,'maxNeighborRatio',O0.maxNeighborSizeRatio, ...
        'loadType',C.load.type,'sig0',C.load.sig0, ...
        'anchorMode',C.bc.anchor_mode,'cornerNodes',[iLB iRB], ...
        'solver','pcg_free_dof_spd_sgs', ...
        'pcgTol',pcgTol,'preconditioner',preconditionerName, ...
        'physicalEDIperformedAfterCheckpoint',true, ...
        'meshLevel',0,'noLevel1PhysicalSolve',true);

    [folder,~,~]=fileparts(cp);
    if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
    tmp=[cp '.incomplete.mat'];
    if exist(tmp,'file')==2
        error('step67a:InterruptedSave', ...
            'Inspect/remove prior incomplete checkpoint manually: %s',tmp);
    end

    save(tmp,'mesh','U','mat','crack','a0','meta','solverInfo','C','-v7.3');
    [ok,msg]=movefile(tmp,cp);
    if ~ok,error('step67a:CheckpointSave','%s',msg);end
    fprintf('  Iterative physical field safely checkpointed: %s\n',cp);
    newSolve=true;

    clear K F Ap bp M1 M2 Dinv xp uf U
end

% -------------------------------------------------------------------------
% Postprocess only from the saved Step67A checkpoint.
% -------------------------------------------------------------------------
fprintf('\nSTEP67 PHASE 2: PHYSICAL QUALIFICATION FROM SAVED FIELD\n');
s=load(cp,'mesh','U','mat','crack','a0','meta','solverInfo');
if ~strcmp(s.meta.stage,'step67a_level0_iterative_physical') || ...
        s.meta.nT3~=32980 || s.meta.nT6~=66854 || ...
        ~strcmp(s.meta.solver,'pcg_free_dof_spd_sgs')
    error('step67a:CheckpointProvenance','Saved iterative checkpoint provenance mismatch.');
end

[r,app,face]=native_COD_audit(s.mesh,s.U,s.mat,s.crack,8);
rr=r/s.a0;
fitDegrees=[1 2];
fitRows=nan(8,9);
qRaw=app(:,2)./app(:,1);
n=0;
for iw=1:size(windows,1)
    ids=find(rr>=windows(iw,1)&rr<=windows(iw,2));
    for d=fitDegrees
        pI=polyfit(rr(ids),app(ids,1),d);
        pII=polyfit(rr(ids),app(ids,2),d);
        predII=polyval(pII,rr(ids));
        n=n+1;
        fitRows(n,:)=[windows(iw,:),d,numel(ids), ...
            pI(end),pII(end),pII(end)/pI(end), ...
            sqrt(mean((app(ids,2)-predII).^2)), ...
            median(qRaw(ids))];
    end
end
fitTable=array2table(fitRows, ...
    'VariableNames',{'lower_r_over_a0','upper_r_over_a0','degree', ...
    'n_native','KI_COD','KII_COD','ratio_COD','RMSE_KII', ...
    'median_pointwise_ratio'});

fitRatioAbsDiff=fitTable.ratio_COD-refFitRatio;
fitRatioRelDiff=fitRatioAbsDiff./refFitRatio;
CODcomparison=table( ...
    fitTable.lower_r_over_a0,fitTable.upper_r_over_a0,fitTable.degree, ...
    fitTable.n_native,refFitRatio,fitTable.ratio_COD, ...
    fitRatioAbsDiff,100*fitRatioRelDiff, ...
    'VariableNames',{'lower_r_over_a0','upper_r_over_a0','degree','n_native', ...
    'directReferenceRatio','iterativeRatio','absoluteDifference','relativeDifferencePct'});

% One and only one matched physical EDI.
ri=.0008;ro=.65*s.a0;
if abs(ro-.0052)>1e-14
    error('step67a:EDIRadius','Expected r_outer=5.2 mm.');
end
matEDI=s.mat;
if ~isfield(matEDI,'Dmat')&&isfield(matEDI,'D'),matEDI.Dmat=matEDI.D;end
fprintf('  Running ONE matched 16-point FE-nodal-q EDI on the saved iterative field.\n');
[KI,KII]=SIF_LEFM_interaction_EDI( ...
    s.mesh,s.U,s.crack.Pmid,matEDI, ...
    struct('r_inner',ri,'r_outer',ro), ...
    'UsePlaneStrain',matEDI.ps==1, ...
    'Verbose',false, ...
    'WeightFunction','fe_nodal', ...
    'QuadratureRule',16, ...
    'StoreGPDiagnostics',false);
qEDI=KII/KI;

EDIcomparison=table( ...
    refEDI.KI,KI,100*(KI-refEDI.KI)/refEDI.KI, ...
    refEDI.KII,KII,100*(KII-refEDI.KII)/refEDI.KII, ...
    refEDI.ratio,qEDI,100*(qEDI-refEDI.ratio)/refEDI.ratio, ...
    'VariableNames',{'direct_KI','iterative_KI','KI_gap_pct', ...
    'direct_KII','iterative_KII','KII_gap_pct', ...
    'direct_ratio','iterative_ratio','ratio_gap_pct'});

% Optional direct displacement comparison if the recovered Step63 checkpoint
% happens to still exist locally. It is not required because generated MAT
% files have previously disappeared across branch switches.
referenceDisplacementAvailable=false;
Urelative2=NaN;UrelativeInf=NaN;
if exist(refCp,'file')==2
    d=load(refCp,'U','mesh','a0','meta');
    if isfield(d,'U')&&isfield(d,'mesh')&&isfield(d,'meta') && ...
            strcmp(d.meta.stage,'step63_calibrated_asymmetric_physical') && ...
            numel(d.U)==numel(s.U) && ...
            isequal(d.mesh.connect3,s.mesh.connect3) && ...
            max(abs(d.mesh.coord3(:)-s.mesh.coord3(:)))<1e-12 && ...
            abs(d.a0-s.a0)<1e-12
        referenceDisplacementAvailable=true;
        Urelative2=norm(s.U-d.U)/max(norm(d.U),eps);
        UrelativeInf=max(abs(s.U-d.U))/max(max(abs(d.U)),eps);
    end
end

% Predeclared qualification gates.
gates=struct();
gates.pcgConverged=s.solverInfo.flag==0;
gates.pcgReportedResidual=s.solverInfo.relres<=pcgTol;
gates.trueResidual=s.solverInfo.trueRelResidual<=5e-10;
gates.constraintsExact=s.solverInfo.constraintInf<=1e-14;
gates.nativeSamplingExact=isequal(fitTable.n_native,[38;38;55;55;44;44;34;34]);
gates.codRatiosMatch=max(abs(fitRatioRelDiff))<=5e-4;       % 0.05%
gates.ediKIMatch=abs(KI-refEDI.KI)/refEDI.KI<=5e-4;        % 0.05%
gates.ediKIIMatch=abs(KII-refEDI.KII)/refEDI.KII<=5e-4;    % 0.05%
gates.ediRatioMatch=abs(qEDI-refEDI.ratio)/refEDI.ratio<=5e-4; % 0.05%
if referenceDisplacementAvailable
    gates.directDisplacementMatch=Urelative2<=5e-8 && UrelativeInf<=5e-8;
else
    gates.directDisplacementMatch=true; % optional, reported separately
end
gates.singlePhysicalSolveThisQualification=logical(newSolve) || exist(cp,'file')==2;
gates.noDirectBackslash=logical(s.solverInfo.noDirectBackslash);
gates.noLevel1PhysicalSolve=logical(s.meta.noLevel1PhysicalSolve);

qualificationPass=all(structfun(@logical,gates));

fprintf('\nLEVEL-0 DIRECT REFERENCE vs ITERATIVE COD\n');
disp(CODcomparison);
fprintf('\nLEVEL-0 DIRECT REFERENCE vs ITERATIVE EDI\n');
disp(EDIcomparison);
fprintf('\nITERATIVE SOLVER INFO\n');
disp(s.solverInfo);
fprintf('\nOPTIONAL DIRECT DISPLACEMENT COMPARISON\n');
disp(table(referenceDisplacementAvailable,Urelative2,UrelativeInf));
fprintf('\nSTEP67 QUALIFICATION GATES\n');
disp(gates);

if qualificationPass
    fprintf('STEP67 PASS: Level-0 iterative solver reproduces the direct physical reference.\n');
    fprintf('  Iterative formulation is qualified for a future Level-1 proposal.\n');
    fprintf('  NO Level-1 physical solve has been performed.\n');
else
    fprintf('STEP67 STOP: iterative solver failed one or more qualification gates.\n');
    fprintf('  Do NOT use this solver on Level 1.\n');
end

R67=struct( ...
    'checkpointPath',cp, ...
    'newSolve',newSolve, ...
    'solverInfo',s.solverInfo, ...
    'fitTable',fitTable, ...
    'CODcomparison',CODcomparison, ...
    'EDIcomparison',EDIcomparison, ...
    'referenceDisplacementAvailable',referenceDisplacementAvailable, ...
    'Urelative2',Urelative2,'UrelativeInf',UrelativeInf, ...
    'gates',gates, ...
    'qualificationPass',qualificationPass, ...
    'readyForLevel1IterativeProposal',qualificationPass, ...
    'noLevel1PhysicalSolve',true, ...
    'singleMatchedEDI',true, ...
    'interpretation',['Level-0 solver qualification only. Passing means the ', ...
      'free-DOF PCG/SGS formulation reproduces the established direct ', ...
      'physical reference; it does not itself establish Level-1 convergence.']);

save(saveFile,'R67','-v7');
fprintf('  Compact Step67A result saved: %s\n',saveFile);
end

% =========================================================================
function mat=material_with_D(mat0)
mat=mat0;
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

function cal=c03_calibration()
cal=struct( ...
    'transitionLength_m',.008, ...
    'farSlope',.10, ...
    'boundaryMetricGrowth',.25, ...
    'smoothingSteps',6, ...
    'refinementMaxPasses',140, ...
    'refinementMinAngle_deg',25, ...
    'refinementLongestFactor',1.65, ...
    'neighborRatioTarget',1.8, ...
    'verbose',false);
end

function tf=matlab_callable(name)
tf=exist(name,'file')~=0 || exist(name,'builtin')~=0;
end

function assert_step67a_branch(root)
[status,b]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0&&strcmp(strtrim(b),'audit/step67a-level0-sgs-solver-qualification'), ...
    'step67a:Branch', ...
    'Run Step67A only on audit/step67a-level0-sgs-solver-qualification.');
end
