function R68=main_step68_level1_sgs_physical_convergence(varargin)
%MAIN_STEP68_LEVEL1_SGS_PHYSICAL_CONVERGENCE
% One explicitly authorized Level-1 physical solve using the SGS-PCG
% formulation qualified on Level 0 in Step67A.
%
% DEFAULT IS SAFE:
%   AllowSolve = false
%
% This is the first direct mesh-family convergence experiment:
%   Level 0: accepted C03 physical reference (Step63R/64/67A)
%   Level 1: same C03 family, structured scale 0.5 (Step65)
%
% Fixed numerical protocol:
%   - deterministic Level-1 C03 reconstruction from committed Step62 source;
%   - same material, unit remote-y traction, minimal anchoring;
%   - homogeneous essential constraints imposed on free DOFs;
%   - symamd ordering;
%   - parameter-free symmetric Gauss-Seidel preconditioner;
%   - PCG tol=1e-10, maxit=5000;
%   - immediate solved-field checkpoint;
%   - same four COD windows, polynomial degrees 1 and 2;
%   - exactly one 16-point FE-nodal-q EDI on r=[0.8,5.2] mm.
%
% There is NO radius sweep, solver tuning, or Level-0 physical solve here.
%
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');

ip=inputParser;
addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'SourceCandidateFile', ...
    fullfile(vdir,'step62_structured_graded_mesh_candidate_T3.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'CheckpointFile', ...
    fullfile(vdir,'step68_level1_sgs_physical_solved.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'SaveFile', ...
    fullfile(vdir,'step68_level1_sgs_convergence_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'Level0ReferenceFile', ...
    fullfile(vdir,'step67a_level0_sgs_qualification_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;

addpath(genpath(root));
assert_step68_branch(root);

sourceFile=char(opt.SourceCandidateFile);
cp=char(opt.CheckpointFile);
saveFile=char(opt.SaveFile);
level0File=char(opt.Level0ReferenceFile);

if exist(sourceFile,'file')~=2
    error('step68:MissingSource', ...
        'Missing committed archived Step62 candidate: %s',sourceFile);
end
for name={'pcg','symamd'}
    if ~matlab_callable(name{1})
        error('step68:MissingIterativeTool', ...
            'Step68 requires callable MATLAB %s.',name{1});
    end
end

fprintf('\n============================================================\n');
fprintf('STEP 68: LEVEL-1 SGS-PCG PHYSICAL CONVERGENCE SOLVE\n');
fprintf('============================================================\n');
fprintf('  ONE Level-1 physical solve only; no Level-0 solve.\n');
fprintf('  Same C03 family, scale 0.5, same physical geometry/support.\n');
fprintf('  SGS-PCG exactly as qualified in Step67A.\n');
fprintf('  Same COD windows; one matched 0.8-5.2 mm EDI.\n');
fprintf('  NO radius sweep. NO solver tuning.\n');

% -------------------------------------------------------------------------
% Deterministic exact Level-1 C03 reconstruction, mesh-only.
% -------------------------------------------------------------------------
cal=c03_calibration();
O1=main_step62_structured_graded_mesh( ...
    'SourceCandidateFile',sourceFile, ...
    'SavePrefix',fullfile(vdir,'step68_internal_level1'), ...
    'Visible','off', ...
    'RunSynthetic',false, ...
    'ExteriorCalibration',cal, ...
    'WriteArtifacts',false, ...
    'ReturnCandidate',true, ...
    'Verbose',false, ...
    'Level',1);

expectedMaxRatio=1.79938154;
if ~O1.gates.structuralPass || ...
        O1.summary.newT3~=122691 || O1.summary.newT6Nodes~=246701 || ...
        abs(O1.design.scale-.5)>1e-14 || ...
        abs(O1.maxNeighborSizeRatio-expectedMaxRatio)>5e-7 || ...
        abs(O1.summary.pairedRadius_mm-6)>1e-12 || ...
        ~O1.gates.completeT3Pairing || ~O1.gates.completeT6Pairing || ...
        ~O1.gates.candidateSupportInsidePairedRegion || ...
        ~O1.gates.exteriorOutsideSupport || ...
        ~O1.gates.tipFanSixTriangles
    error('step68:Level1Mismatch', ...
        'Deterministic Level-1 C03 reconstruction no longer matches Step65.');
end

candidate=O1.candidate;
level1MeshSummary=struct( ...
    'nT3',O1.summary.newT3, ...
    'nT6',O1.summary.newT6Nodes, ...
    'maxNeighborRatio',O1.maxNeighborSizeRatio, ...
    'tipMedian_mm',O1.summary.newTipMedian_mm, ...
    'pairedRadius_mm',O1.summary.pairedRadius_mm, ...
    'patchMinAngle_deg',O1.summary.patchMinAngle_deg, ...
    'exteriorMinAngle_deg',O1.summary.exteriorMinAngle_deg);
clear O1

P=candidate.p; T=candidate.t; crack=candidate.crack; mat0=candidate.mat;
clear candidate
[P6,T6]=T3toT6_fast(P,T);
mesh=struct('coord3',P,'connect3',T,'coord',P6,'connect',T6);
a0=norm(diff(crack.Pmid,1,1));

if size(T,1)~=122691 || size(P6,1)~=246701
    error('step68:MeshCount','Unexpected exact Level-1 mesh counts.');
end
if abs(a0-.008)>1e-12 || ...
        norm(crack.Pmid(1,:)-[.199988872196,-.0208170339049])>2e-12 || ...
        norm(crack.Pmid(end,:)-[.207985904782,-.0210349096129])>2e-12
    error('step68:CrackGeometry','Qualified crack geometry changed.');
end

% Native sampling gate before solve.
zeroU=zeros(2*size(P6,1),1);
[rZero,~,face0]=native_COD_audit(mesh,zeroU,mat0,crack,8);
rrZero=rZero/a0;
windows=[.04 .20;.04 .30;.08 .30;.12 .30];
sampleN=zeros(4,1);
for k=1:4
    sampleN(k)=nnz(rrZero>=windows(k,1)&rrZero<=windows(k,2));
end
clear zeroU rZero rrZero
if face0.nUpper~=272 || face0.nLower~=272 || ...
        face0.gridMismatch>1e-12 || any(sampleN~=[74;108;86;67])
    error('step68:NativeSampling', ...
        'Expected exact Level-1 native sampling 74/108/86/67.');
end

% -------------------------------------------------------------------------
% Same physical setup as Level 0.
% -------------------------------------------------------------------------
xmin=min(P(:,1)); xmax=max(P(:,1));
ymin=min(P(:,2)); ymax=max(P(:,2));
A=xmax; B=max(abs([ymin,ymax]));
if abs(xmin)>1e-12 || abs(A-.30)>1e-12 || ...
        abs(ymin+.10)>1e-12 || abs(ymax-.10)>1e-12
    error('step68:PlateBoundary','Plate boundary changed.');
end

Cref=cfg_hole_initiation();
if abs(Cref.E-mat0.E)>1e-12*max(1,abs(mat0.E)) || ...
        abs(Cref.nu-mat0.nu)>1e-14 || Cref.ps~=mat0.ps || ...
        ~strcmp(Cref.load.type,'remote_tension_y') || ...
        abs(Cref.load.sig0-1)>1e-14 || ...
        ~strcmp(Cref.bc.anchor_mode,'minimal')
    error('step68:PhysicalConfig', ...
        'Canonical material/loading/anchoring differs from audited physics.');
end

mat=material_with_D(mat0);
[nip2,xip2,w2,Nextr]=integr();
quad=struct('nip2',nip2,'xip2',xip2,'w2',w2,'Nextr',Nextr);

[~,iLB]=min(sum((P-[0,-B]).^2,2));
[~,iRB]=min(sum((P-[A,-B]).^2,2));
if norm(P(iLB,:)-[0,-B])>1e-12 || norm(P(iRB,:)-[A,-B])>1e-12
    error('step68:Corners','Could not identify exact plate bottom corners.');
end
fixvar=unique([2*iLB-1;2*iLB;2*iRB]);
ndof=2*size(P6,1);
freeMask=true(ndof,1); freeMask(fixvar)=false; free=find(freeMask);
clear freeMask

C=struct('A',A,'B',B,'E',mat.E,'nu',mat.nu,'ps',mat.ps,'a0',a0, ...
    'load',Cref.load,'bc',Cref.bc);

pcgTol=1e-10;
pcgMaxIt=5000;
preconditionerName='symmetric_gauss_seidel';

fprintf('  Level 1: T3=%d, T6=%d, DOF=%d, free DOF=%d.\n', ...
    size(T,1),size(P6,1),ndof,numel(free));
fprintf('  PCG tol=%.1e, maxit=%d; SGS preconditioner, no tuning parameter.\n', ...
    pcgTol,pcgMaxIt);

% -------------------------------------------------------------------------
% Resolve Level-0 physical references.
% Exact Step67A result is preferred if its local small-data file survived.
% Otherwise use the recorded audit fingerprints, explicitly marked rounded
% for EDI and high-precision for the eight COD ratios.
% -------------------------------------------------------------------------
fallbackFitRatio=[ ...
    0.0001052647913;
    0.0001056239885;
    0.0001048946546;
    0.0001057078877;
    0.0001045631315;
    0.0001058187639;
    0.0001041556245;
    0.0001058252201];
fallbackEDI=struct('KI',0.43785,'KII',4.6547e-05,'ratio',0.00010631);

level0ReferenceSource="embedded_audit_fingerprint";
level0FitRatio=fallbackFitRatio;
level0EDI=fallbackEDI;
level0ReferencePrecision="COD high-precision audit fingerprint; EDI displayed/rounded";

if exist(level0File,'file')==2
    z=load(level0File,'R67a');
    if isfield(z,'R67a') && ...
            isfield(z.R67a,'qualificationPass') && z.R67a.qualificationPass && ...
            isfield(z.R67a,'fitTable') && height(z.R67a.fitTable)==8 && ...
            isfield(z.R67a,'EDIcomparison') && height(z.R67a.EDIcomparison)==1
        level0FitRatio=z.R67a.fitTable.ratio_COD;
        level0EDI=struct( ...
            'KI',z.R67a.EDIcomparison.iterative_KI, ...
            'KII',z.R67a.EDIcomparison.iterative_KII, ...
            'ratio',z.R67a.EDIcomparison.iterative_ratio);
        level0ReferenceSource="exact_saved_step67a";
        level0ReferencePrecision="exact saved Step67A iterative result";
    end
end

% -------------------------------------------------------------------------
% Existing valid Step68 checkpoint is reused. Otherwise exactly one solve.
% -------------------------------------------------------------------------
newSolve=false;
if exist(cp,'file')==2
    s0=load(cp,'meta','mesh','U','mat','crack','a0','solverInfo');
    if ~isfield(s0,'meta') || ...
            ~strcmp(s0.meta.stage,'step68_level1_sgs_physical') || ...
            s0.meta.nT3~=122691 || s0.meta.nT6~=246701 || ...
            abs(s0.meta.maxNeighborRatio-expectedMaxRatio)>5e-7 || ...
            abs(s0.a0-a0)>1e-12 || ...
            ~isequal(s0.mesh.connect3,T) || ...
            max(abs(s0.mesh.coord3(:)-P(:)))>1e-12 || ...
            numel(s0.U)~=ndof || any(~isfinite(s0.U))
        error('step68:ExistingCheckpointMismatch', ...
            'Existing Step68 checkpoint is not the exact current Level-1 field.');
    end
    fprintf('  Reusing existing Step68 checkpoint; NO new physical solve.\n');
else
    if ~opt.AllowSolve
        error('step68:ExplicitSolveApprovalRequired', ...
            ['Step68 is prepared but guarded. Rerun with ''AllowSolve'',true ', ...
             'only after explicit investigator authorization.']);
    end

    fprintf('\nSTEP68 PHASE 1: ASSEMBLE SYMMETRIC LEVEL-1 SYSTEM\n');
    fprintf('  stif_assem called with fixvar=[]; no row clamping.\n');
    [availBefore,totalPhysical]=available_memory_gib();
    if isfinite(availBefore)
        fprintf('  MATLAB-reported available physical memory before assembly: %.3f GiB.\n',availBefore);
    end

    K=stif_assem(mesh,mat,quad,[]);
    if size(K,1)~=ndof || size(K,2)~=ndof
        error('step68:StiffnessSize','Unexpected stiffness dimensions.');
    end
    symErr=norm(K-K.',1)/max(1,norm(K,1));
    if ~isfinite(symErr) || symErr>5e-13
        error('step68:StiffnessSymmetry', ...
            'Unclamped stiffness symmetry error %.3e exceeds gate.',symErr);
    end

    % Identical unit remote-y loading.
    F=zeros(ndof,1);
    eps1=1e-8*max(A,2*B);
    elod=edge_loads_T6(mesh.coord,B,eps1);
    qy=C.load.sig0;
    for k=1:size(elod,1)
        node=elod(k,1); w=elod(k,2);
        F(2*node)=F(2*node)+qy*w;
    end
    F(fixvar)=0;

    Kff=K(free,free);
    Ff=F(free);
    Kff=(Kff+Kff.')/2;
    if any(diag(Kff)<=0)
        error('step68:NonpositiveDiagonal','Free-DOF K has nonpositive diagonal.');
    end

    fprintf('  K symmetry gate: %.3e. Applying symamd permutation.\n',symErr);
    p=symamd(Kff);
    Ap=Kff(p,p);
    bp=Ff(p);

    % Important memory release before SGS construction.
    nnzK=nnz(K);
    nnzAp=nnz(Ap);
    clear K Kff F Ff

    d=diag(Ap);
    if any(~isfinite(d)) || any(d<=0)
        error('step68:NonpositiveSGSDiagonal', ...
            'SGS requires a finite positive diagonal.');
    end

    [availPre,totalPhysical2]=available_memory_gib();
    if isfinite(totalPhysical2),totalPhysical=totalPhysical2;end

    tPre=tic;
    M1=tril(Ap); % D+L
    Dinv=spdiags(1./d,0,numel(d),numel(d));
    M2=Dinv*M1.'; % D^{-1}(D+U)
    preconditionerSeconds=toc(tPre);
    nnzM1=nnz(M1); nnzM2=nnz(M2);
    preconditionerGiB=(16*(nnzM1+nnzM2)+ ...
        8*(size(M1,2)+size(M2,2)+2))/1024^3;

    [availAfterPre,totalPhysical3]=available_memory_gib();
    if isfinite(totalPhysical3),totalPhysical=totalPhysical3;end

    fprintf('  SGS complete: nnz(M1)+nnz(M2)=%d, approx sparse storage %.4f GiB, %.2f s.\n', ...
        nnzM1+nnzM2,preconditionerGiB,preconditionerSeconds);
    if isfinite(availAfterPre)
        fprintf('  Available physical memory after SGS construction: %.3f GiB.\n',availAfterPre);
    end
    fprintf('  Starting exactly ONE explicitly authorized Level-1 PCG physical solve.\n');

    tSolve=tic;
    [xp,flag,relres,iter,resvec]=pcg(Ap,bp,pcgTol,pcgMaxIt,M1,M2);
    solveSeconds=toc(tSolve);

    if flag~=0 || ~isfinite(relres) || relres>pcgTol || any(~isfinite(xp))
        error('step68:PCGFailed', ...
            ['Level-1 PCG failed with fixed qualified settings: ', ...
             'flag=%d, relres=%.3e, iter=%d. No retuning in Step68.'], ...
            flag,relres,iter);
    end

    trueRelResidual=norm(Ap*xp-bp)/max(norm(bp),eps);
    if trueRelResidual>5e-10
        error('step68:ResidualGate', ...
            'True permuted free-system residual %.3e exceeds gate.',trueRelResidual);
    end

    uf=zeros(numel(free),1);
    uf(p)=xp;
    U=zeros(ndof,1);
    U(free)=uf;
    U(fixvar)=0;
    constraintInf=max(abs(U(fixvar)));
    if constraintInf>1e-14
        error('step68:ConstraintGate','Constraint infinity norm %.3e.',constraintInf);
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
        'nnzK',nnzK,'nnzAp',nnzAp, ...
        'nnzM1',nnzM1,'nnzM2',nnzM2, ...
        'preconditionerGiBApprox',preconditionerGiB, ...
        'preconditionerSeconds',preconditionerSeconds, ...
        'solveSeconds',solveSeconds, ...
        'availableGiBBeforeAssembly',availBefore, ...
        'availableGiBBeforePreconditioner',availPre, ...
        'availableGiBAfterPreconditioner',availAfterPre, ...
        'totalPhysicalGiB',totalPhysical, ...
        'onePhysicalLinearSolve',true, ...
        'noDirectBackslash',true, ...
        'noSolverTuning',true);

    meta=struct( ...
        'stage','step68_level1_sgs_physical', ...
        'source','one explicitly authorized Level-1 SGS-PCG physical convergence solve', ...
        'selectedCandidate','C03', ...
        'meshLevel',1,'scale',.5, ...
        'nT3',size(T,1),'nT6',size(P6,1),'ndof',ndof, ...
        'a0',a0,'maxNeighborRatio',expectedMaxRatio, ...
        'loadType',C.load.type,'sig0',C.load.sig0, ...
        'anchorMode',C.bc.anchor_mode,'cornerNodes',[iLB iRB], ...
        'solver','pcg_free_dof_spd_sgs', ...
        'pcgTol',pcgTol,'pcgMaxIt',pcgMaxIt, ...
        'preconditioner',preconditionerName, ...
        'physicalEDIperformedAfterCheckpoint',true, ...
        'singleMatchedEDI',true, ...
        'noRadiusSweep',true, ...
        'noLevel0PhysicalSolve',true);

    [folder,~,~]=fileparts(cp);
    if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
    tmp=[cp '.incomplete.mat'];
    if exist(tmp,'file')==2
        error('step68:InterruptedSave', ...
            'Inspect/remove prior incomplete checkpoint manually: %s',tmp);
    end

    save(tmp,'mesh','U','mat','crack','a0','meta','solverInfo','C', ...
        'level1MeshSummary','-v7.3');
    [ok,msg]=movefile(tmp,cp);
    if ~ok,error('step68:CheckpointSave','%s',msg);end
    fprintf('  Level-1 physical field safely checkpointed: %s\n',cp);
    newSolve=true;

    clear Ap bp M1 M2 Dinv xp uf U d p
end

% -------------------------------------------------------------------------
% Phase 2: postprocess only from safely saved Level-1 field.
% -------------------------------------------------------------------------
fprintf('\nSTEP68 PHASE 2: LEVEL-1 COD + SINGLE MATCHED EDI\n');
s=load(cp,'mesh','U','mat','crack','a0','meta','solverInfo','level1MeshSummary');
if ~strcmp(s.meta.stage,'step68_level1_sgs_physical') || ...
        s.meta.nT3~=122691 || s.meta.nT6~=246701 || ...
        s.meta.meshLevel~=1 || abs(s.meta.scale-.5)>1e-14 || ...
        ~strcmp(s.meta.solver,'pcg_free_dof_spd_sgs')
    error('step68:CheckpointProvenance','Saved Level-1 checkpoint provenance mismatch.');
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

if ~isequal(fitTable.n_native,[74;74;108;108;86;86;67;67]) || ...
        any(~isfinite(fitTable.ratio_COD))
    error('step68:CODGate','Level-1 COD sampling/results invalid.');
end

CODconvergence=table( ...
    fitTable.lower_r_over_a0,fitTable.upper_r_over_a0,fitTable.degree, ...
    fitTable.n_native,level0FitRatio,fitTable.ratio_COD, ...
    fitTable.ratio_COD-level0FitRatio, ...
    100*(fitTable.ratio_COD-level0FitRatio)./level0FitRatio, ...
    'VariableNames',{'lower_r_over_a0','upper_r_over_a0','degree','level1_n_native', ...
    'level0_ratio','level1_ratio','delta_L1_minus_L0','relative_change_pct'});

% Exactly one matched physical EDI.
ri=.0008; ro=.65*s.a0;
if abs(ro-.0052)>1e-14
    error('step68:EDIRadius','Expected r_outer=5.2 mm.');
end
matEDI=s.mat;
if ~isfield(matEDI,'Dmat')&&isfield(matEDI,'D'),matEDI.Dmat=matEDI.D;end
fprintf('  Running ONE matched 16-point FE-nodal-q EDI on saved Level-1 field.\n');
[KI1,KII1]=SIF_LEFM_interaction_EDI( ...
    s.mesh,s.U,s.crack.Pmid,matEDI, ...
    struct('r_inner',ri,'r_outer',ro), ...
    'UsePlaneStrain',matEDI.ps==1, ...
    'Verbose',false, ...
    'WeightFunction','fe_nodal', ...
    'QuadratureRule',16, ...
    'StoreGPDiagnostics',false);
q1=KII1/KI1;
if ~all(isfinite([KI1,KII1,q1])) || KI1<=0
    error('step68:InvalidEDI','Level-1 matched EDI returned invalid values.');
end

EDI=table(KI1,KII1,q1,ri,ro,ro/s.a0, ...
    'VariableNames',{'KI','KII','ratio','r_inner','r_outer','r_outer_over_a0'});

EDIConvergence=table( ...
    level0EDI.KI,KI1,100*(KI1-level0EDI.KI)/level0EDI.KI, ...
    level0EDI.KII,KII1,100*(KII1-level0EDI.KII)/level0EDI.KII, ...
    level0EDI.ratio,q1,100*(q1-level0EDI.ratio)/level0EDI.ratio, ...
    'VariableNames',{'level0_KI','level1_KI','KI_change_pct', ...
    'level0_KII','level1_KII','KII_change_pct', ...
    'level0_ratio','level1_ratio','ratio_change_pct'});

% Cross-extractor consistency on Level 1.
codMin=min(fitTable.ratio_COD);
codMax=max(fitTable.ratio_COD);
codMean=mean(fitTable.ratio_COD);
codMedian=median(fitTable.ratio_COD);
CrossExtractor=table( ...
    q1,codMin,codMax,codMean,codMedian, ...
    q1>=codMin&&q1<=codMax, ...
    100*(codMean-q1)/q1,100*(codMedian-q1)/q1, ...
    'VariableNames',{'EDI_ratio','COD_min','COD_max','COD_mean','COD_median', ...
    'EDI_inside_COD_range','COD_mean_gap_pct','COD_median_gap_pct'});

% Numerical validity gates only. We deliberately do not gate the unknown
% physical Level-1 convergence outcome against an expected answer.
gates=struct();
gates.level1MeshExact=s.meta.nT3==122691 && s.meta.nT6==246701;
gates.pcgConverged=s.solverInfo.flag==0;
gates.pcgReportedResidual=s.solverInfo.relres<=pcgTol;
gates.trueResidual=s.solverInfo.trueRelResidual<=5e-10;
gates.constraintsExact=s.solverInfo.constraintInf<=1e-14;
gates.nativeSamplingExact=isequal(fitTable.n_native,[74;74;108;108;86;86;67;67]);
gates.codFinite=all(isfinite(fitTable.ratio_COD));
gates.ediFinite=all(isfinite([KI1,KII1,q1]))&&KI1>0;
gates.singleMatchedEDI=logical(s.meta.singleMatchedEDI);
gates.noRadiusSweep=logical(s.meta.noRadiusSweep);
gates.noDirectBackslash=logical(s.solverInfo.noDirectBackslash);
gates.noSolverTuning=logical(s.solverInfo.noSolverTuning);
gates.noLevel0PhysicalSolve=logical(s.meta.noLevel0PhysicalSolve);

numericalPass=all(structfun(@logical,gates));

fprintf('\nLEVEL-0 -> LEVEL-1 COD CONVERGENCE\n');
disp(CODconvergence);
fprintf('\nLEVEL-0 -> LEVEL-1 MATCHED EDI CONVERGENCE\n');
disp(EDIConvergence);
fprintf('\nLEVEL-1 CROSS-EXTRACTOR CONSISTENCY\n');
disp(CrossExtractor);
fprintf('\nLEVEL-1 ITERATIVE SOLVER INFO\n');
disp(s.solverInfo);
fprintf('\nSTEP68 NUMERICAL GATES\n');
disp(gates);
fprintf('  Level-0 reference source: %s\n',level0ReferenceSource);
fprintf('  Level-0 reference precision: %s\n',level0ReferencePrecision);

if numericalPass
    fprintf('STEP68 COMPLETE: Level-1 physical field and fixed postprocessing are numerically valid.\n');
    fprintf('  Interpret convergence from the reported Level-0 -> Level-1 changes.\n');
else
    fprintf('STEP68 STOP: numerical validity gate failed. Do not interpret convergence.\n');
end

R68=struct( ...
    'checkpointPath',cp, ...
    'newSolve',newSolve, ...
    'solverInfo',s.solverInfo, ...
    'level1MeshSummary',s.level1MeshSummary, ...
    'fitTable',fitTable, ...
    'CODconvergence',CODconvergence, ...
    'EDI',EDI, ...
    'EDIConvergence',EDIConvergence, ...
    'CrossExtractor',CrossExtractor, ...
    'level0ReferenceSource',level0ReferenceSource, ...
    'level0ReferencePrecision',level0ReferencePrecision, ...
    'gates',gates, ...
    'numericalPass',numericalPass, ...
    'oneLevel1PhysicalSolve',true, ...
    'singleMatchedEDI',true, ...
    'noRadiusSweep',true, ...
    'noSolverTuning',true, ...
    'interpretation',['First physical Level-0/Level-1 convergence comparison ', ...
      'within the same qualified C03 structured mesh family using the ', ...
      'Level-0-qualified SGS-PCG solver. Physical convergence is reported, ', ...
      'not forced by an a priori answer gate.']);

save(saveFile,'R68','-v7');
fprintf('  Compact Step68 result saved: %s\n',saveFile);
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

function [availGiB,totalGiB]=available_memory_gib()
availGiB=NaN;totalGiB=NaN;
if ~ispc,return,end
GiB=1024^3;
try
    [u,s]=memory; %#ok<ASGLU>
    if isfield(s,'PhysicalMemory')
        pm=s.PhysicalMemory;
        if isfield(pm,'Available')&&isfinite(pm.Available)
            availGiB=double(pm.Available)/GiB;
        end
        if isfield(pm,'Total')&&isfinite(pm.Total)
            totalGiB=double(pm.Total)/GiB;
        end
    end
    if isnan(availGiB)&&isfield(u,'MemAvailableAllArrays')&& ...
            isfinite(u.MemAvailableAllArrays)
        availGiB=double(u.MemAvailableAllArrays)/GiB;
    end
catch
end
end

function assert_step68_branch(root)
[status,b]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0&&strcmp(strtrim(b),'audit/step68-level1-sgs-physical-convergence'), ...
    'step68:Branch', ...
    'Run Step68 only on audit/step68-level1-sgs-physical-convergence.');
end
