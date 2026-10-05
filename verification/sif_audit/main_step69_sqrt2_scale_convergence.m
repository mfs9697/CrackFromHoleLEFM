function R69=main_step69_sqrt2_scale_convergence(varargin)
%MAIN_STEP69_SQRT2_SCALE_CONVERGENCE
% Additional physical mesh at s=1/sqrt(2), between s=1 and s=1/2.
%
% Purpose:
%   create a geometrically spaced three-mesh family
%       s0 = 1, sm = 1/sqrt(2), s1 = 1/2
%   so an observed order p can be estimated where the three physical
%   results are monotone.
%
% DEFAULT IS SAFE:
%   AllowSolve = false
%
% The run:
%   1) reconstructs/qualifies the continuous-scale mesh, including
%      prescribed pure-I/pure-II/mixed controls;
%   2) performs exactly one SGS-PCG physical solve if authorized;
%   3) checkpoints the field immediately;
%   4) computes the same eight COD fits and one matched 0.8-5.2 mm EDI;
%   5) forms three-level convergence estimates.
%
% No radius sweep, no solver tuning, no other physical mesh is solved.
%
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');

ip=inputParser;
addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'SourceCandidateFile', ...
    fullfile(vdir,'step62_structured_graded_mesh_candidate_T3.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'CheckpointFile', ...
    fullfile(vdir,'step69_sqrt2_scale_sgs_physical_solved.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'SaveFile', ...
    fullfile(vdir,'step69_sqrt2_scale_convergence_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'Level0ReferenceFile', ...
    fullfile(vdir,'step67a_level0_sgs_qualification_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'Level1ReferenceFile', ...
    fullfile(vdir,'step68_level1_sgs_convergence_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;

addpath(genpath(root));
assert_step69_branch(root);

sourceFile=char(opt.SourceCandidateFile);
cp=char(opt.CheckpointFile);
saveFile=char(opt.SaveFile);
level0File=char(opt.Level0ReferenceFile);
level1File=char(opt.Level1ReferenceFile);

if exist(sourceFile,'file')~=2
    error('step69:MissingSource','Missing archived Step62 source: %s',sourceFile);
end
for name={'pcg','symamd'}
    if ~matlab_callable(name{1})
        error('step69:MissingTool','Step69 requires callable MATLAB %s.',name{1});
    end
end

s0=1;
sm=1/sqrt(2);
s1=.5;
r=sqrt(2);

fprintf('\n============================================================\n');
fprintf('STEP 69: THREE-SCALE CONVERGENCE, MIDPOINT s=1/sqrt(2)\n');
fprintf('============================================================\n');
fprintf('  Existing physical scales: s0=1, s1=0.5.\n');
fprintf('  Additional scale: sm=1/sqrt(2)=%.12g.\n',sm);
fprintf('  Constant refinement ratio r=sqrt(2)=%.12g.\n',r);
fprintf('  One physical midpoint solve only; no radius sweep/tuning.\n');

% -------------------------------------------------------------------------
% Mesh qualification at explicit continuous scale sm.
% -------------------------------------------------------------------------
cal=c03_calibration();
Om=main_step62_structured_graded_mesh( ...
    'SourceCandidateFile',sourceFile, ...
    'SavePrefix',fullfile(vdir,'step69_sqrt2_mesh'), ...
    'Visible','off', ...
    'RunSynthetic',true, ...
    'ExteriorCalibration',cal, ...
    'WriteArtifacts',false, ...
    'ReturnCandidate',true, ...
    'Verbose',true, ...
    'Level',0, ...
    'Scale',sm);

if ~Om.gates.structuralPass || ~Om.synthetic.performed || ~Om.synthetic.passed || ...
        abs(Om.design.scale-sm)>1e-14 || ...
        abs(Om.summary.pairedRadius_mm-6)>1e-12 || ...
        Om.maxNeighborSizeRatio>1.8+5e-12 || ...
        ~Om.gates.completeT3Pairing || ~Om.gates.completeT6Pairing || ...
        ~Om.gates.candidateSupportInsidePairedRegion || ...
        ~Om.gates.exteriorOutsideSupport || ...
        ~Om.gates.tipFanSixTriangles
    error('step69:MidMeshQualification', ...
        's=1/sqrt(2) mesh failed the predeclared family qualification.');
end

candidate=Om.candidate;
midMeshSummary=struct( ...
    'scale',Om.design.scale, ...
    'targetTip_mm',1e3*Om.targetTipEdge_m, ...
    'actualTipMedian_mm',Om.summary.newTipMedian_mm, ...
    'nT3',Om.summary.newT3, ...
    'nT6',Om.summary.newT6Nodes, ...
    'pairedRadius_mm',Om.summary.pairedRadius_mm, ...
    'patchMinAngle_deg',Om.summary.patchMinAngle_deg, ...
    'exteriorMinAngle_deg',Om.summary.exteriorMinAngle_deg, ...
    'maxNeighborRatio',Om.maxNeighborSizeRatio, ...
    'primarySupport',Om.summary.newPrimarySupport, ...
    'literalSupport',Om.candidateLiteralSupportCount, ...
    'samplingTable',Om.samplingTable, ...
    'syntheticTable',Om.synthetic.table);
clear Om

P=candidate.p;T=candidate.t;crack=candidate.crack;mat0=candidate.mat;
clear candidate
[P6,T6]=T3toT6_fast(P,T);
mesh=struct('coord3',P,'connect3',T,'coord',P6,'connect',T6);
a0=norm(diff(crack.Pmid,1,1));
if abs(a0-.008)>1e-12
    error('step69:CrackLength','Expected a0=8 mm.');
end

% Native sampling pre-solve.
zeroU=zeros(2*size(P6,1),1);
[rZero,~,face0]=native_COD_audit(mesh,zeroU,mat0,crack,8);
rrZero=rZero/a0;
windows=[.04 .20;.04 .30;.08 .30;.12 .30];
sampleN=zeros(4,1);
for k=1:4
    sampleN(k)=nnz(rrZero>=windows(k,1)&rrZero<=windows(k,2));
end
clear zeroU rZero rrZero
if face0.nUpper~=face0.nLower || face0.gridMismatch>1e-12 || any(sampleN<12)
    error('step69:NativeSampling','Midpoint native crack-face sampling invalid.');
end

fprintf('\nMIDPOINT MESH QUALIFICATION\n');
fprintf('  scale=%.12g; target tip=%.9g mm; actual tip=%.9g mm.\n', ...
    sm,midMeshSummary.targetTip_mm,midMeshSummary.actualTipMedian_mm);
fprintf('  T3=%d, T6=%d, max neighbor ratio=%.9g.\n', ...
    midMeshSummary.nT3,midMeshSummary.nT6,midMeshSummary.maxNeighborRatio);
fprintf('  Native COD windows: [%s].\n',sprintf('%d ',sampleN));
fprintf('  Prescribed-field qualification: PASS.\n');

% -------------------------------------------------------------------------
% Same physical setup and Step67A-qualified SGS-PCG solver.
% -------------------------------------------------------------------------
xmin=min(P(:,1));xmax=max(P(:,1));ymin=min(P(:,2));ymax=max(P(:,2));
A=xmax;B=max(abs([ymin,ymax]));
if abs(xmin)>1e-12 || abs(A-.30)>1e-12 || ...
        abs(ymin+.10)>1e-12 || abs(ymax-.10)>1e-12
    error('step69:PlateBoundary','Physical plate boundary changed.');
end

Cref=cfg_hole_initiation();
if abs(Cref.E-mat0.E)>1e-12*max(1,abs(mat0.E)) || ...
        abs(Cref.nu-mat0.nu)>1e-14 || Cref.ps~=mat0.ps || ...
        ~strcmp(Cref.load.type,'remote_tension_y') || ...
        abs(Cref.load.sig0-1)>1e-14 || ...
        ~strcmp(Cref.bc.anchor_mode,'minimal')
    error('step69:PhysicalConfig','Physical configuration changed.');
end

mat=material_with_D(mat0);
[nip2,xip2,w2,Nextr]=integr();
quad=struct('nip2',nip2,'xip2',xip2,'w2',w2,'Nextr',Nextr);

[~,iLB]=min(sum((P-[0,-B]).^2,2));
[~,iRB]=min(sum((P-[A,-B]).^2,2));
fixvar=unique([2*iLB-1;2*iLB;2*iRB]);
ndof=2*size(P6,1);
freeMask=true(ndof,1);freeMask(fixvar)=false;free=find(freeMask);
clear freeMask

C=struct('A',A,'B',B,'E',mat.E,'nu',mat.nu,'ps',mat.ps,'a0',a0, ...
    'load',Cref.load,'bc',Cref.bc);

pcgTol=1e-10;
pcgMaxIt=5000;
preconditionerName='symmetric_gauss_seidel';

% -------------------------------------------------------------------------
% Endpoint references.
% COD fallbacks are high precision at s=1 and reconstructed from the
% Step68 printed deltas at s=0.5. EDI fallbacks are displayed/rounded.
% Exact saved endpoint files override these fallbacks when available.
% -------------------------------------------------------------------------
qCOD0=[ ...
    0.0001052647913;
    0.0001056239885;
    0.0001048946546;
    0.0001057078877;
    0.0001045631315;
    0.0001058187639;
    0.0001041556245;
    0.0001058252201];
qCOD1=qCOD0+[ ...
    4.3927e-7;
    5.5086e-7;
    3.9598e-7;
    4.8944e-7;
    3.3970e-7;
    4.0703e-7;
    3.3389e-7;
    3.7156e-7];

EDI0=struct('KI',0.43785,'KII',4.6547e-5,'ratio',0.00010631);
EDI1=struct('KI',0.43786,'KII',4.6643e-5,'ratio',0.00010652);
codEndpointSource="embedded_audit_fingerprints";
ediEndpointSource="displayed_rounded_audit_fingerprints";
ediFormalOrder=false;

if exist(level0File,'file')==2
    z=load(level0File,'R67a');
    if isfield(z,'R67a')&&z.R67a.qualificationPass&&height(z.R67a.fitTable)==8
        qCOD0=z.R67a.fitTable.ratio_COD;
        EDI0=struct('KI',z.R67a.EDIcomparison.iterative_KI, ...
            'KII',z.R67a.EDIcomparison.iterative_KII, ...
            'ratio',z.R67a.EDIcomparison.iterative_ratio);
        codEndpointSource="saved_step67a_plus_audit_step68";
        ediEndpointSource="saved_step67a_plus_rounded_step68";
    end
end
if exist(level1File,'file')==2
    z=load(level1File,'R68');
    if isfield(z,'R68')&&z.R68.numericalPass&&height(z.R68.fitTable)==8
        qCOD1=z.R68.fitTable.ratio_COD;
        EDI1=struct('KI',z.R68.EDI.KI,'KII',z.R68.EDI.KII, ...
            'ratio',z.R68.EDI.ratio);
        if contains(codEndpointSource,"saved_step67a")
            codEndpointSource="exact_saved_step67a_and_step68";
        else
            codEndpointSource="audit_step67a_plus_saved_step68";
        end
        if exist(level0File,'file')==2
            ediEndpointSource="exact_saved_step67a_and_step68";
            ediFormalOrder=true;
        else
            ediEndpointSource="rounded_step67a_plus_saved_step68";
        end
    end
end

% -------------------------------------------------------------------------
% One midpoint physical solve, checkpoint first.
% -------------------------------------------------------------------------
newSolve=false;
if exist(cp,'file')==2
    sOld=load(cp,'meta','mesh','U','mat','crack','a0','solverInfo');
    if ~isfield(sOld,'meta') || ...
            ~strcmp(sOld.meta.stage,'step69_sqrt2_scale_sgs_physical') || ...
            abs(sOld.meta.scale-sm)>1e-14 || ...
            sOld.meta.nT3~=size(T,1) || sOld.meta.nT6~=size(P6,1) || ...
            abs(sOld.a0-a0)>1e-12 || ...
            ~isequal(sOld.mesh.connect3,T) || ...
            max(abs(sOld.mesh.coord3(:)-P(:)))>1e-12 || ...
            numel(sOld.U)~=ndof || any(~isfinite(sOld.U))
        error('step69:ExistingCheckpointMismatch', ...
            'Existing Step69 checkpoint belongs to another field.');
    end
    fprintf('  Reusing existing midpoint checkpoint; NO new physical solve.\n');
else
    if ~opt.AllowSolve
        error('step69:ExplicitSolveApprovalRequired', ...
            ['Midpoint mesh passed qualification. Physical solve is guarded; ', ...
             'rerun with ''AllowSolve'',true only after authorization.']);
    end

    fprintf('\nSTEP69 PHASE 1: MIDPOINT SGS-PCG PHYSICAL SOLVE\n');
    [availBefore,totalPhysical]=available_memory_gib();
    if isfinite(availBefore)
        fprintf('  Available physical memory before assembly: %.3f GiB.\n',availBefore);
    end

    K=stif_assem(mesh,mat,quad,[]);
    symErr=norm(K-K.',1)/max(1,norm(K,1));
    if ~isfinite(symErr)||symErr>5e-13
        error('step69:Symmetry','Stiffness symmetry gate failed: %.3e.',symErr);
    end

    F=zeros(ndof,1);
    eps1=1e-8*max(A,2*B);
    elod=edge_loads_T6(mesh.coord,B,eps1);
    for k=1:size(elod,1)
        F(2*elod(k,1))=F(2*elod(k,1))+C.load.sig0*elod(k,2);
    end
    F(fixvar)=0;

    Kff=K(free,free);
    Ff=F(free);
    Kff=(Kff+Kff.')/2;
    if any(diag(Kff)<=0)
        error('step69:Diagonal','Free-DOF K has nonpositive diagonal.');
    end

    p=symamd(Kff);
    Ap=Kff(p,p);
    bp=Ff(p);
    nnzK=nnz(K);nnzAp=nnz(Ap);
    clear K Kff F Ff

    d=diag(Ap);
    if any(~isfinite(d))||any(d<=0)
        error('step69:SGSDiagonal','SGS diagonal invalid.');
    end
    tPre=tic;
    M1=tril(Ap);
    Dinv=spdiags(1./d,0,numel(d),numel(d));
    M2=Dinv*M1.';
    preconditionerSeconds=toc(tPre);
    nnzM1=nnz(M1);nnzM2=nnz(M2);
    preconditionerGiB=(16*(nnzM1+nnzM2)+ ...
        8*(size(M1,2)+size(M2,2)+2))/1024^3;
    [availAfterPre,totalPhysical2]=available_memory_gib();
    if isfinite(totalPhysical2),totalPhysical=totalPhysical2;end

    fprintf('  Symmetry error %.3e; SGS storage %.4f GiB.\n',symErr,preconditionerGiB);
    fprintf('  Starting exactly ONE authorized s=1/sqrt(2) PCG solve.\n');

    tSolve=tic;
    [xp,flag,relres,iter,resvec]=pcg(Ap,bp,pcgTol,pcgMaxIt,M1,M2);
    solveSeconds=toc(tSolve);
    if flag~=0||~isfinite(relres)||relres>pcgTol||any(~isfinite(xp))
        error('step69:PCGFailed', ...
            'Fixed SGS-PCG failed: flag=%d, relres=%.3e, iter=%d.',flag,relres,iter);
    end

    trueRelResidual=norm(Ap*xp-bp)/max(norm(bp),eps);
    if trueRelResidual>5e-10
        error('step69:Residual','True residual %.3e exceeds gate.',trueRelResidual);
    end

    uf=zeros(numel(free),1);uf(p)=xp;
    U=zeros(ndof,1);U(free)=uf;U(fixvar)=0;
    constraintInf=max(abs(U(fixvar)));
    if constraintInf>1e-14
        error('step69:Constraint','Constraint gate failed.');
    end

    solverInfo=struct( ...
        'method','pcg_free_dof_spd_sgs', ...
        'scale',sm,'pcgTol',pcgTol,'pcgMaxIt',pcgMaxIt, ...
        'preconditioner',preconditionerName, ...
        'flag',flag,'relres',relres,'iter',iter, ...
        'trueRelResidual',trueRelResidual,'constraintInf',constraintInf, ...
        'resvecFinal',resvec(end),'resvecLength',numel(resvec), ...
        'symmetryError',symErr,'nnzK',nnzK,'nnzAp',nnzAp, ...
        'nnzM1',nnzM1,'nnzM2',nnzM2, ...
        'preconditionerGiBApprox',preconditionerGiB, ...
        'preconditionerSeconds',preconditionerSeconds, ...
        'solveSeconds',solveSeconds, ...
        'availableGiBBeforeAssembly',availBefore, ...
        'availableGiBAfterPreconditioner',availAfterPre, ...
        'totalPhysicalGiB',totalPhysical, ...
        'onePhysicalLinearSolve',true,'noDirectBackslash',true, ...
        'noSolverTuning',true);

    meta=struct( ...
        'stage','step69_sqrt2_scale_sgs_physical', ...
        'scale',sm,'refinementRatio',r, ...
        'nT3',size(T,1),'nT6',size(P6,1),'ndof',ndof, ...
        'a0',a0,'maxNeighborRatio',midMeshSummary.maxNeighborRatio, ...
        'solver','pcg_free_dof_spd_sgs','pcgTol',pcgTol, ...
        'singleMatchedEDI',true,'noRadiusSweep',true,'noSolverTuning',true);

    [folder,~,~]=fileparts(cp);
    if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
    tmp=[cp '.incomplete.mat'];
    if exist(tmp,'file')==2
        error('step69:InterruptedSave','Inspect/remove incomplete file: %s',tmp);
    end
    save(tmp,'mesh','U','mat','crack','a0','meta','solverInfo','C', ...
        'midMeshSummary','-v7.3');
    [ok,msg]=movefile(tmp,cp);
    if ~ok,error('step69:CheckpointSave','%s',msg);end
    fprintf('  Midpoint physical field safely checkpointed: %s\n',cp);
    newSolve=true;
    clear Ap bp M1 M2 Dinv xp uf U d p
end

% -------------------------------------------------------------------------
% Postprocess from saved midpoint field.
% -------------------------------------------------------------------------
fprintf('\nSTEP69 PHASE 2: MIDPOINT COD + SINGLE MATCHED EDI\n');
s=load(cp,'mesh','U','mat','crack','a0','meta','solverInfo','midMeshSummary');
if ~strcmp(s.meta.stage,'step69_sqrt2_scale_sgs_physical') || ...
        abs(s.meta.scale-sm)>1e-14
    error('step69:CheckpointProvenance','Midpoint checkpoint provenance mismatch.');
end

[rCOD,app,face]=native_COD_audit(s.mesh,s.U,s.mat,s.crack,8);
rr=rCOD/s.a0;
fitRows=nan(8,9);n=0;qRaw=app(:,2)./app(:,1);
for iw=1:4
    ids=find(rr>=windows(iw,1)&rr<=windows(iw,2));
    for degree=[1 2]
        pI=polyfit(rr(ids),app(ids,1),degree);
        pII=polyfit(rr(ids),app(ids,2),degree);
        predII=polyval(pII,rr(ids));
        n=n+1;
        fitRows(n,:)=[windows(iw,:),degree,numel(ids), ...
            pI(end),pII(end),pII(end)/pI(end), ...
            sqrt(mean((app(ids,2)-predII).^2)),median(qRaw(ids))];
    end
end
fitTable=array2table(fitRows,'VariableNames', ...
    {'lower_r_over_a0','upper_r_over_a0','degree','n_native', ...
     'KI_COD','KII_COD','ratio_COD','RMSE_KII','median_pointwise_ratio'});

ri=.0008;ro=.65*s.a0;
matEDI=s.mat;if ~isfield(matEDI,'Dmat')&&isfield(matEDI,'D'),matEDI.Dmat=matEDI.D;end
fprintf('  Running ONE matched 16-point FE-nodal-q EDI on midpoint field.\n');
[KIm,KIIm]=SIF_LEFM_interaction_EDI( ...
    s.mesh,s.U,s.crack.Pmid,matEDI, ...
    struct('r_inner',ri,'r_outer',ro), ...
    'UsePlaneStrain',matEDI.ps==1,'Verbose',false, ...
    'WeightFunction','fe_nodal','QuadratureRule',16,'StoreGPDiagnostics',false);
qmEDI=KIIm/KIm;
EDI=table(KIm,KIIm,qmEDI,ri,ro, ...
    'VariableNames',{'KI','KII','ratio','r_inner','r_outer'});

% -------------------------------------------------------------------------
% Three-level convergence.
% -------------------------------------------------------------------------
nfit=height(fitTable);
pCOD=nan(nfit,1);qInfCOD=nan(nfit,1);orderOK=false(nfit,1);
for k=1:nfit
    [pCOD(k),qInfCOD(k),orderOK(k)]=three_level_order( ...
        qCOD0(k),fitTable.ratio_COD(k),qCOD1(k),r);
end

ThreeScaleCOD=table( ...
    fitTable.lower_r_over_a0,fitTable.upper_r_over_a0,fitTable.degree, ...
    repmat(s0,nfit,1),repmat(sm,nfit,1),repmat(s1,nfit,1), ...
    qCOD0,fitTable.ratio_COD,qCOD1, ...
    100*(fitTable.ratio_COD-qCOD0)./qCOD0, ...
    100*(qCOD1-fitTable.ratio_COD)./fitTable.ratio_COD, ...
    pCOD,qInfCOD,orderOK, ...
    'VariableNames',{'lower_r_over_a0','upper_r_over_a0','degree', ...
    's_coarse','s_mid','s_fine','q_coarse','q_mid','q_fine', ...
    'coarse_to_mid_pct','mid_to_fine_pct','observedOrder_p', ...
    'richardson_q_infinity','orderAdmissible'});

[pEDI,qInfEDI,ediOrderAdmissible]=three_level_order( ...
    EDI0.ratio,qmEDI,EDI1.ratio,r);
ThreeScaleEDI=table( ...
    s0,sm,s1,EDI0.KI,KIm,EDI1.KI, ...
    EDI0.KII,KIIm,EDI1.KII, ...
    EDI0.ratio,qmEDI,EDI1.ratio, ...
    100*(qmEDI-EDI0.ratio)/EDI0.ratio, ...
    100*(EDI1.ratio-qmEDI)/qmEDI, ...
    pEDI,qInfEDI,ediOrderAdmissible,ediFormalOrder, ...
    'VariableNames',{'s_coarse','s_mid','s_fine', ...
    'KI_coarse','KI_mid','KI_fine','KII_coarse','KII_mid','KII_fine', ...
    'q_coarse','q_mid','q_fine','coarse_to_mid_pct','mid_to_fine_pct', ...
    'observedOrder_p','richardson_q_infinity','orderAdmissible', ...
    'endpointReferencesExact'});

% Numerical integrity gates only.
gates=struct();
gates.midMeshQualified=true;
gates.syntheticPassed=true;
gates.pcgConverged=s.solverInfo.flag==0;
gates.trueResidual=s.solverInfo.trueRelResidual<=5e-10;
gates.constraintsExact=s.solverInfo.constraintInf<=1e-14;
gates.codFinite=all(isfinite(fitTable.ratio_COD));
gates.ediFinite=all(isfinite([KIm,KIIm,qmEDI]))&&KIm>0;
gates.singleMatchedEDI=true;
gates.noRadiusSweep=true;
gates.noSolverTuning=true;
numericalPass=all(structfun(@logical,gates));

fprintf('\nTHREE-SCALE COD CONVERGENCE\n');
disp(ThreeScaleCOD);
fprintf('\nTHREE-SCALE MATCHED EDI CONVERGENCE\n');
disp(ThreeScaleEDI);
fprintf('\nMIDPOINT SOLVER INFO\n');
disp(s.solverInfo);
fprintf('\nSTEP69 NUMERICAL GATES\n');
disp(gates);
fprintf('  COD endpoint source: %s\n',codEndpointSource);
fprintf('  EDI endpoint source: %s\n',ediEndpointSource);
if ~ediFormalOrder
    fprintf(['  NOTE: EDI p/Richardson values use at least one rounded endpoint ', ...
        'reference and are approximate, not formal.\n']);
end

R69=struct( ...
    'scale',sm,'refinementRatio',r, ...
    'checkpointPath',cp,'newSolve',newSolve, ...
    'midMeshSummary',s.midMeshSummary, ...
    'solverInfo',s.solverInfo,'fitTable',fitTable,'EDI',EDI, ...
    'ThreeScaleCOD',ThreeScaleCOD,'ThreeScaleEDI',ThreeScaleEDI, ...
    'codEndpointSource',codEndpointSource, ...
    'ediEndpointSource',ediEndpointSource, ...
    'ediFormalOrder',ediFormalOrder, ...
    'gates',gates,'numericalPass',numericalPass, ...
    'interpretation',['Three geometrically spaced mesh scales. Observed ', ...
      'orders are reported only where consecutive changes have the same ', ...
      'sign and decrease consistently; EDI order is marked approximate ', ...
      'unless exact saved endpoint references are available.']);

save(saveFile,'R69','-v7');
fprintf('  Compact Step69 result saved: %s\n',saveFile);
end

% =========================================================================
function [p,qInf,ok]=three_level_order(q0,qm,q1,r)
p=NaN;qInf=NaN;ok=false;
d0=q0-qm;d1=qm-q1;
scale=max(abs([q0,qm,q1]));
tol=max(1e-15,1e-10*scale);
if ~all(isfinite([q0,qm,q1,r]))||r<=1||abs(d0)<=tol||abs(d1)<=tol
    return
end
ratio=d0/d1;
if ratio<=1
    return
end
p=log(ratio)/log(r);
if ~isfinite(p)||p<=0
    p=NaN;return
end
den=r^p-1;
if abs(den)<=1e-12,return,end
qInf=q1+(q1-qm)/den;
ok=isfinite(qInf);
end

function mat=material_with_D(mat0)
mat=mat0;E=mat.E;nu=mat.nu;
if mat.ps==1
    coef=E/((1+nu)*(1-2*nu));
    D=coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
else
    coef=E/(1-nu^2);
    D=coef*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
end
mat.D=D;mat.Dmat=D;
end

function cal=c03_calibration()
cal=struct('transitionLength_m',.008,'farSlope',.10, ...
    'boundaryMetricGrowth',.25,'smoothingSteps',6, ...
    'refinementMaxPasses',140,'refinementMinAngle_deg',25, ...
    'refinementLongestFactor',1.65,'neighborRatioTarget',1.8, ...
    'verbose',false);
end

function tf=matlab_callable(name)
tf=exist(name,'file')~=0||exist(name,'builtin')~=0;
end

function [availGiB,totalGiB]=available_memory_gib()
availGiB=NaN;totalGiB=NaN;if ~ispc,return,end
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
    if isnan(availGiB)&&isfield(u,'MemAvailableAllArrays')&&isfinite(u.MemAvailableAllArrays)
        availGiB=double(u.MemAvailableAllArrays)/GiB;
    end
catch
end
end

function assert_step69_branch(root)
[status,b]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0&&strcmp(strtrim(b),'audit/step69-sqrt2-scale-convergence'), ...
    'step69:Branch', ...
    'Run Step69 only on audit/step69-sqrt2-scale-convergence.');
end
