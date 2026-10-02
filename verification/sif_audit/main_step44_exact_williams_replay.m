function O44=main_step44_exact_williams_replay(O38,varargin)
%MAIN_STEP44_EXACT_WILLIAMS_REPLAY
% Exact Williams-field COD test on SAVED Step-38 T6 geometry.
% Zero new meshes, FEM solutions, stiffness matrices, or real-field EDI.
% Cases: unit I, unit II, and strongly Mode-I-dominated mixed mode.
% The COD algorithm is EXACTLY shared with the Step-38 postprocessor.
% First require agreement with the previously stored real-field COD arrays.
% Optional RunEDI=true: TWO costly 16-point FE-nodal-q integrations on
% the exact unit fields at ONE prescribed original annulus. Unit-case
% results reconstruct the mixed case by linearity (not a third integral).
% Individual EDI cases are cached immediately to survive interruptions.
%
% Usage:
%   load('step38_tip_refined_solved_results.mat','O38')
%   O44=main_step44_exact_williams_replay(O38);
%   % Only after reviewing COD results:
%   O44=main_step44_exact_williams_replay(O38,'RunEDI',true);
%
% Small output: profiles, fits, numerical gates; no mesh, U or stiffness.
% These exact-field tests verify extraction, NOT the accuracy of FEM U.

ip=inputParser;
addParameter(ip,'CheckpointPath',O38.checkpointPath, ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'RunEDI',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'EDIRouterOverA0',0.65, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'FitWindows',[.04 .30;.08 .30;.12 .30], ...
    @(x)isnumeric(x)&&size(x,2)==2&&all(isfinite(x(:)))&& ...
        all(x(:,1)>0)&all(x(:,1)<x(:,2))&all(x(:,2)<.7));
addParameter(ip,'FitDegrees',[1 2], ...
    @(x)isnumeric(x)&&isvector(x)&&all(ismember(x,[1 2])));
addParameter(ip,'SavePrefix',fullfile(pwd,'step44_exact_williams'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'MaxPureLeakageOverKI',1e-8, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'MaxMixedKIIErrorRel',1e-4, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
parse(ip,varargin{:});
opt=ip.Results;
need(O38,'nativeR');need(O38,'nativeApparent');
need(O38,'rInner');need(O38,'rOuterOverA0');
need(O38,'KI');need(O38,'KII');
cp=char(opt.CheckpointPath);
if exist(cp,'file')~=2
    error('step44:CheckpointAbsent','Solved Step-38 checkpoint missing: %s',cp);
end
d=dir(cp);
checkpointKey=sprintf('%s|%d|%.16g',cp,d.bytes,d.datenum);
prefix=char(opt.SavePrefix);
[folder,~,~]=fileparts(prefix);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end

% Deliberately require the solved checkpoint, not the small O38 arrays
% alone. Compact arrays are insufficient to replay exact fields on T6.
s=load(cp,'mesh','U','mat','crack','a0');
for f={'mesh','U','mat','crack','a0'}
    need(s,f{1});
end
mesh=s.mesh;U=s.U;mat=s.mat;crack=s.crack;a0=s.a0;
clear s
n=size(mesh.coord,1);
if numel(U)~=2*n || any(~isfinite(U)) || ...
        abs(norm(diff(crack.Pmid,1,1))/a0-1)>1e-10
    error('step44:CheckpointMismatch','Saved T6, U or crack is inconsistent.');
end
fprintf('\n============================================================\n');
fprintf('STEP 44: EXACT WILLIAMS REPLAY ON SAVED STEP-38 MESH\n');
fprintf('============================================================\n');
fprintf('  T6 nodes=%d, crack length=%.8g m. No new FEM solve.\n',n,a0);

% Regression gate: the extracted shared routine MUST reproduce the
% original Step-38 native physical COD before synthetic testing.
[rActual,appActual,face]=native_COD_audit(mesh,U,mat,crack,8,true);
rOld=O38.nativeR(:);
appOld=O38.nativeApparent;
if ~isequal(size(rActual),size(rOld)) || ...
        ~isequal(size(appActual),size(appOld))
    error('step44:RefactorCount','Step-38 COD count/shape changed.');
end
dr=max(abs(rActual-rOld));
dk=max(abs(appActual(:)-appOld(:)));
passRegression=dr<=1e-12*a0 && dk<=1e-12*max(1,max(abs(appOld(:))));
if ~passRegression
    error('step44:RefactorRegression', ...
        'Shared COD differs from stored Step38: dr=%g, dK=%g.',dr,dk);
end
fprintf('  REGRESSION PASS: n=%d; max dr=%g m, max dK=%g\n', ...
    numel(rActual),dr,dk);
faceSide=face.faceSide;
clear U appActual face
if ~isfield(mat,'ps') || ~isfield(mat,'E') || ~isfield(mat,'nu')
    error('step44:Material','Checkpoint missing isotropic elasticity fields.');
end

% Same local axes as the Step-38 native COD extractor.
tip=crack.Pmid(end,:);
e1=(tip-crack.Pmid(1,:))/a0;
e2=[-e1(2),e1(1)];
R=[e1(:),e2(:)];
xl=(mesh.coord-tip)*R;
upper=find(faceSide==+1);
lower=find(faceSide==-1);
faceIDs=find(faceSide~=0);
if isempty(upper)||isempty(lower) || ...
        any(ismember(upper,lower))
    error('step44:FaceIDs','Crack-face labels invalid.');
end
% Real nodal face geometry was already accepted by the SAME Step-38
% on-face tolerance. For the ideal analytical branch cut, snap ONLY
% negligible transverse roundoff in evaluation coordinates (not mesh).
faceResidual=max(abs(xl(faceIDs,2)));
faceTol=max(1e-12,1e-8*a0);
if faceResidual>faceTol || any(xl(faceIDs,1)>=-faceTol)
    error('step44:FaceGeometry', ...
        'Synthetic evaluation requires actual classified crack-face nodes.');
end
xFace=xl(faceIDs,:);
xFace(:,2)=0;
fprintf('  max transverse face-coordinate roundoff=%.3e m\n',faceResidual);
% No real-field SIF is assumed correct: tiny mixed-mode amplitudes
% merely define a known synthetic signal near the regime of interest.
rat=O38.rOuterOverA0(:);
[~,iref]=min(abs(rat-opt.EDIRouterOverA0));
testK=[1 0;0 1;O38.KI(iref) O38.KII(iref)];
if ~all(isfinite(testK(:))) || testK(3,1)<=0 || testK(3,2)==0
    error('step44:TestAmplitudes','Invalid stored small mixed-mode amplitudes.');
end
caseNames={'Exact_I','Exact_II','Exact_tiny_mixed'};
profile=cell(3,1);
radial=cell(3,1);
maxRawDeviation=zeros(3,2);
fits=nan(3*size(opt.FitWindows,1)*numel(opt.FitDegrees), ...
    11);
nr=0;
for icase=1:3
    % COD depends ONLY on displacements at crack-face T6 nodes. Keep
    % interior displacement identically zero to avoid unnecessary
    % expensive whole-mesh exact-field generation for these COD tests.
    ul=exact_williams_displacement_audit(xFace, ...
        testK(icase,1),testK(icase,2),mat.E,mat.nu,mat.ps, ...
        'UpperFaceIDs',find(faceSide(faceIDs)==1), ...
        'LowerFaceIDs',find(faceSide(faceIDs)==-1));
    uv=reshape(ul,2,[]).'*R.';
    synthetic=zeros(2*n,1);
    synthetic(2*faceIDs-1)=uv(:,1);
    synthetic(2*faceIDs)=uv(:,2);
    clear ul uv
    [r,p]=native_COD_audit(mesh,synthetic,mat,crack,8);
    clear synthetic
    if numel(r)~=numel(rActual) || max(abs(r-rActual))>1e-12*a0
        error('step44:RadialMismatch','Exact COD changed native sampling locations.');
    end
    profile{icase}=p;
    radial{icase}=r;
    target=testK(icase,:);
    maxRawDeviation(icase,:)=max(abs(bsxfun(@minus,p,target)),[],1);
    fprintf('  %-17s imposed=[%+.8g,%+.8g] raw max abs=[%.3e,%.3e]\n', ...
        caseNames{icase},target,maxRawDeviation(icase,:));
    for iw=1:size(opt.FitWindows,1)
        w=opt.FitWindows(iw,:);
        ix=(r/a0>=w(1))&(r/a0<=w(2));
        for degree=sort(unique(opt.FitDegrees(:).'))
            if nnz(ix)<4*(degree+1)
                error('step44:TooFewNativeNodes', ...
                    'Exact-field fit [%g,%g], degree %d has %d nodes.', ...
                    w,degree,nnz(ix));
            end
            coefI=polyfit(r(ix)/a0,p(ix,1),degree);
            coefII=polyfit(r(ix)/a0,p(ix,2),degree);
            ki=coefI(end);kii=coefII(end);
            normScale=max(abs(target));
            nr=nr+1;
            fits(nr,:)=[icase,w,degree,nnz(ix), ...
                ki,kii,(ki-target(1))/normScale, ...
                (kii-target(2))/normScale, ...
                sqrt(mean((polyval(coefII,r(ix)/a0)-p(ix,2)).^2)), ...
                (kii/ki)]; % pure-II ratio may be Inf; diagnostic only
        end
    end
end
fits=fits(1:nr,:);
fitTable=array2table(fits,'VariableNames', { ...
    'case_index','lower_r_over_a0','upper_r_over_a0','degree','n_native', ...
    'KI_COD','KII_COD','deltaKI_over_Kscale', ...
    'deltaKII_over_Kscale','fit_RMSE_KII','ratio_if_defined'});
% Median unit-field response estimates the signed recovery matrix.
Mcod=[median(profile{1}(:,1)) median(profile{2}(:,1)); ...
      median(profile{1}(:,2)) median(profile{2}(:,2))];
mixedPred=(Mcod*testK(3,:).').';
mixedRaw=profile{3};
linearity=max(abs(mixedRaw - ...
    (testK(3,1)*profile{1}+testK(3,2)*profile{2})),[],1);
pureLeak=max([max(abs(profile{1}(:,2))), ...
              max(abs(profile{2}(:,1)))]);
mixedKIIrel=max(abs(mixedRaw(:,2)-testK(3,2)))/abs(testK(3,2));
maxFitMixedII=max(abs(fitTable.deltaKII_over_Kscale( ...
    fitTable.case_index==3))) * max(abs(testK(3,:))) / abs(testK(3,2));
passCOD=pureLeak<=opt.MaxPureLeakageOverKI && ...
    mixedKIIrel<=opt.MaxMixedKIIErrorRel && ...
    maxFitMixedII<=opt.MaxMixedKIIErrorRel && ...
    norm(Mcod-eye(2),'fro')<=1e-8;
fprintf('\nCOD EXACT-FIELD GATES\n');
fprintf('  ||M_COD-I||F=%.3e, max cross-mode leakage=%.3e\n', ...
    norm(Mcod-eye(2),'fro'),pureLeak);
fprintf('  mixed relative KII raw=%.3e, fits=%.3e, superposition=%.3e\n', ...
    mixedKIIrel,maxFitMixedII,max(linearity));
fprintf('  COD acceptance: %s\n',mat2str(passCOD));

O44=struct('settings',opt,'checkpointPath',cp, ...
    'regression',struct('passed',passRegression,'maxDr',dr,'maxDK',dk), ...
    'caseNames',{caseNames},'testK',testK,'rNative',{radial}, ...
    'apparent',{profile},'maxRawDeviation',maxRawDeviation, ...
    'fitTable',fitTable,'CODmatrix',Mcod,'mixedPred',mixedPred, ...
    'mixedLinearityResidual',linearity, ...
    'CODgates',struct('passed',passCOD,'maxPureLeakage',pureLeak, ...
      'mixedKIIrel',mixedKIIrel,'maxFitMixedKIIrel',maxFitMixedII), ...
    'EDI',struct('performed',false));
codFile=[prefix '_cod_small_data.mat'];
save(codFile,'O44');
fprintf('  COD small result saved: %s\n',codFile);

if opt.RunEDI
    % TWO large exact-field integrals only. Cache after EACH one.
    % The tiny-mixed EDI result is reconstructed by strict linearity:
    % M_EDI*[KI;KII]; it is NOT an independently evaluated third case.
    ri=O38.rInner;
    ro=rat(iref)*a0;
    if ~(ri>=2*O38.refinedTip && ro>ri)
        error('step44:BadAnnulus','Existing EDI annulus not admissible.');
    end
    progressFile=[prefix '_edi_progress.mat'];
    progress=struct('key',checkpointKey,'ri',ri,'ro',ro, ...
        'K',nan(2,2),'done',false(1,2));
    if exist(progressFile,'file')==2
        cpSaved=load(progressFile,'progress');
        prior=cpSaved.progress;
        if ~strcmp(prior.key,checkpointKey) || prior.ri~=ri || ...
                prior.ro~=ro || ~isequal(size(prior.K),[2 2])
            error('step44:StaleEDIProgress', ...
                'EDI progress belongs to a different checkpoint/domain.');
        end
        progress=prior;
    end
    if ~isfield(mat,'Dmat') && isfield(mat,'D')
        mat.Dmat=mat.D;
    end
    fprintf('\nEXACT-FIELD EDI (TWO UNIT INTEGRATIONS, CACHED)\n');
    for j=1:2
        if progress.done(j)
            fprintf('  Unit case %d cached.\n',j);
            continue
        end
        % Full-mesh exact field needed only for the EDI extraction.
        % This is EVALUATION of analytical fields, not an FEM solve.
        xEval=xl;
        xEval(faceIDs,2)=0; % same roundoff-only branch-cut evaluation
        uu=exact_williams_displacement_audit(xEval,j==1,j==2, ...
            mat.E,mat.nu,mat.ps, ...
            'UpperFaceIDs',upper,'LowerFaceIDs',lower);
        clear xEval
        % Form global displacements before entering the existing EDI.
        uv=reshape(uu,2,[]).'*R.';
        clear uu
        synthetic=reshape(uv.',[],1);
        clear uv
        [ki,kii]=SIF_LEFM_interaction_EDI( ...
            mesh,synthetic,crack.Pmid,mat, ...
            struct('r_inner',ri,'r_outer',ro), ...
            'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
            'WeightFunction','fe_nodal','QuadratureRule',16, ...
            'StoreGPDiagnostics',false);
        clear synthetic
        progress.K(:,j)=[ki;kii];
        progress.done(j)=true;
        tmp=[progressFile '.incomplete.mat'];
        save(tmp,'progress');
        [ok,msg]=movefile(tmp,progressFile,'f');
        if ~ok,error('step44:ProgressSave','%s',msg);end
        fprintf('  unit %d complete: recovered=[%+.10e,%+.10e]\n',j,ki,kii);
    end
    Medi=progress.K;
    inferredMix=Medi*testK(3,:).';
    O44.EDI=struct('performed',true,'rInner',ri,'rOuter',ro, ...
        'recoveryMatrix',Medi, ...
        'matrixMisfit',norm(Medi-eye(2),'fro'), ...
        'pureILeakageOverKI',abs(Medi(2,1)), ...
        'inferredMixed',inferredMix, ...
        'inferredMixedKIIrel', ...
            (inferredMix(2)-testK(3,2))/abs(testK(3,2)), ...
        'isLinearityReconstruction',true);
    fprintf('  ||M_EDI-I||F = %.6g; pure-I KII leakage = %.6g\n', ...
        O44.EDI.matrixMisfit,O44.EDI.pureILeakageOverKI);
    fprintf('  INFERRED mixed KII relative difference = %+.6g\n', ...
        O44.EDI.inferredMixedKIIrel);
end
save(codFile,'O44'); % compact, without mesh, synthetic U or stiffness
fprintf('STEP 44 finished. No new FEM solve, mesh or stiffness matrix.\n');
end

function need(S,field)
if ~isstruct(S)||~isfield(S,field)||isempty(S.(field))
    error('step44:MissingInput','Required field %s is absent.',field);
end
end
