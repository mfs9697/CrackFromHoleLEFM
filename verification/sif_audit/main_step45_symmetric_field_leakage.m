function O45=main_step45_symmetric_field_leakage(P45,varargin)
%MAIN_STEP45_SYMMETRIC_FIELD_LEAKAGE
% Measure spurious signed KII in an ACTUAL symmetric zero-theta FEM field.
% The checkpoint is created by main_step45_prepare_symmetric_checkpoint.
% Default COD-only phase: NO FEM solve, mesh generation or EDI integration.
% Explicit RunEDI=true evaluates the SAME saved FEM field with 16-point
% FE-nodal q interaction EDI; progress is cached after EACH outer radius.
% This is not Step44's exact prescribed Williams field.
%
% Examples:
%   O45=main_step45_symmetric_field_leakage(P45);
%   O45=main_step45_symmetric_field_leakage(P45,'RunEDI',true);
%   O45=main_step45_symmetric_field_leakage(P45, ...
%       'RunEDI',true,'ROuterOverA0',[.50 .65 .80]);
%
% The 1e-6 ratio is a proposed VERIFICATION TARGET, not a confidence band.
if isstruct(P45)
    if ~isfield(P45,'checkpointPath')
        error('step45:BadInput','P45 must contain checkpointPath.');
    end
    cp=char(P45.checkpointPath);
else
    cp=char(P45);
end
ip=inputParser;
addParameter(ip,'RunEDI',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'ROuterOverA0',0.65, ...
    @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&& ...
        all(isfinite(x))&&all(x>0)&&all(x<1));
addParameter(ip,'InnerRadius',[], ...
    @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&& ...
        isfinite(x)&&x>0));
addParameter(ip,'FitWindows',[.04 .30;.08 .30;.12 .30], ...
    @(x)isnumeric(x)&&size(x,2)==2&&all(isfinite(x(:)))&& ...
        all(x(:,1)>0)&all(x(:,1)<x(:,2))&&all(x(:,2)<.7));
addParameter(ip,'FitDegrees',[1 2], ...
    @(x)isnumeric(x)&&isvector(x)&& ...
        all(ismember(x,[1 2])));
addParameter(ip,'MinFitNodes',8, ...
    @(x)isnumeric(x)&&isscalar(x)&&x==round(x)&&x>=6);
addParameter(ip,'TargetLeakageRatio',1e-6, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'TinyRatioReference',1.0795940665e-4, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'SavePrefix',fullfile(pwd, ...
    'step45_symmetric_field_leakage'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});
opt=ip.Results;
if exist(cp,'file')~=2
    error('step45:CheckpointAbsent', ...
        'Saved symmetric FEM checkpoint does not exist: %s',cp);
end
f=dir(cp);
key=sprintf('%s|%d|%.16g',cp,f.bytes,f.datenum);
s=load(cp,'mesh','U','mat','crack','a0','meta');
for field={'mesh','U','mat','crack','a0','meta'}
    need(s,field{1});
end
mesh=s.mesh;U=s.U;mat=s.mat;crack=s.crack;
a0=s.a0;meta=s.meta;
clear s
if ~isfield(meta,'caseType') || ...
        ~strcmp(meta.caseType,'step45_centered_half_theta0') || ...
        ~isfield(meta,'zeroAngleDeg') || meta.zeroAngleDeg~=0
    error('step45:WrongCheckpoint', ...
        'Not the approved centered-hole theta=0 FEM checkpoint.');
end
if numel(U)~=2*size(mesh.coord,1) || ...
        any(~isfinite(U)) || ...
        abs(norm(diff(crack.Pmid,1,1))/a0-1)>1e-10
    error('step45:DataMismatch', ...
        'Saved displacement or crack geometry mismatch.');
end
prefix=char(opt.SavePrefix);
[folder,~,~]=fileparts(prefix);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
fprintf('\n============================================================\n');
fprintf('STEP 45: GENUINE SYMMETRIC FEM MODE-II LEAKAGE\n');
fprintf('============================================================\n');
fprintf('  Source: %s; Npoly=%d, T6 nodes=%d\n', ...
    meta.source,meta.Npoly,size(mesh.coord,1));
fprintf('  Physical symmetry requires KII=0. No new FEM solve.\n');

[r,apparent,face]=native_COD_audit(mesh,U,mat,crack,8,true);
faceSide=face.faceSide;
if face.nUpper<8 || face.nLower<8
    error('step45:FaceCount','Insufficient native face samples.');
end
if ~isfinite(face.gridMismatch) || ...
        face.gridMismatch>max(1e-12,1e-8*a0)
    warning('step45:FaceMismatch', ...
        ['Upper/lower radial grids differ: opposite-face interpolation ', ...
         'may bias very small tangential jumps.']);
end
x=r/a0;
q=apparent(:,2)./apparent(:,1);
if any(~isfinite(q)) || ~any(apparent(:,1)>0)
    error('step45:InvalidCOD','COD normal component or ratios invalid.');
end
% Bands are descriptive. Very first tip nodes may be dominated by
% conventional T6 singular-field approximation error.
bands=[0 .04;.04 .08;.08 .12;.12 .20;.20 .30];
raw=nan(size(bands,1),6);
for i=1:size(bands,1)
    idx=(x>=bands(i,1)&x<bands(i,2));
    raw(i,1:3)=[bands(i,:),nnz(idx)];
    if nnz(idx)>0
        raw(i,4:6)=[median(q(idx)),min(q(idx)),max(q(idx))];
    end
end
rawTable=array2table(raw,'VariableNames',{ ...
    'lower_r_over_a0','upper_r_over_a0','n_native', ...
    'median_raw_ratio','min_raw_ratio','max_raw_ratio'});

win=opt.FitWindows;
degrees=sort(unique(opt.FitDegrees(:).'));
rows=nan(size(win,1)*numel(degrees),9);
n=0;
for iw=1:size(win,1)
    ids=(x>=win(iw,1)&x<=win(iw,2));
    for deg=degrees
        n=n+1;
        rows(n,1:4)=[win(iw,:),deg,nnz(ids)];
        if nnz(ids)<max(opt.MinFitNodes,4*(deg+1))
            fprintf('  COD [%g,%g] degree=%d SKIP (only %d native nodes)\n', ...
                win(iw,:),deg,nnz(ids));
            continue
        end
        a=polyfit(x(ids),apparent(ids,1),deg);
        b=polyfit(x(ids),apparent(ids,2),deg);
        ki=a(end);kii=b(end);
        if ki<=0 || ~isfinite(ki)
            warning('step45:BadKI','Nonpositive fitted KI in one COD window.');
            continue
        end
        rows(n,5:9)=[ki,kii,kii/ki, ...
            abs(kii/ki)/opt.TinyRatioReference, ...
            sqrt(mean((polyval(b,x(ids))-apparent(ids,2)).^2))];
    end
end
fitTable=array2table(rows,'VariableNames',{ ...
    'lower_r_over_a0','upper_r_over_a0','degree','n_native', ...
    'KI_COD','KII_COD','ratio_COD', ...
    'abs_leakage_as_fraction_of_tiny_signal','RMSE_KII'});
accepted=fitTable.ratio_COD(isfinite(fitTable.ratio_COD));
if isempty(accepted)
    maxCOD=NaN;
    passCOD=false;
    fprintf('  COD gate NOT EVALUATED: insufficient qualified fitting windows.\n');
else
    maxCOD=max(abs(accepted));
    passCOD=maxCOD<=opt.TargetLeakageRatio;
end
fprintf('\nRAW SYMMETRIC NATIVE COD BY BANDS\n');disp(rawTable);
fprintf('\nSYMMETRIC FEM COD INTERCEPTS\n');disp(fitTable);
fprintf('  max finite fitted |KII/KI| = %.8e\n',maxCOD);
fprintf('  leakage target %.2e met (COD, available fits): %s\n', ...
    opt.TargetLeakageRatio,mat2str(passCOD));
if isfield(meta,'previousStep18EDI') && ...
        meta.previousStep18EDI.available
    fprintf(['  Note: Step18 historical 7-point EDI reference is ', ...
        'stored in O45.meta.previousStep18EDI, not treated ', ...
        'as a matched Step45 EDI result.\n']);
end
H=local_tip_scale(mesh,crack.Pmid(end,:));
hTip=H.median;
fprintf(['  measured tip-edge median=%.8g m, htip/a0=%.6g, ', ...
    'tip-adjacent T3 elements=%d (above/below=%d/%d)\n'], ...
    hTip,hTip/a0,H.nTipElements,H.nAbove,H.nBelow);
O45=struct('settings',opt,'checkpointPath',cp,'meta',meta, ...
    'tipEdgeMedian',hTip,'tipTopology',H,'nativeR',r,'nativeApparent',apparent, ...
    'faceCounts',[face.nUpper,face.nLower], ...
    'faceGridMismatch',face.gridMismatch, ...
    'rawBands',rawTable,'fitTable',fitTable, ...
    'CODgates',struct('maxAbsRatio',maxCOD, ...
      'targetRatio',opt.TargetLeakageRatio, ...
      'passed',passCOD, ...
      'evaluated',~isempty(accepted)), ...
    'EDI',struct('performed',false));
clear mesh U mat face faceSide
smallFile=[prefix '_small_data.mat'];
save(smallFile,'O45'); % only compact arrays, not FEM field
fprintf('  Saved compact COD results: %s\n',smallFile);

if ~opt.RunEDI
    fprintf(['  COD-only phase done. To compare both extractors ', ...
        'run with ''RunEDI'',true.\n']);
    return
end
% Only when explicitly requested, reload exactly the same saved field.
s=load(cp,'mesh','U','mat','crack');
mesh=s.mesh;U=s.U;mat=s.mat;crack=s.crack;
clear s
if ~isfield(mat,'Dmat') && isfield(mat,'D')
    mat.Dmat=mat.D;
end
rRat=sort(unique(opt.ROuterOverA0(:).'));
ri=zeros(size(rRat));
for i=1:numel(rRat)
    if isempty(opt.InnerRadius)
        % Matches original Step18 radius selection. Integration order,
        % unlike Step18, is explicitly fixed to 16 points.
        ri(i)=max(0.10*rRat(i)*a0,2*hTip);
    else
        ri(i)=opt.InnerRadius;
    end
    if ri(i)<2*hTip || ri(i)>=rRat(i)*a0
        error('step45:BadEDIDomain', ...
            'At outer/a0=%.2f need 2*htip<=r_inner<r_outer.',rRat(i));
    end
end
progressFile=[prefix '_edi_progress.mat'];
progress=struct('key',key,'rat',rRat,'inner',ri,'rule',16, ...
    'method','fe_nodal','KI',nan(size(rRat)), ...
    'KII',nan(size(rRat)),'done',false(size(rRat)));
if exist(progressFile,'file')==2
    pre=load(progressFile,'progress');
    prev=pre.progress;
    if ~strcmp(prev.key,key) || ~isequal(prev.rat,rRat) || ...
            ~isequal(prev.inner,ri) || prev.rule~=16 || ...
            ~strcmp(prev.method,'fe_nodal')
        error('step45:StaleEDIProgress', ...
            'Saved EDI progress uses different field or quadrature settings.');
    end
    if any(~isfinite(prev.KI(prev.done))) || ...
            any(~isfinite(prev.KII(prev.done)))
        error('step45:CorruptedEDIProgress', ...
            'Cached completed EDI entries are not finite.');
    end
    progress=prev;
end
fprintf('\n16-POINT SYMMETRIC FEM EDI (SAVED FIELD)\n');
for i=1:numel(rRat)
    if progress.done(i)
        fprintf('  EDI r_outer/a0 %.2f reused from cache.\n',rRat(i));
        continue
    end
    [ki,kii]=SIF_LEFM_interaction_EDI( ...
        mesh,U,crack.Pmid,mat, ...
        struct('r_inner',ri(i),'r_outer',rRat(i)*a0), ...
        'UsePlaneStrain',mat.ps==1,'Verbose',false, ...
        'WeightFunction','fe_nodal','QuadratureRule',16, ...
        'StoreGPDiagnostics',false);
    if ~isfinite(ki) || ~isfinite(kii) || ki<=0
        error('step45:InvalidEDI', ...
            'The actual symmetric FEM EDI returned invalid KI or KII.');
    end
    progress.KI(i)=ki;
    progress.KII(i)=kii;
    progress.done(i)=true;
    tmp=[progressFile '.incomplete.mat'];
    save(tmp,'progress');
    [ok,msg]=movefile(tmp,progressFile,'f');
    if ~ok,error('step45:EDISave','%s',msg);end
    fprintf(['  EDI r_outer/a0=%.2f: KI=%.9e, KII=%+.9e, ', ...
        'KII/KI=%+.9e (cached)\n'], ...
        rRat(i),ki,kii,kii/ki);
end
qEDI=progress.KII./progress.KI;
resultsTable=table(rRat(:),ri(:),progress.KI(:), ...
    progress.KII(:),qEDI(:), ...
    abs(qEDI(:))/opt.TinyRatioReference, ...
    'VariableNames',{'r_outer_over_a0','r_inner', ...
    'KI','KII','ratio','abs_leakage_as_fraction_of_tiny_signal'});
maxEDI=max(abs(qEDI));
O45.EDI=struct('performed',true,'table',resultsTable, ...
    'maxAbsRatio',maxEDI, ...
    'targetRatio',opt.TargetLeakageRatio, ...
    'passed',maxEDI<=opt.TargetLeakageRatio, ...
    'domainRatioSpread',max(qEDI)-min(qEDI), ...
    'progressFile',progressFile);
fprintf('\n16-POINT SYMMETRIC FEM EDI RESULTS\n');disp(resultsTable);
fprintf('  EDI target met across requested domains: %s\n', ...
    mat2str(O45.EDI.passed));
fprintf('  COD target met across evaluable windows: %s\n', ...
    mat2str(passCOD));
fprintf(['  These are actual-FEM symmetry residuals, not ', ...
    'uncertainty bounds on the asymmetric problem.\n']);
clear mesh U mat
save(smallFile,'O45');
fprintf('STEP 45 completed on saved FEM field. No new FEM solve.\n');
end

function H=local_tip_scale(mesh,tip)
P=mesh.coord3;
T=mesh.connect3;
tol=max(1e-12,1e-8*max(1,max(abs(P(:)))));
radius=hypot(P(:,1)-tip(1),P(:,2)-tip(2));
ids=find(radius<=min(radius)+tol);
near=T(any(ismember(T,ids),2),:);
if isempty(near)
    error('step45:TipGeometry','No tip-adjacent T3 elements.');
end
p1=P(near(:,1),:);p2=P(near(:,2),:);p3=P(near(:,3),:);
L=[hypot(p1(:,1)-p2(:,1),p1(:,2)-p2(:,2)); ...
   hypot(p2(:,1)-p3(:,1),p2(:,2)-p3(:,2)); ...
   hypot(p3(:,1)-p1(:,1),p3(:,2)-p1(:,2))];
L=L(isfinite(L)&L>tol);
if isempty(L),error('step45:TipScale','No valid tip edges.');end
cen=(p1+p2+p3)/3;
H=struct('median',median(L), ...
    'nTipElements',size(near,1), ...
    'nAbove',nnz(cen(:,2)>tip(2)+tol), ...
    'nBelow',nnz(cen(:,2)<tip(2)-tol), ...
    'nOnAxis',nnz(abs(cen(:,2)-tip(2))<=tol));
end

function need(s,f)
if ~isstruct(s)||~isfield(s,f)||isempty(s.(f))
    error('step45:MissingField','Required field %s is absent.',f);
end
end
