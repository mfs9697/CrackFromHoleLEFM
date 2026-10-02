function R53=main_step53_common_grid_parity(varargin)
%MAIN_STEP53_COMMON_GRID_PARITY  Paired-mesh, IDENTICAL-POINT parity audit.
% Step52 used UPPER MESH NODES as sample points separately per mesh.
% Step53 removes that sampling difference: the SAME predeclared spatial
% (r,theta) positions above and below the crack are evaluated on BOTH
% previously saved FEM meshes with their T6 shape functions.
%
% Four checks/measurements:
%   - affine even-x/odd-y field (exact T6 polynomial self-check)
%   - direct analytical pure-I mirror symmetry (branch convention)
%   - T6-interpolated pure-I exact NODAL synthetic field (mesh artifact)
%   - already solved FEM field (numerical reflection parity)
%
% No new mesh, FEM solve, EDI, polynomial SIF fit or boundary changes.
% Compare ONLY shared physically located grid points; use the SAME
% physical native opening reference radii and shared FEM denominator.
%
%   addpath(genpath(pwd));
%   R53=main_step53_common_grid_parity();
%   disp(R53.summary);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');
ip=inputParser;
addParameter(ip,'Step48File',fullfile(vdir, ...
    'step48_refined_matched_edi_comparison_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'SaveFile',fullfile(vdir, ...
    'step53_common_grid_parity_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;
addpath(genpath(root));
x=load(char(opt.Step48File),'R48');
if ~isfield(x,'R48') || ~isfield(x.R48,'comparison') || ...
        ~isfield(x.R48,'baselineSource') || ...
        ~isfield(x.R48,'refined')
    error('step53:MissingMatchedResult', ...
        'The previously completed Step48 saved R48 is required.');
end
R48=x.R48;
T=R48.comparison;
if ~istable(T) || height(T)~=2 || ...
        ~all(ismember({'r_inner','r_outer'}, ...
                       T.Properties.VariableNames))
    error('step53:MissingMatchedAnnulus', ...
        'Step48 must contain its matched original/refined EDI annulus.');
end
ri=T.r_inner(1);ro=T.r_outer(1);
if ~all(isfinite([ri,ro])) || ri<=0 || ro<=ri || ...
        max(abs(T.r_inner-ri))>1e-12 || ...
        max(abs(T.r_outer-ro))>1e-12
    error('step53:AnnuliDiffer','Original/refined absolute annuli differ.');
end
o=load(char(R48.baselineSource),'O45');
if ~isfield(o,'O45') || ~isfield(o.O45,'checkpointPath') || ...
        ~isfield(R48.refined,'checkpointPath')
    error('step53:MissingCheckpointReference', ...
        'Step48 must refer to both previously solved checkpoints.');
end
paths={char(o.O45.checkpointPath), ...
       char(R48.refined.checkpointPath)};
names={'Original Step45','Refined Step47'};
regions={'Common COD disk','Common EDI annulus'};
rRatios={[.08 .12 .16 .20 .24 .28], ...
         [.24 .30 .36 .42 .48 .54 .60]};
anglesDeg=[25 45 65 90 115 135 155]; % avoid slit, all positive y
refRatios=[.14 .18 .22 .26]; % SAME fixed radii on both native faces
S=cell(2,1); referenceOpening=nan(2,1);
syntheticOpening=nan(2,1);
tip0=[];
a0=0;
fprintf('\n============================================================\n');
fprintf('STEP 53: IDENTICAL PHYSICAL POINTS ON TWO SOLVED FEM MESHES\n');
fprintf('============================================================\n');
fprintf('  NO FEM solve, mesh generation, EDI or SIF fit.\n');
fprintf('  Step48 matched annulus: [%.12g, %.12g] m.\n',ri,ro);
for imesh=1:2
    cp=paths{imesh};
    if exist(cp,'file')~=2
        error('step53:CheckpointMissing','Saved FEM file absent: %s',cp);
    end
    d=load(cp,'mesh','U','mat','crack','a0','meta');
    mustfields={'mesh','U','mat','crack','a0','meta'};
    for j=1:numel(mustfields)
        if ~isfield(d,mustfields{j})
            error('step53:BadCheckpoint', ...
                'Missing saved field %s.',mustfields{j});
        end
    end
    if ~isfield(d.meta,'caseType') || ...
            ~strcmp(d.meta.caseType,'step45_centered_half_theta0') || ...
            d.meta.Npoly~=240 || abs(d.a0-.004)>1e-12 || ...
            size(d.mesh.connect,2)~=6 || ...
            ~isequal(d.mesh.connect(:,1:3),d.mesh.connect3) || ...
            numel(d.U)~=2*size(d.mesh.coord,1) || ...
            any(~isfinite(d.U))
        error('step53:NotApprovedControl','Unexpected FEM checkpoint.');
    end
    tip=d.crack.Pmid(end,:);
    if norm(diff(d.crack.Pmid,1,1)-[d.a0 0])>1e-12
        error('step53:NotHorizontal','Control crack must be horizontal.');
    end
    if imesh==1
        tip0=tip;
        a0=d.a0;
        mat0=d.mat;
        Pmid0=d.crack.Pmid;
    else
        if ~isfield(d.meta,'stage') || ...
                ~strcmp(d.meta.stage,'step47_refined_control') || ...
                norm(d.crack.Pmid-Pmid0,'fro')>1e-12 || ...
                abs(d.mat.E-mat0.E)>1e-8 || ...
                abs(d.mat.nu-mat0.nu)>1e-12 || ...
                d.mat.ps~=mat0.ps
            error('step53:PhysicalProblemChanged', ...
                'The two checkpoints are not the same physical control.');
        end
    end
    [faceR,app,side]=native_COD_audit( ...
        d.mesh,d.U,d.mat,d.crack,8,true);
    X=d.mesh.coord;
    Xlocal=X-tip;
    faceIDs=find(side.faceSide~=0);
    Xlocal(faceIDs,2)=0; % analytical evaluation only, never move FEM nodes
    Usyn=exact_williams_displacement_audit( ...
        Xlocal,1,0,d.mat.E,d.mat.nu,d.mat.ps, ...
        'UpperFaceIDs',find(side.faceSide==1), ...
        'LowerFaceIDs',find(side.faceSide==-1));
    [rCheck,appCheck]=native_COD_audit( ...
        d.mesh,Usyn,d.mat,d.crack,8);
    if numel(rCheck)~=numel(faceR) || ...
            max(abs(rCheck-faceR))>1e-12 || ...
            max(abs(appCheck(:,1)-1))>1e-8 || ...
            max(abs(appCheck(:,2)))>1e-8
        error('step53:ExactCODSelfCheck', ...
            'Synthetic pure-I nodal field failed COD exact self-check.');
    end
    mu=d.mat.E/(2*(1+d.mat.nu));
    if d.mat.ps==1
        kappa=3-4*d.mat.nu;
    else
        kappa=(3-d.mat.nu)/(1+d.mat.nu);
    end
    scale=mu/(kappa+1)*sqrt(2*pi./faceR);
    normalJump=app(:,1)./scale;
    sameReferenceR=refRatios*a0;
    if min(faceR)>min(sameReferenceR) || ...
            max(faceR)<max(sameReferenceR)
        error('step53:FaceReferenceOutOfRange', ...
            'Same physical crack-opening reference not spanned by native nodes.');
    end
    openingAtRef=interp1(faceR,normalJump,sameReferenceR,'pchip');
    if any(~isfinite(openingAtRef)) || any(openingAtRef<=0)
        error('step53:InvalidOpening','Reference mode-I opening invalid.');
    end
    referenceOpening(imesh)=mean(openingAtRef);
    syntheticOpening(imesh)=mean(1./( ...
        mu/(kappa+1)*sqrt(2*pi./sameReferenceR)));
    tri=triangulation(d.mesh.connect3,d.mesh.coord3);
    uActual=reshape(d.U,2,[]).';
    uSyn=reshape(Usyn,2,[]).';
    uAffine=[X(:,1)-tip(1),X(:,2)-tip(2)];
    out=cell(1,2);
    for ir=1:2
        [rr,theta]=ndgrid(rRatios{ir}*a0,anglesDeg);
        rr=rr(:);
        theta=theta(:);
        xQ=tip(1)+rr.*cosd(theta);
        y=rr.*sind(theta);
        if any(y<=0)
            error('step53:WrongGrid','Every query must be strictly off-face.');
        end
        qPlus=[xQ,tip(2)+y];
        qMinus=[xQ,tip(2)-y];
        [tPlus,bPlus]=pointLocation(tri,qPlus);
        [tMinus,bMinus]=pointLocation(tri,qMinus);
        valid=isfinite(tPlus)&isfinite(tMinus);
        aPlus=local_eval(d.mesh,uActual,tPlus,bPlus);
        aMinus=local_eval(d.mesh,uActual,tMinus,bMinus);
        sPlus=local_eval(d.mesh,uSyn,tPlus,bPlus);
        sMinus=local_eval(d.mesh,uSyn,tMinus,bMinus);
        fPlus=local_eval(d.mesh,uAffine,tPlus,bPlus);
        fMinus=local_eval(d.mesh,uAffine,tMinus,bMinus);
        out{ir}=struct('valid',valid,'actualP',aPlus,'actualM',aMinus, ...
            'synP',sPlus,'synM',sMinus, ...
            'affP',fPlus,'affM',fMinus, ...
            'qPlus',qPlus,'qMinus',qMinus);
    end
    S{imesh}=out;
    fprintf('  %s: native face reference opening %.9e m\n', ...
        names{imesh},referenceOpening(imesh));
end
if ~all(isfinite([referenceOpening;syntheticOpening])) || ...
        min([referenceOpening;syntheticOpening])<=0 || ...
        max(abs(syntheticOpening-syntheticOpening(1)))>1e-12
    error('step53:ReferenceMismatch', ...
        'Same geometry and material must use compatible references.');
end
% ONE SHARED PHYSICAL normalization for both measured FEM fields:
normAct=mean(referenceOpening);
% Another ONE SHARED exact-analytical normalization for prescribed KI=1:
normSyn=mean(syntheticOpening);
rows=cell(4,1);
gridInfo=cell(2,1);
k=0;
for ir=1:2
    A=S{1}{ir};
    B=S{2}{ir};
    if any(abs(A.qPlus(:)-B.qPlus(:))>1e-12) || ...
            any(abs(A.qMinus(:)-B.qMinus(:))>1e-12)
        error('step53:NotCommonPoints','Physical sample grids differ.');
    end
    joint=A.valid&B.valid;
    nJoint=nnz(joint);
    nTotal=numel(joint);
    fprintf('  %s: common located points %d/%d', ...
        regions{ir},nJoint,nTotal);
    fprintf(' (original %d, refined %d).\n', ...
        nnz(A.valid),nnz(B.valid));
    if nJoint<.8*nTotal
        error('step53:LowCommonCoverage', ...
            'Cannot compare meshes with less than 80%% common query coverage.');
    end
    qP=A.qPlus(joint,:);
    qM=A.qMinus(joint,:);
    directP=exact_williams_displacement_audit( ...
        qP-tip0,1,0,mat0.E,mat0.nu,mat0.ps);
    directM=exact_williams_displacement_audit( ...
        qM-tip0,1,0,mat0.E,mat0.nu,mat0.ps);
    directP=reshape(directP,2,[]).';
    directM=reshape(directM,2,[]).';
    directError=max(abs([directP(:,1)-directM(:,1); ...
                         directP(:,2)+directM(:,2)]));
    if directError>1e-10*normSyn
        error('step53:AnalyticParityFailed', ...
            'Direct pure-I Williams mirror symmetry failed.');
    end
    gridInfo{ir}=struct('region',regions{ir}, ...
        'r_over_a0',rRatios{ir},'anglesDeg',anglesDeg, ...
        'nPredeclared',nTotal, ...
        'nCommonLocated',nJoint, ...
        'commonUpperPoints',qP, ...
        'commonLowerPoints',qM);
    for imesh=1:2
        if imesh==1,Z=A;else,Z=B;end
        aP=Z.actualP(joint,:);aM=Z.actualM(joint,:);
        sP=Z.synP(joint,:);sM=Z.synM(joint,:);
        fP=Z.affP(joint,:);fM=Z.affM(joint,:);
        affineError=max(abs([fP(:,1)-fM(:,1); ...
                             fP(:,2)+fM(:,2)]));
        if affineError>1e-10*a0
            error('step53:T6AffineFailed', ...
                'T6 barycentric interpolation failed exact affine parity.');
        end
        ax=aP(:,1)-aM(:,1);
        ay=aP(:,2)+aM(:,2);
        sy=sP(:,2)+sM(:,2);
        sx=sP(:,1)-sM(:,1);
        ay=ay-median(ay);  % remove shared physical y-translation gauge
        sy=sy-median(sy);
        k=k+1;
        rows{k}=struct( ...
            'mesh',names{imesh},'region',regions{ir}, ...
            'nFixedGrid',nTotal,'nCommonCompared',nJoint, ...
            'maxAffineError_m',affineError, ...
            'maxAnalyticalPureIParityError',directError, ...
            'rmsActualEvenUx_m',sqrt(mean(ax.^2)), ...
            'rmsActualGaugeOddUy_m',sqrt(mean(ay.^2)), ...
            'rmsActualEvenUx_overSharedOpenY', ...
                         sqrt(mean(ax.^2))/normAct, ...
            'rmsActualGaugeOddUy_overSharedOpenY', ...
                         sqrt(mean(ay.^2))/normAct, ...
            'rmsPureINodalEvenUx_overExactOpenY', ...
                         sqrt(mean(sx.^2))/normSyn, ...
            'rmsPureINodalGaugeOddUy_overExactOpenY', ...
                         sqrt(mean(sy.^2))/normSyn);
    end
end
summary=struct2table(vertcat(rows{:}));
fprintf('\nSTEP 53: SAME-LOCATION TWO-MESH REFLECTION PARITY\n');
disp(summary);
fprintf('  Shared actual FEM normal opening reference = %.9e m.\n',normAct);
fprintf('  Shared synthetic KI=1 normal opening reference = %.9e m.\n',normSyn);
fprintf(['  This is a fixed-grid displacement-parity comparison, ', ...
    'NOT a KII uncertainty bound or correction to EDI.\n']);
R53=struct('summary',summary, ...
    'rInner',ri,'rOuter',ro, ...
    'referenceRadiiOverA0',refRatios, ...
    'referenceActualByMesh',referenceOpening, ...
    'sharedActualOpenY',normAct,'sharedSyntheticOpenY',normSyn, ...
    'sampleGrids',{gridInfo}, ...
    'noNewFEM',true,'noEDI',true);
path=char(opt.SaveFile);
[folder,~,~]=fileparts(path);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
save(path,'R53'); % summary + small fixed point grid only
fprintf('  Compact fixed-grid report: %s\n',path);
end

function interp=local_eval(mesh,field,tri,bary)
n=numel(tri);
interp=nan(n,2);
valid=isfinite(tri);
if ~any(valid),return;end
idx=find(valid);
a=bary(valid,1);b=bary(valid,2);c=bary(valid,3);
N=[a.*(2*a-1),b.*(2*b-1),c.*(2*c-1), ...
    4*a.*b,4*b.*c,4*c.*a];
ids=mesh.connect(tri(valid),:);
for j=1:numel(idx)
    interp(idx(j),:)=N(j,:)*field(ids(j,:),:);
end
end
