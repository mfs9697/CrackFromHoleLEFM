function [P63,O63]=main_step63_calibrated_asymmetric_physical_solve(varargin)
%MAIN_STEP63_CALIBRATED_ASYMMETRIC_PHYSICAL_SOLVE
% Guarded ONE-SOLVE driver for the exact Step62B calibrated asymmetric mesh.
%
% DEFAULT IS SAFE: AllowSolve=false.
%
% The function:
%   1) validates the exact locally selected C03 Step62B candidate;
%   2) NEVER regenerates/remeshes the candidate;
%   3) performs exactly ONE explicitly authorized physical FEM solve;
%   4) saves a unique solved-field checkpoint before postprocessing;
%   5) runs native COD only on the saved field.
%
% NO physical EDI is performed here. EDI is intentionally reserved for a
% separate later step after the physical displacement/COD result is reviewed.
%
% Safe/default:
%   [P63,O63]=main_step63_calibrated_asymmetric_physical_solve();
%
% Authorized:
%   [P63,O63]=main_step63_calibrated_asymmetric_physical_solve( ...
%       'AllowSolve',true);
%
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');

ip=inputParser;
addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'CandidateFile', ...
    fullfile(vdir,'step62b_calibrated_mesh_selected_candidate_T3.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'CalibrationFile', ...
    fullfile(vdir,'step62b_mesh_calibration_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'CheckpointFile', ...
    fullfile(vdir,'step63_calibrated_asymmetric_physical_solved.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'CODSavePrefix', ...
    fullfile(vdir,'step63_calibrated_asymmetric_cod'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;

addpath(genpath(root));
assert_step63_branch(root);

candFile=char(opt.CandidateFile);
calFile=char(opt.CalibrationFile);
cp=char(opt.CheckpointFile);
codPrefix=char(opt.CODSavePrefix);

if exist(candFile,'file')~=2
    error('step63:MissingCandidate', ...
        'Missing exact Step62B selected candidate: %s',candFile);
end
if exist(calFile,'file')~=2
    error('step63:MissingCalibrationReport', ...
        'Missing completed Step62B compact report: %s',calFile);
end

a=load(candFile,'candidate');
b=load(calFile,'R');
if ~isfield(a,'candidate') || ~isfield(b,'R')
    error('step63:BadInputFile', ...
        'Expected candidate and R variables in Step62B saved files.');
end
candidate=a.candidate;
R62b=b.R;

for f={'p','t','crack','mat','gates','provenance', ...
        'structuredDesign','exteriorDesign','synthetic', ...
        'scientificallyReadyForFEMProposal'}
    need(candidate,f{1});
end
for f={'selectedCandidate','selectedCalibration','selectedStep62', ...
        'syntheticTable','noPhysicalU','noPhysicalFEM'}
    need(R62b,f{1});
end

% -------------------------------------------------------------------------
% Exact Step62B scientific-selection provenance.
% -------------------------------------------------------------------------
if ~logical(candidate.scientificallyReadyForFEMProposal) || ...
        ~logical(candidate.gates.structuralPass) || ...
        ~logical(candidate.synthetic.passed) || ...
        ~R62b.noPhysicalU || ~R62b.noPhysicalFEM || ...
        ~strcmp(char(R62b.selectedCandidate),'C03') || ...
        ~logical(R62b.selectedStep62.readyForOneAsymmetricFEMProposal)
    error('step63:CandidateNotAuthorized', ...
        'The saved file is not the qualified Step62B C03 mesh.');
end

cal=R62b.selectedCalibration;
requiredCal={'transitionLength_m','farSlope','boundaryMetricGrowth', ...
    'neighborRatioTarget'};
for k=1:numel(requiredCal),need(cal,requiredCal{k});end
if abs(cal.transitionLength_m-.008)>1e-14 || ...
        abs(cal.farSlope-.10)>1e-14 || ...
        abs(cal.boundaryMetricGrowth-.25)>1e-14 || ...
        abs(cal.neighborRatioTarget-1.8)>1e-14
    error('step63:WrongCalibration', ...
        'Expected selected C03 calibration L=8 mm, farSlope=0.10.');
end

P=candidate.p;
T=candidate.t;
crack=candidate.crack;
mat0=candidate.mat;
if size(P,2)~=2 || size(T,2)~=3 || ...
        size(T,1)~=32980 || ...
        any(~isfinite(P(:))) || any(T(:)<1) || any(T(:)>size(P,1))
    error('step63:CandidateCountOrGeometry', ...
        'Exact selected C03 T3 candidate has unexpected size/entries.');
end
[P6,T6]=T3toT6_fast(P,T);
if size(P6,1)~=66854
    error('step63:T6Count','Expected selected C03 T6 node count 66854.');
end
if any(triangle_signed(P,T)<=0)
    error('step63:T3Orientation','Selected candidate contains nonpositive T3 area.');
end

candidateHash=sha256(candFile);
calibrationHash=sha256(calFile);

% Calibration/design checks independent of expected physical SIF.
maxRatio=max_neighbor_ratio(T,longest_edges(P,T));
if maxRatio>1.8+5e-12 || ...
        abs(maxRatio-R62b.selectedStep62.maxNeighborSizeRatio)>5e-10
    error('step63:GradingMismatch', ...
        'Selected candidate no longer reproduces Step62B grading result.');
end
if candidate.structuredDesign.tipIncidentTopologicalEdges~=7 || ...
        candidate.structuredDesign.tipTriangles~=6
    error('step63:TipDesignChanged','Expected six-triangle/seven-edge tip fan.');
end
if ~isfield(candidate.exteriorDesign,'calibration') || ...
        abs(candidate.exteriorDesign.calibration.transitionLength_m-.008)>1e-14 || ...
        abs(candidate.exteriorDesign.calibration.farSlope-.10)>1e-14 || ...
        abs(candidate.exteriorDesign.calibration.neighborRatioTarget-1.8)>1e-14
    error('step63:EmbeddedCalibrationMismatch', ...
        'Candidate embedded exterior calibration differs from C03.');
end

% Crack topology: distinct coincident faces, one shared tip.
up=unique(crack.upperNodes(:),'stable');
lo=unique(crack.lowerNodes(:),'stable');
tipID=crack.tipNode;
shared=intersect(up,lo);
if numel(shared)~=1 || shared~=tipID
    error('step63:CrackTopology','Crack faces must share only the tip node.');
end
a0=norm(diff(crack.Pmid,1,1));
if abs(a0-.008)>1e-12 || ...
        norm(crack.Pmid(1,:)-[.199988872196,-.0208170339049])>2e-12 || ...
        norm(crack.Pmid(end,:)-[.207985904782,-.0210349096129])>2e-12
    error('step63:CrackGeometry','Physical crack geometry changed.');
end

% Native crack-face sampling is checked BEFORE the solve using zero U.
meshTry=struct('coord3',P,'connect3',T,'coord',P6,'connect',T6);
zeroU=zeros(2*size(P6,1),1);
[r0,~,face0]=native_COD_audit(meshTry,zeroU,mat0,crack,8);
x0=r0/a0;
windows=[.04 .20;.04 .30;.08 .30;.12 .30];
sampleN=zeros(size(windows,1),1);
for k=1:size(windows,1)
    sampleN(k)=nnz(x0>=windows(k,1)&x0<=windows(k,2));
end
if face0.nUpper~=138 || face0.nLower~=138 || ...
        face0.gridMismatch>1e-12 || any(sampleN<[12;12;12;12])
    error('step63:NativeSampling','Selected candidate native sampling changed.');
end

% -------------------------------------------------------------------------
% Physical setup. Geometry comes exclusively from the selected mesh.
% Historical Stage-II physics used unit remote-y traction and minimal
% anchoring. No nominal hole geometry or mesher is invoked here.
% -------------------------------------------------------------------------
xmin=min(P(:,1)); xmax=max(P(:,1));
ymin=min(P(:,2)); ymax=max(P(:,2));
A=xmax;
B=max(abs([ymin,ymax]));
if abs(xmin)>1e-12 || abs(A-.30)>1e-12 || ...
        abs(ymin+.10)>1e-12 || abs(ymax-.10)>1e-12
    error('step63:PlateBoundary','Selected candidate plate boundary changed.');
end
if ~isfield(mat0,'E')||~isfield(mat0,'nu')||~isfield(mat0,'ps')
    error('step63:Material','Candidate material is incomplete.');
end

C=struct();
C.A=A;
C.B=B;
C.E=mat0.E;
C.nu=mat0.nu;
C.ps=mat0.ps;
C.a0=a0;
C.load=struct('type','remote_tension_y','sig0',1.0);
C.bc=struct('anchor_mode','minimal');
C.solver=struct('linear_solver','backslash','verbose',0);

[~,iLB]=min(sum((P-[0,-B]).^2,2));
[~,iRB]=min(sum((P-[A,-B]).^2,2));
if norm(P(iLB,:)-[0,-B])>1e-12 || ...
        norm(P(iRB,:)-[A,-B])>1e-12
    error('step63:Corners','Could not identify exact plate bottom corners.');
end
G=struct('p',P,'t',T);
G.edgeSets=struct('corners',struct( ...
    'left_bottom',iLB,'right_bottom',iRB));
G.meta=struct('A',A,'B',B);

fprintf('\n============================================================\n');
fprintf('STEP 63: CALIBRATED ASYMMETRIC PHYSICAL FIELD — ONE SOLVE\n');
fprintf('============================================================\n');
fprintf('  Exact selected Step62B candidate: %s\n',candFile);
fprintf('  C03: L=8 mm, farSlope=0.10, max neighbor ratio %.9g.\n',maxRatio);
fprintf('  T3=%d, T6=%d, crack a0=%.6g m, native faces=%d/%d.\n', ...
    size(T,1),size(P6,1),a0,face0.nUpper,face0.nLower);
fprintf('  Physical loading: unit remote tension y; minimal rigid-body anchoring.\n');
fprintf('  NO remesh. NO physical EDI in Step63.\n');

% Existing valid checkpoint is reused, never overwritten.
if exist(cp,'file')==2
    old=load(cp,'meta','a0');
    if ~isfield(old,'meta') || ...
            ~strcmp(old.meta.stage,'step63_calibrated_asymmetric_physical') || ...
            ~strcmp(old.meta.candidateSHA256,candidateHash) || ...
            ~strcmp(old.meta.calibrationSHA256,calibrationHash) || ...
            old.meta.nT3~=size(T,1) || old.meta.nT6~=size(P6,1) || ...
            abs(old.a0-a0)>1e-12
        error('step63:ExistingCheckpointMismatch', ...
            'Existing Step63 checkpoint does not belong to this exact candidate.');
    end
    fprintf('  Reusing existing Step63 checkpoint; NO new solve.\n');
    P63=struct('checkpointPath',cp,'newSolve',false, ...
        'reusedFile',true,'meta',old.meta);
else
    if ~opt.AllowSolve
        error('step63:ExplicitSolveApprovalRequired', ...
            ['The calibrated mesh is qualified, but the physical solve is ', ...
             'guarded. Rerun with ''AllowSolve'',true only after explicit ', ...
             'investigator authorization.']);
    end

    fprintf('  Starting exactly ONE explicitly authorized asymmetric FEM solve.\n');
    S=solve_cracked_LEFM(C,G,'lambda',1.0);

    if ~isfield(S,'mesh') || ~isfield(S,'U') || ~isfield(S,'mat') || ...
            size(S.mesh.connect3,1)~=size(T,1) || ...
            size(S.mesh.coord3,1)~=size(P,1) || ...
            size(S.mesh.coord,1)~=size(P6,1) || ...
            numel(S.U)~=2*size(P6,1) || any(~isfinite(S.U))
        error('step63:SolveOutput','Solver returned incompatible field dimensions.');
    end
    if ~isequal(S.mesh.connect3,T) || ...
            max(abs(S.mesh.coord3(:)-P(:)))>1e-12 || ...
            ~isequal(S.mesh.connect,T6) || ...
            max(abs(S.mesh.coord(:)-P6(:)))>1e-12
        error('step63:SolvedMeshChanged', ...
            'Solver did not use the exact selected C03 T3/T6 mesh.');
    end
    if abs(S.mat.E-mat0.E)>1e-12*max(1,abs(mat0.E)) || ...
            abs(S.mat.nu-mat0.nu)>1e-14 || S.mat.ps~=mat0.ps
        error('step63:SolvedMaterialChanged','Solver material differs from candidate.');
    end

    mesh=S.mesh;
    U=S.U;
    mat=S.mat;

    meta=struct( ...
        'stage','step63_calibrated_asymmetric_physical', ...
        'source','one explicitly authorized physical FEM solve on exact Step62B C03', ...
        'candidateFile',candFile, ...
        'candidateSHA256',candidateHash, ...
        'calibrationFile',calFile, ...
        'calibrationSHA256',calibrationHash, ...
        'selectedCandidate','C03', ...
        'transitionLength_m',cal.transitionLength_m, ...
        'farSlope',cal.farSlope, ...
        'neighborRatioTarget',cal.neighborRatioTarget, ...
        'maxNeighborRatio',maxRatio, ...
        'nT3',size(T,1),'nT6',size(P6,1), ...
        'a0',a0, ...
        'loadType',C.load.type,'sig0',C.load.sig0, ...
        'anchorMode',C.bc.anchor_mode, ...
        'cornerNodes',[iLB iRB], ...
        'physicalEDIperformed',false, ...
        'meshRegenerated',false);

    [folder,~,~]=fileparts(cp);
    if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
    tmp=[cp '.incomplete.mat'];
    if exist(tmp,'file')==2
        error('step63:InterruptedSave', ...
            'Inspect/remove prior incomplete checkpoint manually: %s',tmp);
    end
    save(tmp,'mesh','U','mat','crack','a0','meta','C','-v7.3');
    [ok,msg]=movefile(tmp,cp);
    if ~ok,error('step63:CheckpointSave','%s',msg);end
    fprintf('  Solved field safely checkpointed: %s\n',cp);

    P63=struct('checkpointPath',cp,'newSolve',true, ...
        'reusedFile',false,'meta',meta);

    % Release large transient solver outputs before COD.
    clear S G mesh U mat
end

% -------------------------------------------------------------------------
% COD-only physical postprocessing from SAVED checkpoint.
% -------------------------------------------------------------------------
fprintf('\nSTEP 63 PHASE 2: NATIVE COD ONLY — NO PHYSICAL EDI\n');
s=load(cp,'mesh','U','mat','crack','a0','meta');
if ~strcmp(s.meta.candidateSHA256,candidateHash) || ...
        ~strcmp(s.meta.calibrationSHA256,calibrationHash)
    error('step63:CheckpointProvenance','Checkpoint provenance mismatch.');
end

[r,app,face]=native_COD_audit(s.mesh,s.U,s.mat,s.crack,8);
rr=r/s.a0;
fitDegrees=[1 2];
fitRows=nan(size(windows,1)*numel(fitDegrees),9);
qRaw=app(:,2)./app(:,1);
n=0;
for iw=1:size(windows,1)
    ids=find(rr>=windows(iw,1)&rr<=windows(iw,2));
    for d=fitDegrees
        if numel(ids)<max(8,2*(d+1)),continue,end
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
fitTable=array2table(fitRows(1:n,:), ...
    'VariableNames',{'lower_r_over_a0','upper_r_over_a0','degree', ...
    'n_native','KI_COD','KII_COD','ratio_COD','RMSE_KII', ...
    'median_pointwise_ratio'});

bands=[0 .04;.04 .08;.08 .12;.12 .20;.20 .30];
rawRows=nan(size(bands,1),6);
for k=1:size(bands,1)
    ids=rr>=bands(k,1)&rr<bands(k,2);
    if k==size(bands,1),ids=rr>=bands(k,1)&rr<=bands(k,2);end
    rawRows(k,:)=[bands(k,:),nnz(ids), ...
        median(qRaw(ids),'omitnan'),min(qRaw(ids),[],'omitnan'), ...
        max(qRaw(ids),[],'omitnan')];
end
rawTable=array2table(rawRows,'VariableNames',{ ...
    'lower_r_over_a0','upper_r_over_a0','n_native', ...
    'median_raw_ratio','min_raw_ratio','max_raw_ratio'});

fprintf('\nRAW PHYSICAL NATIVE COD RATIOS\n');
disp(rawTable);
fprintf('\nPHYSICAL COD INTERCEPTS — NO EDI YET\n');
disp(fitTable);

O63=struct( ...
    'checkpointPath',cp, ...
    'candidateFile',candFile, ...
    'candidateSHA256',candidateHash, ...
    'calibrationFile',calFile, ...
    'calibrationSHA256',calibrationHash, ...
    'rawTable',rawTable, ...
    'fitTable',fitTable, ...
    'nativeR',r, ...
    'nativeApparent',app, ...
    'face',face, ...
    'physicalEDIperformed',false, ...
    'newFEMThisCall',P63.newSolve, ...
    'interpretation',['First physical field on calibrated C03 mesh; ', ...
      'COD only. No claim of physical KII validation before matched EDI.']);

save([codPrefix '_small_data.mat'],'O63');
fprintf('  Compact COD result saved: %s\n',[codPrefix '_small_data.mat']);
fprintf('  STOP HERE. Inspect physical COD before any EDI extraction.\n');
end

% =========================================================================
function L=longest_edges(P,T)
a=P(T(:,1),:);b=P(T(:,2),:);c=P(T(:,3),:);
L=max([vecnorm(b-c,2,2),vecnorm(a-c,2,2),vecnorm(a-b,2,2)],[],2);
end

function m=max_neighbor_ratio(T,L)
n=size(T,1);
E=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
which=repmat((1:n)',3,1);
[~,~,g]=unique(E,'rows');
cnt=accumarray(g,1);
lo=accumarray(g,L(which),[],@min);
hi=accumarray(g,L(which),[],@max);
m=max(hi(cnt==2)./lo(cnt==2));
end

function A=triangle_signed(P,T)
a=P(T(:,1),:);b=P(T(:,2),:);c=P(T(:,3),:);
A=.5*((b(:,1)-a(:,1)).*(c(:,2)-a(:,2))- ...
    (c(:,1)-a(:,1)).*(b(:,2)-a(:,2)));
end

function h=sha256(path)
f=fopen(path,'rb');
if f<0,error('step63:HashOpen','Cannot open %s.',path);end
guard=onCleanup(@()fclose(f));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(f)
    bytes=fread(f,1024*1024,'*uint8');
    md.update(typecast(bytes,'int8'));
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
clear guard
end

function assert_step63_branch(root)
[status,b]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0&&strcmp(strtrim(b),'audit/step63-calibrated-physical-solve'), ...
    'step63:Branch', ...
    'Run Step63 only on audit/step63-calibrated-physical-solve.');
end

function need(s,f)
if ~isstruct(s)||~isfield(s,f)||isempty(s.(f))
    error('step63:MissingField','Required field %s is missing.',f);
end
end
