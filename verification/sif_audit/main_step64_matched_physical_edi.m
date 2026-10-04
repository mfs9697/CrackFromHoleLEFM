function R64=main_step64_matched_physical_edi(varargin)
%MAIN_STEP64_MATCHED_PHYSICAL_EDI
% Exactly ONE matched physical interaction-EDI extraction on the already
% solved Step63 calibrated asymmetric FEM field.
%
% NO FEM solve, NO mesh generation/remeshing, NO radius sweep.
%
% Fixed domain:
%   r_inner = 0.8 mm
%   r_outer = 5.2 mm = 0.65*a0
%
% Extraction:
%   16-point Dunavant quadrature
%   WeightFunction = 'fe_nodal'
%
% The result is compared with the Step63 native-COD fits on the SAME
% physical displacement field. A completed matching Step64 result is
% reused rather than recomputed.
%
% Usage:
%   addpath(genpath(pwd));
%   R64=main_step64_matched_physical_edi();
%   disp(R64.EDI);
%   disp(R64.CODcomparison);
%
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');

ip=inputParser;
addParameter(ip,'Step63Checkpoint', ...
    fullfile(vdir,'step63_calibrated_asymmetric_physical_solved.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'Step63CODFile', ...
    fullfile(vdir,'step63_calibrated_asymmetric_cod_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'SaveFile', ...
    fullfile(vdir,'step64_matched_physical_edi_small_data.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;

addpath(genpath(root));
assert_step64_branch(root);

cp=char(opt.Step63Checkpoint);
codFile=char(opt.Step63CODFile);
saveFile=char(opt.SaveFile);

if exist(cp,'file')~=2
    error('step64:MissingStep63Checkpoint', ...
        'Missing saved physical Step63 checkpoint: %s',cp);
end
if exist(codFile,'file')~=2
    error('step64:MissingStep63COD', ...
        'Missing Step63 COD result: %s',codFile);
end

s=load(cp,'mesh','U','mat','crack','a0','meta');
c=load(codFile,'O63');
for f={'mesh','U','mat','crack','a0','meta'},need(s,f{1});end
if ~isfield(c,'O63'),error('step64:BadCODFile','Expected O63 in COD file.');end
O63=c.O63;

% -------------------------------------------------------------------------
% Strong provenance: exact calibrated physical field from Step63.
% -------------------------------------------------------------------------
m=s.meta;
requiredMeta={'stage','selectedCandidate','nT3','nT6','a0', ...
    'loadType','sig0','anchorMode','candidateSHA256', ...
    'calibrationSHA256','maxNeighborRatio','physicalEDIperformed', ...
    'meshRegenerated'};
for k=1:numel(requiredMeta),need(m,requiredMeta{k});end
if ~strcmp(m.stage,'step63_calibrated_asymmetric_physical') || ...
        ~strcmp(m.selectedCandidate,'C03') || ...
        m.nT3~=32980 || m.nT6~=66854 || ...
        abs(s.a0-.008)>1e-12 || abs(m.a0-s.a0)>1e-12 || ...
        ~strcmp(m.loadType,'remote_tension_y') || abs(m.sig0-1)>1e-14 || ...
        ~strcmp(m.anchorMode,'minimal') || ...
        abs(m.maxNeighborRatio-1.79678451)>5e-7 || ...
        logical(m.physicalEDIperformed) || logical(m.meshRegenerated)
    error('step64:WrongStep63Field', ...
        'Checkpoint is not the exact unswept physical C03 Step63 field.');
end
if size(s.mesh.connect3,1)~=m.nT3 || size(s.mesh.coord,1)~=m.nT6 || ...
        numel(s.U)~=2*m.nT6 || any(~isfinite(s.U))
    error('step64:CheckpointDimensions', ...
        'Step63 checkpoint mesh/U dimensions are inconsistent.');
end
if ~isfield(O63,'physicalEDIperformed') || O63.physicalEDIperformed || ...
        ~isfield(O63,'fitTable') || ~istable(O63.fitTable) || ...
        height(O63.fitTable)~=8 || ...
        ~isfield(O63,'rawTable') || ~istable(O63.rawTable) || ...
        ~strcmp(O63.candidateSHA256,m.candidateSHA256) || ...
        ~strcmp(O63.calibrationSHA256,m.calibrationSHA256)
    error('step64:CODProvenance', ...
        'Step63 COD result is missing or belongs to another field.');
end

ri=.0008;
ro=.65*s.a0;
if abs(ro-.0052)>1e-14
    error('step64:MatchedRadius','Expected r_outer=5.2 mm.');
end
hTip=tip_edge_median(s.mesh.coord3,s.mesh.connect3,s.crack.Pmid(end,:));
if ri<2*hTip || ri>=ro
    error('step64:AnnulusInadmissible', ...
        'Matched annulus is not admissible for the Step63 tip scale.');
end

% Hashes make a previously completed Step64 result safely reusable.
cpHash=sha256(cp);
codHash=sha256(codFile);
if exist(saveFile,'file')==2
    old=load(saveFile,'R64');
    if isfield(old,'R64') && ...
            isfield(old.R64,'checkpointSHA256') && ...
            isfield(old.R64,'CODSHA256') && ...
            strcmp(old.R64.checkpointSHA256,cpHash) && ...
            strcmp(old.R64.CODSHA256,codHash) && ...
            isfield(old.R64,'singleEDIDomain') && old.R64.singleEDIDomain && ...
            abs(old.R64.rInner-ri)<1e-14 && abs(old.R64.rOuter-ro)<1e-14
        fprintf('\nSTEP64: reusing completed matched physical EDI; NO repeated integration.\n');
        R64=old.R64;
        disp(R64.EDI);
        disp(R64.CODcomparison);
        return
    else
        error('step64:StaleSavedResult', ...
            'Existing Step64 result belongs to another Step63 field.');
    end
end

fprintf('\n============================================================\n');
fprintf('STEP 64: SINGLE MATCHED PHYSICAL EDI ON SAVED STEP63 FIELD\n');
fprintf('============================================================\n');
fprintf('  Saved physical mesh: T3=%d, T6=%d; no solve/remesh.\n',m.nT3,m.nT6);
fprintf('  Exact annulus: r_inner=%.12g m, r_outer=%.12g m (0.65*a0).\n',ri,ro);
fprintf('  Tip median edge=%.12g m; 2*hTip=%.12g m.\n',hTip,2*hTip);
fprintf('  ONE 16-point FE-nodal-q interaction EDI. NO radius sweep.\n');

mat=s.mat;
if ~isfield(mat,'Dmat') && isfield(mat,'D'),mat.Dmat=mat.D;end
[KI,KII]=SIF_LEFM_interaction_EDI( ...
    s.mesh,s.U,s.crack.Pmid,mat, ...
    struct('r_inner',ri,'r_outer',ro), ...
    'UsePlaneStrain',mat.ps==1, ...
    'Verbose',false, ...
    'WeightFunction','fe_nodal', ...
    'QuadratureRule',16, ...
    'StoreGPDiagnostics',false);

q=KII/KI;
if ~all(isfinite([KI,KII,q])) || KI<=0
    error('step64:NonfiniteEDI','Physical EDI returned invalid SIFs.');
end

EDI=table(KI,KII,q,ri,ro,ro/s.a0, ...
    'VariableNames',{'KI','KII','ratio','r_inner','r_outer','r_outer_over_a0'});

F=O63.fitTable;
delta=F.ratio_COD-q;
relPct=100*delta/q;
absRelPct=abs(relPct);
CODcomparison=table( ...
    F.lower_r_over_a0,F.upper_r_over_a0,F.degree,F.n_native, ...
    F.KI_COD,F.KII_COD,F.ratio_COD, ...
    repmat(q,height(F),1),delta,relPct,absRelPct, ...
    'VariableNames',{'lower_r_over_a0','upper_r_over_a0','degree','n_native', ...
    'KI_COD','KII_COD','ratio_COD','ratio_EDI','delta_COD_minus_EDI', ...
    'relative_gap_pct','absolute_relative_gap_pct'});

codMin=min(F.ratio_COD);
codMax=max(F.ratio_COD);
codMean=mean(F.ratio_COD);
codMedian=median(F.ratio_COD);
ediInsideCODRange=(q>=codMin && q<=codMax);
meanGapPct=100*(codMean-q)/q;
medianGapPct=100*(codMedian-q)/q;

Summary=table(KI,KII,q,codMin,codMax,codMean,codMedian, ...
    ediInsideCODRange,meanGapPct,medianGapPct, ...
    'VariableNames',{'KI_EDI','KII_EDI','ratio_EDI', ...
    'COD_ratio_min','COD_ratio_max','COD_ratio_mean','COD_ratio_median', ...
    'EDI_inside_COD_fit_range','COD_mean_gap_pct','COD_median_gap_pct'});

fprintf('\nPHYSICAL STEP63 COD vs SINGLE MATCHED EDI\n');
disp(EDI);
disp(Summary);
disp(CODcomparison);

R64=struct( ...
    'checkpointPath',cp, ...
    'checkpointSHA256',cpHash, ...
    'CODFile',codFile, ...
    'CODSHA256',codHash, ...
    'EDI',EDI, ...
    'Summary',Summary, ...
    'CODcomparison',CODcomparison, ...
    'rInner',ri,'rOuter',ro,'hTip',hTip, ...
    'singleEDIDomain',true, ...
    'noNewFEM',true,'noMesh',true,'noRemesh',true,'noRadiusSweep',true, ...
    'interpretation',['One physical EDI extraction on the same calibrated ', ...
      'Step63 field used by COD. Agreement is cross-extractor evidence on ', ...
      'one physical mesh, not yet a mesh-convergence proof.']);

save(saveFile,'R64','-v7');
fprintf('  Compact Step64 result saved: %s\n',saveFile);
fprintf('STEP64 complete: one physical EDI, zero FEM solves.\n');
end

% =========================================================================
function h=tip_edge_median(P,T,tip)
T=T(:,1:3);
r=hypot(P(:,1)-tip(1),P(:,2)-tip(2));
tol=max(1e-12,1e-8*max(1,max(abs(P(:)))));
ids=find(r<=min(r)+tol);
tri=T(any(ismember(T,ids),2),:);
p1=P(tri(:,1),:);p2=P(tri(:,2),:);p3=P(tri(:,3),:);
e=[hypot(p1(:,1)-p2(:,1),p1(:,2)-p2(:,2)); ...
   hypot(p2(:,1)-p3(:,1),p2(:,2)-p3(:,2)); ...
   hypot(p3(:,1)-p1(:,1),p3(:,2)-p1(:,2))];
e=e(isfinite(e)&e>tol);
if isempty(e),error('step64:TipEdges','No tip-adjacent edges found.');end
h=median(e);
end

function h=sha256(path)
f=fopen(path,'rb');
if f<0,error('step64:HashOpen','Cannot open %s.',path);end
guard=onCleanup(@()fclose(f));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(f)
    bytes=fread(f,1024*1024,'*uint8');
    md.update(typecast(bytes,'int8'));
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
clear guard
end

function assert_step64_branch(root)
[status,b]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0&&strcmp(strtrim(b),'audit/step64-matched-physical-edi'), ...
    'step64:Branch','Run Step64 only on audit/step64-matched-physical-edi.');
end

function need(s,f)
if ~isstruct(s)||~isfield(s,f)||isempty(s.(f))
    error('step64:MissingField','Required field %s is missing.',f);
end
end
