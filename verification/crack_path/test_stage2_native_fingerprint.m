function Report=test_stage2_native_fingerprint(varargin)
% Investigator-local regression of the isolated P1 failure; no physical solve.
ip=inputParser;
addParameter(ip,'StudyDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'WorkDir',tempname,@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(genpath(root));
study=char(ip.Results.StudyDir);work=char(ip.Results.WorkDir);mkdir(work);
assert(~isempty(study),'nativefp:StudyRequired','Pass the saved isolated study directory.');
z=load(fullfile(root,'paper','data','accepted_stage1_source.mat'),'R0');R0=z.R0;
z=load(fullfile(study,'independent_core_2_exterior_1','p1_seed','P1_candidate.mat'),'candidate');c=z.candidate;
mesh=make_mesh(c);fp=stage2_native_sampling_fingerprint(c,mesh,c.mat);
assert(fp.nUpper==80&&fp.nLower==80&&fp.coreFaceNodes==60&& ...
    isequal(fp.nativePoints,[19;28;23;18]));count=1;
ordered=c;ordered.crack.upperNodes=flipud(ordered.crack.upperNodes(:));
ordered.crack.lowerNodes=flipud(ordered.crack.lowerNodes(:));
fp2=stage2_native_sampling_fingerprint(ordered,mesh,c.mat);assert(isequaln(fp,fp2));count=count+1;
% Reordering every T3 node changes T6 IDs but not the geometric fingerprint.
renumbered=c;n=size(c.p,1);renumbered.p=flipud(c.p);renumbered.t=n+1-c.t;
renumbered.crack.upperNodes=n+1-c.crack.upperNodes;
renumbered.crack.lowerNodes=n+1-c.crack.lowerNodes;renumbered.crack.tipNode=n+1-c.crack.tipNode;
fp2=stage2_native_sampling_fingerprint(renumbered,make_mesh(renumbered),c.mat);
assert(fp2.nUpper==fp.nUpper&&isequal(fp2.nativePoints,fp.nativePoints)&& ...
    max(abs(fp2.fullFaceR_m-fp.fullFaceR_m))<=1e-12);count=count+1;
% Shift paired midsides together: counts/grid pairing remain unchanged, but
% the independent T3 face-chain check must detect the altered T6 positions.
[~,~,diag]=native_COD_audit(mesh,zeros(2*size(mesh.coord,1),1),c.mat,c.crack,8,true);
direction=diff(c.crack.Pmid);direction=direction/norm(direction);
r=-(mesh.coord-c.crack.Pmid(end,:))*direction(:);
u=find(diag.faceSide==1&r>.0004&r<.0008& ...
    (1:size(mesh.coord,1))'>size(c.p,1),1);
lower=find(diag.faceSide==-1);[~,j]=min(abs(r(lower)-r(u)));l=lower(j);
bad=mesh;bad.coord([u,l],:)=bad.coord([u,l],:)+1e-6*direction;
expect(@()stage2_native_sampling_fingerprint(c,bad,c.mat),'stage2phys:NativeSampling');count=count+1;
% T3 and T6 agree after a paired core-face shift, so the canonical core check
% is needed independently of the complete-face representation.
bad=c;if isfield(bad,'nativeSamplingFingerprint'),bad=rmfield(bad,'nativeSamplingFingerprint');end
u=c.crack.upperNodes;r3=-(c.p(u,:)-c.crack.Pmid(end,:))*direction(:);
u=u(find(r3>.0004&r3<.0008,1));lower=c.crack.lowerNodes;
[~,j]=min(vecnorm(c.p(lower,:)-c.p(u,:),2,2));l=lower(j);
bad.p([u,l],:)=bad.p([u,l],:)+1e-6*direction;
expect(@()stage2_native_sampling_fingerprint(bad,make_mesh(bad),c.mat),'stage2phys:NativeSampling');count=count+1;
bad=c;bad.nativeSamplingFingerprint=fp;bad.nativeSamplingFingerprint.nativePoints(1)=20;
expect(@()stage2_native_sampling_fingerprint(bad,mesh,c.mat),'stage2phys:NativeSampling');count=count+1;
% The public solver now gets past the fingerprint and reaches its explicit
% solve guard, proving the saved qualified candidate is directly usable.
expect(@()main_stage2_theta0_physical_solve('FrozenState',R0,'Candidate',c, ...
    'AlternativeQualifiedCandidate',true,'AllowSolve',false, ...
    'CheckpointFile',fullfile(work,'absent.mat'),'SaveFile',fullfile(work,'unused.mat')), ...
    'stage2phys:ExplicitSolveApprovalRequired');count=count+1;
legacy=fullfile(root,'verification','crack_path','tip_2h0_independent_run', ...
    'p1_seed','tip2h0_P1_candidate.mat');
if exist(legacy,'file')==2
    z=load(legacy,'candidate');old=z.candidate;
    legacyFP=stage2_native_sampling_fingerprint(old,make_mesh(old),old.mat);
    assert(legacyFP.nUpper==70&&legacyFP.coreFaceNodes==60&& ...
        isequal(legacyFP.nativePoints,fp.nativePoints));count=count+1;
end
cacheChecks=0;
for scale=[2,.5]
    if scale==2,label='core_2_exterior_1';else,label='core_0p5_exterior_1';end
    for k=[17,21,22,23]
        base=fullfile(study,'fixed_geometry',label);
        z=load(fullfile(base,sprintf('step_%03d_candidate.mat',k)),'candidate');candidate=z.candidate;
        z=load(fullfile(base,sprintf('step_%03d_physical_small.mat',k)),'R');R=z.R;
        cp=fullfile(base,sprintf('step_%03d_physical_solved.mat',k));
        validate_isolated_fixed_result_cache(R,candidate,R0,cp);cacheChecks=cacheChecks+1;
        if k==17&&scale==2
            mutant=R;mutant.exteriorMeshControls.exteriorScale=2;
            expect(@()validate_isolated_fixed_result_cache(mutant,candidate,R0,cp), ...
                'isolatedtip:FixedReuseMismatch');count=count+1;
        end
    end
end
Report=struct('pass',true,'checks',count,'fixedCacheChecks',cacheChecks,'physicalSolves',0);
save(fullfile(work,'native_fingerprint_tests.mat'),'Report','fp','-v7');
fprintf('NATIVE FINGERPRINT REGRESSION PASS: %d checks, %d accepted fixed-cache checks; zero solves.\n',count,cacheChecks);
end

function mesh=make_mesh(c)
[p6,t6]=T3toT6_fast(c.p,c.t);
mesh=struct('coord3',c.p,'connect3',c.t,'coord',p6,'connect',t6);
end

function expect(fn,id)
try
    fn();
catch ME
    assert(strcmp(ME.identifier,id),'nativefp:WrongError','Expected %s, got %s: %s',id,ME.identifier,ME.message);
    return
end
error('nativefp:MissingError','Expected %s.',id);
end
