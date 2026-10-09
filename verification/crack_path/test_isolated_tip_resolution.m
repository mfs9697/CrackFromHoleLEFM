function Report=test_isolated_tip_resolution(varargin)
% Small deterministic mesh fixtures and metadata checks. No physical FEM solve.
ip=inputParser;addParameter(ip,'WorkDir',tempname,@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});work=char(ip.Results.WorkDir);mkdir(work);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(genpath(root));
count=0;da=.004;
for scale=[.5,1,2]
    assert(resolve_crack_exterior_scale(scale)==scale);
    assert(resolve_crack_exterior_scale(scale,[])==scale);
    assert(resolve_crack_exterior_scale(scale,1)==1);count=count+3;
end
expect(@()resolve_crack_exterior_scale(2,NaN),'crackmesh:ExteriorScale');count=count+1;
reference=build_stage2_scaled_audited_core([0,0],[1,0],da,'Scale',1);
coarse=build_stage2_scaled_audited_core([0,0],[1,0],da,'Scale',2);
assert(size(reference.local.connect3,1)==12678&&size(coarse.local.connect3,1)==3318);
assert(coarse.hTip==2*reference.hTip&&coarse.rCore==reference.rCore);count=count+2;
r=linspace(.75*da,20*da,1001);cal=struct('transitionLength_m',da,'farSlope',.10);
base=reference.design.hBase_m;rp=reference.rCore;slope=reference.design.slope;
hReference=crack_exterior_size_law(r,rp,base+slope*rp,.625*da,1,slope,cal);
es=resolve_crack_exterior_scale(2,1);
hIsolated=crack_exterior_size_law(r,rp,es*(base+slope*rp),.625*da,es,slope,cal);
assert(isequal(hReference,hIsolated));count=count+1;
for scale=[.5,1,2]
    es=resolve_crack_exterior_scale(scale);
    actual=crack_exterior_size_law(r,rp,es*(base+slope*rp),.625*da,es,slope,cal);
    % Original implementation retained here as an independent numeric oracle.
    t=max(0,r-rp);L=cal.transitionLength_m;hRp=scale*(base+slope*rp);cap=.625*da-hRp;
    z=t/L;logcosh=z+log1p(exp(-2*z))-log(2);
    increment=scale*(slope*t+(cal.farSlope-slope)*L*logcosh);
    historical=hRp+cap*tanh(increment/cap);
    assert(isequal(actual,historical));count=count+1;
end
% Tiny cracked-square source: both exterior builders must be bitwise identical
% for an omitted scale and explicit ExteriorScale=CoreScale.
X=da*[-1,-1;1,-1;1,1;-1,1;-1,0;0,0;-1,0];
T=[1,2,6;2,3,6;3,4,6;4,5,6;7,1,6];
cr=struct('Pmid',da*[-1,0;0,0],'upperNodes',5,'lowerNodes',7,'tipNode',6);
design=coarse.design;design.transitionLength_m=da;design.farCap_m=.625*da;
design.exteriorCalibration=struct('verbose',false);
explicit=design;explicit.exteriorScale=2;
for name={'build_stage2_scaled_audited_exterior','build_stage3c_polyline_exterior'}
    builder=str2func(name{1});
    if contains(name{1},'polyline'),cr.Pmid=da*[-1,0;-.875,0;0,0];end
    [p,t,info,ids,edges]=builder(X,T,cr,coarse.local.coord3,coarse.local.connect3,rp,design);
    [p2,t2,info2,ids2,edges2]=builder(X,T,cr,coarse.local.coord3,coarse.local.connect3,rp,explicit);
    assert(isequal(p,p2)&&isequal(t,t2)&&isequaln(info,info2)&& ...
        isequal(ids,ids2)&&isequal(edges,edges2));count=count+1;
    isolated=design;isolated.exteriorScale=1;
    [~,~,iso]=builder(X,T,cr,coarse.local.coord3,coarse.local.connect3,rp,isolated);
    assert(iso.exteriorScale==1&&iso.hRp_m==base+slope*rp);
    radii=[rp;iso.ringTable.radius_m(1:end-1)];
    expected=crack_exterior_size_law(radii,rp,base+slope*rp,.625*da,1,slope,cal);
    assert(isequal(iso.ringTable.targetH_m,expected));count=count+1;
end
legacy=struct('farCapOverA0',.625,'transitionOverA0',1,'calibrationOverride',struct());
candidate=struct('coreMeshControls',struct('scale',2),'exteriorMeshControls',legacy);
wanted=struct('coreScale',2,'exteriorScale',2,'farCapOverIncrement',.625, ...
    'transitionOverIncrement',1,'calibrationOverride',struct());
assert(crack_candidate_mesh_controls_match(candidate,wanted));count=count+1;
wanted.exteriorScale=1;
assert(~crack_candidate_mesh_controls_match(candidate,wanted));count=count+1;
candidate.exteriorMeshControls.exteriorScale=1;
assert(crack_candidate_mesh_controls_match(candidate,wanted));count=count+1;
saved=struct('meta',struct('coreScale',2,'exteriorMeshControls', ...
    crack_exterior_identity(legacy,2)));
expect(@()assert_crack_checkpoint_exterior(saved,candidate,'test:Checkpoint'),'test:Checkpoint');count=count+1;
saved.meta.exteriorMeshControls.exteriorScale=1;
assert_crack_checkpoint_exterior(saved,candidate,'test:Checkpoint');count=count+1;
saved.meta=rmfield(saved.meta,'exteriorMeshControls');
expect(@()assert_crack_checkpoint_exterior(saved,candidate,'test:Checkpoint'),'test:Checkpoint');count=count+1;
candidate.exteriorMeshControls=legacy;
assert_crack_checkpoint_exterior(saved,candidate,'test:Checkpoint');count=count+1;
% Public resume rejects a different exterior before examining path history.
z=load(fullfile(root,'paper','data','accepted_stage1_source.mat'),'R0');
state=struct('vertices',zeros(4,2),'thetaDeg',zeros(3,1), ...
    'completedPhysicalSegments',2,'nextThetaDeg',0,'exteriorMeshControls',wanted);
state.exteriorMeshControls.exteriorScale=2;
expect(@()run_incremental_crack_path('FrozenState',z.R0,'CoreScale',2,'ExteriorScale',1, ...
    'ResumeState',state,'RegressionGates',false,'AllowPhysicalSolves',false, ...
    'OutputDir',fullfile(work,'incompatible_resume')),'pathrun:ResumeMeshFamily');count=count+1;
state.exteriorMeshControls=rmfield(state.exteriorMeshControls,'exteriorScale');
expect(@()run_incremental_crack_path('FrozenState',z.R0,'CoreScale',2,'ExteriorScale',1, ...
    'ResumeState',state,'RegressionGates',false,'AllowPhysicalSolves',false, ...
    'OutputDir',fullfile(work,'incompatible_legacy_resume')),'pathrun:ResumeMeshFamily');count=count+1;
metadataFree=rmfield(state,'exteriorMeshControls');
expect(@()run_incremental_crack_path('FrozenState',z.R0,'CoreScale',1,'ExteriorScale',2, ...
    'ResumeState',metadataFree,'RegressionGates',false,'AllowPhysicalSolves',false, ...
    'OutputDir',fullfile(work,'missing_exterior_resume')),'pathrun:ResumeMeshFamilyMissing');count=count+1;
% P1 compact shortcut rejects the scale before any checkpoint file access.
fp=struct('scale',2,'hTip_m',coarse.hTip,'rInner_m',.1*da,'rOuter_m',.65*da, ...
    'rCore_m',rp,'expectedCoreT3',3318,'expectedEDIElements',2976, ...
    'expectedOptimizedSupport',2700,'expectedNativeSamples',[19;28;23;18]);
candidate.coreMeshControls=fp;candidate.exteriorMeshControls.exteriorScale=1;
cached=struct('pass',true,'gates',struct('pass',true), ...
    'coreMeshControls',fp,'exteriorMeshControls',legacy);
expect(@()validate_increment_study_p1_cache(cached,candidate,z.R0,'unused.mat'), ...
    'incstudy:P1ReuseMismatch');count=count+1;
% The driver must reject an existing directory before any qualification/solve.
expect(@()main_isolated_tip_resolution_study('OutputDir',work), ...
    'isolatedtip:OutputExists');count=count+1;
archiveChecks=0;
archive=fullfile(root,'verification','crack_path','increment_sensitivity_coarse','da_4mm');
if exist(fullfile(archive,'p1_seed','P1_candidate.mat'),'file')==2
    a=load(fullfile(archive,'p1_seed','P1_candidate.mat'),'candidate');
    a.candidate.exteriorMeshControls.exteriorScale=1;
    expect(@()main_stage2_theta0_physical_solve('FrozenState',z.R0,'Candidate',a.candidate, ...
        'AlternativeQualifiedCandidate',true,'AllowSolve',false, ...
        'CheckpointFile',fullfile(archive,'p1_seed','P1_physical_solved.mat'), ...
        'SaveFile',fullfile(work,'never_written_P1.mat')),'stage2phys:CheckpointMismatch');
    archiveChecks=archiveChecks+1;
    a=load(fullfile(archive,'trajectory','step_003_candidate.mat'),'candidate');
    a.candidate.exteriorMeshControls.exteriorScale=1;
    expect(@()solve_incremental_crack_tip(a.candidate,'FrozenState',z.R0,'AllowSolve',false, ...
        'CheckpointFile',fullfile(archive,'trajectory','step_003_physical_solved.mat'), ...
        'SaveFile',fullfile(work,'never_written_P3.mat')),'pathsolve:CheckpointMismatch');
    archiveChecks=archiveChecks+1;
end
Report=struct('pass',true,'checks',count,'physicalSolves',0,'workDir',work);
Report.archiveChecks=archiveChecks;
save(fullfile(work,'test_isolated_tip_resolution.mat'),'Report','-v7');
fprintf('ISOLATED TIP RESOLUTION TESTS: PASS (%d portable, %d archive checks, no physical solves)\n',count,archiveChecks);
end

function expect(fn,id)
try
    fn();
catch ME
    assert(strcmp(ME.identifier,id),'test:WrongError','Expected %s, got %s: %s',id,ME.identifier,ME.message);
    return
end
error('test:MissingError','Expected error %s.',id);
end
