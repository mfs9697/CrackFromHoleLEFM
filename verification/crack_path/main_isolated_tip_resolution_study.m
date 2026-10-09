function D=main_isolated_tip_resolution_study(varargin)
% Fixed reference P17/P21/P22/P23 at core scales 2 and .5, exterior scale 1;
% then an independent core=2/exterior=1 trajectory, seeded by its own P1.
% Default: full mesh/prescribed-field qualification, no physical solves.
% main_isolated_tip_resolution_study('AllowPhysicalSolves',true)
% Fresh invocations require a new directory. Explicit ResumeStudyDir reuses
% qualified candidates and accepted fields from that study after identity checks.
ip=inputParser;
addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'FastEDI',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Plot',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'ResumeStudyDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'ReferenceEvidenceFile','',@(x)ischar(x)||isstring(x));
addParameter(ip,'FrozenStateFile','',@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});opt=ip.Results;
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(genpath(root));
evidence=local_file(root,opt.ReferenceEvidenceFile, ...
    fullfile(root,'paper','data','evidence_exact.mat'));
frozen=local_file(root,opt.FrozenStateFile, ...
    fullfile(root,'paper','data','accepted_stage1_source.mat'));
z=load(evidence,'E','T','Q');e=z.E;reference=z.T;referenceQualification=z.Q;
z=load(frozen,'R0');R0=z.R0;
assert(R0.stage1Pass&&logical(R0.summary.stage1_pass), ...
    'isolatedtip:FrozenState','A passed exact Stage-I state is required.');
increment=R0.summary.a0_reserved_m;
mouth=[R0.summary.x_star_m,R0.summary.y_star_m];
assert(norm(mouth-e.vertices_m(1,:))<=2e-12, ...
    'isolatedtip:ReferenceMismatch','Reference and frozen Stage-I mouths differ.');
segments=[17,21,22,23];scales=[2,.5];
assert(all(ismember(segments,reference.segment))&&size(e.vertices_m,1)>=24, ...
    'isolatedtip:ReferenceStates','Exact accepted reference P17/P21/P22/P23 required.');
for k=segments
    path=e.vertices_m(1:k+1,:);r=reference(reference.segment==k,:);
    q=referenceQualification(referenceQualification.segment==k,:);
    assert(height(r)==1&&height(q)==1&&logical(q.pass)&& ...
        abs(r.crack_length_mm-1e3*sum(vecnorm(diff(path),2,2)))<=2e-10&& ...
        all(abs(vecnorm(diff(path),2,2)-increment)<=2e-12), ...
        'isolatedtip:ReferenceGeometry','Reference P%d is incompatible with frozen increment.',k);
end
defaultDir=fullfile(root,'verification','crack_path', ...
    ['isolated_tip_resolution_study_',char(datetime('now','Format','yyyyMMdd''T''HHmmssSSS'))]);
resuming=~isempty(char(opt.ResumeStudyDir));
if resuming
    assert(isempty(char(opt.OutputDir)),'isolatedtip:OutputConflict', ...
        'Use ResumeStudyDir alone when continuing an existing study.');
    out=local_file(root,opt.ResumeStudyDir,'');
    z=load(fullfile(out,'isolated_tip_resolution_study.mat'),'D');prior=z.D;
    assert(strcmp(prior.profile,'isolated_tip_resolution')&& ...
        isequal(prior.coreScales,scales)&&prior.exteriorScale==1&& ...
        isequal(prior.segments,segments)&&abs(prior.increment_m-increment)<=1e-14, ...
        'isolatedtip:ResumeIdentity','Saved study uses different controls.');
else
    out=local_file(root,opt.OutputDir,defaultDir);
    assert(exist(out,'dir')~=7&&exist(out,'file')~=2, ...
        'isolatedtip:OutputExists','Choose a new directory; study outputs cannot overwrite archives: %s',out);
    mkdir(out);
end
D=struct('profile','isolated_tip_resolution','coreScales',scales,'exteriorScale',1, ...
    'segments',segments,'increment_m',increment,'referenceEvidenceFile',evidence, ...
    'frozenStateFile',frozen,'outputDir',out,'physicalSolvesEnabled',opt.AllowPhysicalSolves, ...
    'fixedQualification',table(),'fixedPhysical',table(),'independentPath',[], ...
    'resumed',resuming,'reusedFixedPhysicalStates',0,'reusedQualifiedCandidates',0);
for scale=scales
    if scale==2,label='core_2_exterior_1';else,label='core_0p5_exterior_1';end
    fixedDir=fullfile(out,'fixed_geometry',label);mkdir(fixedDir);
    for k=segments
        candidateFile=fullfile(fixedDir,sprintf('step_%03d_candidate.mat',k));
        qualFile=fullfile(fixedDir,sprintf('step_%03d_qualification_small.mat',k));
        if resuming&&exist(candidateFile,'file')==2
            Q=local_load_qualified(candidateFile,qualFile,R0,scale,e.vertices_m(1:k+1,:));
            D.reusedQualifiedCandidates=D.reusedQualifiedCandidates+1;
        else
            Q=qualify_incremental_crack_candidate(e.vertices_m(1:k+1,:), ...
                'FrozenState',R0,'CoreScale',scale,'ExteriorScale',1, ...
                'RunSynthetic',true,'FastEDI',opt.FastEDI,'Plot',opt.Plot, ...
                'SaveCandidate',true,'CandidateFile',candidateFile,'SaveCompact',true, ...
                'CompactFile',fullfile(fixedDir,sprintf('step_%03d_qualification_small.mat',k)));
        end
        assert(Q.pass,'isolatedtip:Qualification','Fixed P%d core %.6g qualification failed.',k,scale);
        row=Q.summary;row.reference_segment=k;
        D.fixedQualification=[D.fixedQualification;row];
        if opt.AllowPhysicalSolves
            cp=fullfile(fixedDir,sprintf('step_%03d_physical_solved.mat',k));
            sf=fullfile(fixedDir,sprintf('step_%03d_physical_small.mat',k));
            if resuming&&exist(sf,'file')==2
                z=load(sf,'R');R=z.R;
                validate_isolated_fixed_result_cache(R,Q.candidate,R0,cp);
                D.reusedFixedPhysicalStates=D.reusedFixedPhysicalStates+1;
                fprintf('Reusing accepted fixed P%d core %.6g; no solve or postprocessing.\n',k,scale);
            else
                R=solve_incremental_crack_tip(Q.candidate,'FrozenState',R0,'AllowSolve',true, ...
                    'FastEDI',opt.FastEDI,'CheckpointFile',cp,'SaveFile',sf);
            end
            assert(R.pass,'isolatedtip:Physical','Fixed P%d core %.6g physical gates failed.',k,scale);
            row=R.summary;ref=reference(reference.segment==k,:);
            row.reference_KI_unit=ref.KI_unit;row.reference_KII_unit=ref.KII_unit;
            row.reference_turn_deg=ref.delta_theta_next_deg;
            row.KI_difference=row.KI_unit-ref.KI_unit;
            row.KII_difference=row.KII_unit-ref.KII_unit;
            row.turn_difference_deg=row.delta_theta_next_MTS_deg-ref.delta_theta_next_deg;
            D.fixedPhysical=[D.fixedPhysical;row];
        end
        local_save(out,D);
        clear Q R
    end
end
p1Dir=fullfile(out,'independent_core_2_exterior_1','p1_seed');mkdir(p1Dir);
p1Candidate=fullfile(p1Dir,'P1_candidate.mat');p1Qual=fullfile(p1Dir,'P1_qualification_small.mat');
if resuming&&exist(p1Candidate,'file')==2
    n=[R0.summary.nmat_x,R0.summary.nmat_y];n=n/norm(n);
    F=local_load_qualified(p1Candidate,p1Qual,R0,2,[mouth;mouth+increment*n]);
    D.reusedQualifiedCandidates=D.reusedQualifiedCandidates+1;
else
    F=main_stage2_embed_scaled_core_full_domain_theta0('FrozenState',R0, ...
        'CoreScale',2,'ExteriorScale',1,'Plot',opt.Plot,'ExteriorVerbose',true, ...
        'SaveCandidate',true,'CandidateFile',fullfile(p1Dir,'P1_candidate.mat'), ...
        'SaveCompact',true,'CompactFile',fullfile(p1Dir,'P1_qualification_small.mat'));
end
assert(F.pass,'isolatedtip:P1Qualification','Independent isolated P1 failed qualification.');
D.P1Qualification=F.summary;local_save(out,D);
if ~opt.AllowPhysicalSolves
    fprintf('Isolated-tip qualification complete. No physical solves performed.\nOutputs: %s\n',out);
    return
end
p1Checkpoint=fullfile(p1Dir,'P1_physical_solved.mat');p1Small=fullfile(p1Dir,'P1_physical_small.mat');
if resuming&&exist(p1Small,'file')==2
    z=load(p1Small,'R');R1=z.R;
    validate_increment_study_p1_cache(R1,F.candidate,R0,p1Checkpoint);
    fprintf('Reusing accepted isolated P1; no solve or postprocessing.\n');
else
    R1=main_stage2_theta0_physical_solve('FrozenState',R0,'Candidate',F.candidate, ...
        'AlternativeQualifiedCandidate',true,'AllowSolve',true, ...
        'CheckpointFile',p1Checkpoint,'SaveFile',p1Small);
end
assert(R1.pass,'isolatedtip:P1Physical','Independent isolated P1 failed physical gates.');
D.P1Physical=R1.summary;local_save(out,D);clear F
% Existing alternative-family mode skips only historical reference-value anchors.
% Qualification, synthetic replay, solver/residual, EDI, COD and MTS gates
% and their numerical tolerances remain in the unchanged production routines.
pathDir=fullfile(out,'independent_core_2_exterior_1','trajectory');resumeOptions={};
if resuming&&exist(fullfile(pathDir,'path_run_state.mat'),'file')==2
    resumeOptions={'ResumeStateFile',fullfile(pathDir,'path_run_state.mat'),'ResumeSourceDir',pathDir};
end
D.independentPath=run_incremental_crack_path('FrozenState',R0,'MaxSegments',23, ...
    'CoreScale',2,'ExteriorScale',1,'AllowPhysicalSolves',true,'RunSynthetic',true, ...
    'FastEDI',opt.FastEDI,'RegressionGates',false,'StopAtCoreClearance',true, ...
    'PlotEachStep',opt.Plot,'ReuseCandidates',true, ...
    'MeshFamilyLabel','isolated_core_2_exterior_1', ...
    'SeedKI',R1.EDI.KI_unit,'SeedKII',R1.EDI.KII_unit, ...
    'OutputDir',pathDir,resumeOptions{:});
local_save(out,D);
assert(strcmp(D.independentPath.stopReason,'max_segments_reached')&& ...
    height(D.independentPath.stepTable)==23,'isolatedtip:PathStopped', ...
    'Independent path stopped before P23 under unchanged gates; inspect saved results.');
end

function path=local_file(root,given,default)
path=char(given);
if isempty(path),path=default;
elseif ~(startsWith(path,filesep)||startsWith(path,'\\')|| ...
        ~isempty(regexp(path,'^[A-Za-z]:[\\/]','once')))
    path=fullfile(root,path);
end
end

function local_save(out,D)
save(fullfile(out,'isolated_tip_resolution_study.mat'),'D','-v7');
if ~isempty(D.fixedQualification)
    writetable(D.fixedQualification,fullfile(out,'fixed_qualification.csv'));
end
if ~isempty(D.fixedPhysical)
    writetable(D.fixedPhysical,fullfile(out,'fixed_physical_comparison.csv'));
end
end

function Q=local_load_qualified(candidateFile,qualFile,R0,scale,path)
z=load(candidateFile,'candidate');c=z.candidate;
z=load(qualFile,'Small');Q=z.Small;
requested=struct('coreScale',scale,'exteriorScale',1, ...
    'farCapOverIncrement',.625,'transitionOverIncrement',1,'calibrationOverride',struct());
assert(Q.pass&&all(structfun(@logical,c.gates))&&all(structfun(@logical,c.syntheticGates))&& ...
    crack_candidate_mesh_controls_match(c,requested),'isolatedtip:QualifiedReuseMismatch', ...
    'Saved candidate has different controls or failed qualification.');
if isfield(c,'path'),stored=c.path;else,stored=c.crack.Pmid;end
assert(isequal(size(stored),size(path))&&norm(stored-path,'fro')<=2e-12&& ...
    isequaln(c.mat.E,R0.C.E)&&isequaln(c.mat.nu,R0.C.nu)&&isequaln(c.mat.ps,R0.C.ps), ...
    'isolatedtip:QualifiedReuseMismatch','Saved qualified geometry/physics differs.');
Q.candidate=c;
fprintf('Reusing saved qualified candidate: %s\n',candidateFile);
end
