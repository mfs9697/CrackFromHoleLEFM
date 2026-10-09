function D=main_isolated_tip_resolution_study(varargin)
% Fixed reference P17/P21/P22/P23 at core scales 2 and .5, exterior scale 1;
% then an independent core=2/exterior=1 trajectory, seeded by its own P1.
% Default: full mesh/prescribed-field qualification, no physical solves.
% main_isolated_tip_resolution_study('AllowPhysicalSolves',true)
% Each invocation requires a new output directory, preserving all archives.
ip=inputParser;
addParameter(ip,'AllowPhysicalSolves',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'FastEDI',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Plot',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
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
out=local_file(root,opt.OutputDir,defaultDir);
assert(exist(out,'dir')~=7&&exist(out,'file')~=2, ...
    'isolatedtip:OutputExists','Choose a new directory; study outputs cannot overwrite archives: %s',out);
mkdir(out);
D=struct('profile','isolated_tip_resolution','coreScales',scales,'exteriorScale',1, ...
    'segments',segments,'increment_m',increment,'referenceEvidenceFile',evidence, ...
    'frozenStateFile',frozen,'outputDir',out,'physicalSolvesEnabled',opt.AllowPhysicalSolves, ...
    'fixedQualification',table(),'fixedPhysical',table(),'independentPath',[]);
for scale=scales
    if scale==2,label='core_2_exterior_1';else,label='core_0p5_exterior_1';end
    fixedDir=fullfile(out,'fixed_geometry',label);mkdir(fixedDir);
    for k=segments
        candidateFile=fullfile(fixedDir,sprintf('step_%03d_candidate.mat',k));
        Q=qualify_incremental_crack_candidate(e.vertices_m(1:k+1,:), ...
            'FrozenState',R0,'CoreScale',scale,'ExteriorScale',1, ...
            'RunSynthetic',true,'FastEDI',opt.FastEDI,'Plot',opt.Plot, ...
            'SaveCandidate',true,'CandidateFile',candidateFile,'SaveCompact',true, ...
            'CompactFile',fullfile(fixedDir,sprintf('step_%03d_qualification_small.mat',k)));
        assert(Q.pass,'isolatedtip:Qualification','Fixed P%d core %.6g qualification failed.',k,scale);
        row=Q.summary;row.reference_segment=k;
        D.fixedQualification=[D.fixedQualification;row];
        if opt.AllowPhysicalSolves
            R=solve_incremental_crack_tip(Q.candidate,'FrozenState',R0,'AllowSolve',true, ...
                'FastEDI',opt.FastEDI, ...
                'CheckpointFile',fullfile(fixedDir,sprintf('step_%03d_physical_solved.mat',k)), ...
                'SaveFile',fullfile(fixedDir,sprintf('step_%03d_physical_small.mat',k)));
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
F=main_stage2_embed_scaled_core_full_domain_theta0('FrozenState',R0, ...
    'CoreScale',2,'ExteriorScale',1,'Plot',opt.Plot,'ExteriorVerbose',true, ...
    'SaveCandidate',true,'CandidateFile',fullfile(p1Dir,'P1_candidate.mat'), ...
    'SaveCompact',true,'CompactFile',fullfile(p1Dir,'P1_qualification_small.mat'));
assert(F.pass,'isolatedtip:P1Qualification','Independent isolated P1 failed qualification.');
D.P1Qualification=F.summary;local_save(out,D);
if ~opt.AllowPhysicalSolves
    fprintf('Isolated-tip qualification complete. No physical solves performed.\nOutputs: %s\n',out);
    return
end
R1=main_stage2_theta0_physical_solve('FrozenState',R0,'Candidate',F.candidate, ...
    'AlternativeQualifiedCandidate',true,'AllowSolve',true, ...
    'CheckpointFile',fullfile(p1Dir,'P1_physical_solved.mat'), ...
    'SaveFile',fullfile(p1Dir,'P1_physical_small.mat'));
assert(R1.pass,'isolatedtip:P1Physical','Independent isolated P1 failed physical gates.');
D.P1Physical=R1.summary;local_save(out,D);clear F
% Existing alternative-family mode skips only historical reference-value anchors.
% Qualification, synthetic replay, solver/residual, EDI, COD and MTS gates
% and their numerical tolerances remain in the unchanged production routines.
D.independentPath=run_incremental_crack_path('FrozenState',R0,'MaxSegments',23, ...
    'CoreScale',2,'ExteriorScale',1,'AllowPhysicalSolves',true,'RunSynthetic',true, ...
    'FastEDI',opt.FastEDI,'RegressionGates',false,'StopAtCoreClearance',true, ...
    'PlotEachStep',opt.Plot,'ReuseCandidates',true, ...
    'MeshFamilyLabel','isolated_core_2_exterior_1', ...
    'SeedKI',R1.EDI.KI_unit,'SeedKII',R1.EDI.KII_unit, ...
    'OutputDir',fullfile(out,'independent_core_2_exterior_1','trajectory'));
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
