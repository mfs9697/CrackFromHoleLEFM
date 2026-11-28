% IAM-UA8-F02: deterministic ILLUSTRATIVE regeneration of accepted M0 meshes with the pinned historical generator.
% Generator: clean checkout of cfdf0f110010d688a4e3c48f6d88a00fd17dc698 (sources equal source_manifest.csv pins, CRLF form).
% Inputs (read-only): paper/data/evidence_exact.mat (accepted vertices) and paper/data/accepted_stage1_source.mat (frozen R0).
% qualify_incremental_crack_candidate with RunSynthetic=false (no synthetic EDI), SaveCandidate=false, SaveCompact=false:
% mesh construction and structural gates ONLY. No stiffness assembly, no FEM solve, no SIF/EDI evaluation.
gen='<generator_cfdf0f1>';
wt='<repo>';
out='<scratch>\f02';
restoredefaultpath;addpath(genpath(gen));
fprintf('GEN qualify=%s\n',which('qualify_incremental_crack_candidate'));
[~,h]=system(['git -C "' gen '" rev-parse HEAD']);fprintf('GEN_HEAD %s\n',strtrim(h));
[~,s]=system(['git -C "' gen '" status --porcelain']);fprintf('GEN_STATUS [%s]\n',strtrim(s));
Z=load(fullfile(wt,'paper','data','accepted_stage1_source.mat'),'R0');R0=Z.R0;
Ev=load(fullfile(wt,'paper','data','evidence_exact.mat'));fprintf('EVIDENCE fields: %s\n',strjoin(fieldnames(Ev),','));
E=Ev.E;fprintf('E fields: %s\n',strjoin(fieldnames(E),','));
V=E.vertices_m;fprintf('vertices %dx%d\n',size(V,1),size(V,2));
% preserved early checkpoint (profiling run of the same production path, 2026-10-05) for a determinism cross-check
pres='<investigator_codex_work>\2026-10-04\files-pasted-by-the-user-you\work\incremental-profile\phaseA_baseline\fresh5\step_005_candidate.mat';
P=load(pres);fprintf('PRESERVED vars: %s\n',strjoin(fieldnames(P),','));
c=P.candidate;fprintf('PRESERVED candidate fields: %s\n',strjoin(fieldnames(c),','));
if isfield(c,'path'),fprintf('PRESERVED path rows %d; max |path - evidence P0..P5| = %.3e m\n',size(c.path,1),max(abs(c.path-V(1:size(c.path,1),:)),[],'all'));end
R=struct();
for k=[5 21]
 t0=tic;
 F=qualify_incremental_crack_candidate(V(1:k+1,:),'FrozenState',R0,'RunSynthetic',false,'SaveCandidate',false,'SaveCompact',false,'Plot',false,'ExteriorVerbose',false);
 fprintf('REGEN P%d structuralPass=%d time=%.1fs\n',k,all(structfun(@logical,F.gates)),toc(t0));
 fprintf('  mesh fields: %s\n',strjoin(fieldnames(F.mesh),','));fprintf('  localMesh fields: %s\n',strjoin(fieldnames(F.localMesh),','));
 fprintf('  crackLocal fields: %s\n',strjoin(fieldnames(F.crackLocal),','));
 disp(F.summary);
 save(fullfile(out,sprintf('F02_regen_P%02d.mat',k)),'F','-v7.3');
 R.(sprintf('P%d',k))=struct('gates',F.gates,'summary',table2struct(F.summary));
end
save(fullfile(out,'F02_regen_summary.mat'),'R','-v7');
fprintf('REGEN_DONE\n');
