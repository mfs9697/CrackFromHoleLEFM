function Out=main_step38_solve_tip_checkpoint(O25,O34,O37,varargin)
%MAIN_STEP38_SOLVE_TIP_CHECKPOINT
% Solve the APPROVED Step-37 tip-only T3 mesh exactly ONCE. The output
% saves a compact solved-field checkpoint BEFORE any expensive EDI work.
% Use the separate main_step38_postprocess_tip_checkpoint(checkpoint)
% to extract K_I, K_II and native-face COD after the expensive solve.
%
% Mesh geometry is NEVER regenerated: uses exactly O37.p,O37.t,O37.crack.
% All six original Step-37 gates MUST have passed.
%
% Example:
%   C38=main_step38_solve_tip_checkpoint(O25,O34,O37);
%   % Once saved, postprocess needs only this checkpoint file:
%   O38=main_step38_postprocess_tip_checkpoint(C38.checkpointPath);
%
% The saved file intentionally excludes stiffness K, F, stress recovery,
% and the OLD Step-34 mesh/U; only small baseline numeric results are saved.
%
% Note: MATLAB can still run out of memory DURING the sparse direct solve.
% This checkpoint protects against postprocessing failures, not failures
% that occur before the FEM solution has been computed.

ip=inputParser;
addParameter(ip,'CheckpointFile', ...
    fullfile(pwd,'step38_tip_refined_solved.mat'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Overwrite',false, ...
    @(x)(islogical(x)||isnumeric(x))&&isscalar(x));
parse(ip,varargin{:});
opt=ip.Results;
checkpointPath=char(opt.CheckpointFile);

required(O25,'config','O25');
required(O25,'mouth','O25');
for field={'mesh','crack','hTipNew','U','mat','newKI','newKII', ...
        'rInner','rOuterOverA0','CODTable','nativeR','nativeApparent'}
    required(O34,field{1},'O34');
end
for field={'p','t','crack','passed','gates','meshStats','settings'}
    required(O37,field{1},'O37');
end
if ~logical(O37.passed)||~all(table2array(O37.gates))
    error('step38:UnapprovedMesh', ...
        'The original Step-37 mesh did not pass ALL six gates.');
end
if exist(checkpointPath,'file')==2 && ~logical(opt.Overwrite)
    error('step38:ExistingCheckpoint', ...
        ['Existing checkpoint not overwritten: %s. ', ...
         'Use the postprocessor instead, or explicitly Overwrite=true.'], ...
        checkpointPath);
end

p=O37.p;
t=O37.t;
crack=O37.crack;
if ~isequal(p(1:size(O34.mesh.coord3,1),:), ...
        O34.mesh.coord3) || ...
        norm(crack.Pmid-O34.crack.Pmid,'fro')>1e-12 || ...
        norm(crack.Pmid(1,:)-O25.mouth)>1e-11
    error('step38:MeshGeometryMismatch', ...
        'Refined mesh is not nested within the audited Step-34 geometry.');
end
if any(size(t,2)~=3) || ...
        size(t,1)~=O37.meshStats.new_T3_triangles(1) || ...
        size(p,1)~=O37.meshStats.new_T3_vertices(1)
    error('step38:MeshCounts','Step-37 saved mesh/counts disagree.');
end
oldTip=O37.meshStats.old_tip_h(1);
tipExpected=O37.meshStats.refined_tip_h(1);
if abs(oldTip/O34.hTipNew-1)>1e-9 || ...
        (1-tipExpected/oldTip)<O37.settings.MinTipReduction
    error('step38:TipScaleMismatch', ...
        'Step-37 measured tip reduction is inconsistent with Step 34.');
end

a0=norm(diff(crack.Pmid,1,1));
C=O25.config;
C.a0=a0;
C.solver.verbose=0;
% Consistent minimal-anchoring constraint nodes on SAME original vertex
% IDs. Only these two corner IDs are used by solve_cracked_LEFM.
xOld=O34.mesh.coord3;
[~,LB]=min(sum((xOld-[0,-C.B]).^2,2));
[~,RB]=min(sum((xOld-[C.A,-C.B]).^2,2));
G=struct('p',p,'t',t);
G.edgeSets.corners=struct('left_bottom',LB, ...
    'right_bottom',RB);
G.meta=struct('A',C.A,'B',C.B);

baseline=struct();
baseline.KI=O34.newKI(:).';
baseline.KII=O34.newKII(:).';
baseline.CODTable=O34.CODTable;
baseline.nativeR=O34.nativeR;
baseline.nativeApparent=O34.nativeApparent;
baseline.hTip=O34.hTipNew;
baseline.nT6=size(O34.mesh.coord,1);
baseline.rInner=O34.rInner;
baseline.rOuterOverA0=O34.rOuterOverA0(:).';
baseline.crack=O34.crack.Pmid;
baseline.coreRadius=O37.settings.CoreRadius;
baseline.protectedRadius=O37.settings.ProtectRadius;
baseline.meshStats=O37.meshStats;
baseline.gates=O37.gates;

fprintf('\n============================================================\n');
fprintf('STEP 38 / PHASE 1: ONE TIP-REFINED FEM SOLVE + CHECKPOINT\n');
fprintf('============================================================\n');
fprintf('  Exactly using approved Step-37 T3: %d vertices, %d triangles.\n', ...
    size(p,1),size(t,1));
fprintf('  measured tip size: %.10e -> %.10e m\n',oldTip,tipExpected);
fprintf('  No EDI or COD postprocessing during Phase 1.\n');

% Source solver still returns K, forces and stress recovery. Keep only
% mesh, U, material and RELEASE K before checkpoint/postprocessing.
solved=solve_cracked_LEFM(C,G,'lambda',1.0);
mesh=solved.mesh;
U=solved.U;
mat=solved.mat;
clear solved G

actualTip=local_tip_edge_median(mesh.coord3,mesh.connect3, ...
    crack.Pmid(end,:));
if abs(actualTip/tipExpected-1)>1e-9
    error('step38:TipDidNotMatch', ...
        'Solved mesh tip median %.12g differs from approved %.12g.', ...
        actualTip,tipExpected);
end
if numel(U)~=2*size(mesh.coord,1)||any(~isfinite(U))
    error('step38:BadSolvedU', ...
        'The solved field has incompatible dimensions or nonfinite values.');
end

[folder,~,~]=fileparts(checkpointPath);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
% Protect an existing good checkpoint from partial writes/interruption.
tempPath=[checkpointPath '.incomplete.mat'];
if exist(tempPath,'file')==2
    delete(tempPath);
end
try
    save(tempPath,'mesh','U','mat','crack','baseline', ...
        'actualTip','a0','-v7.3');
    % Move the complete file into place. Avoid deleting an old checkpoint
    % BEFORE a complete replacement is available.
    [ok,msg]=movefile(tempPath,checkpointPath,'f');
    if ~ok
        error('step38:MoveCheckpoint','%s',msg);
    end
catch ME
    warning('step38:CheckpointSaveFailed', ...
        'Solve completed but checkpoint save failed: %s',ME.message);
    rethrow(ME);
end

Out=struct('checkpointPath',checkpointPath, ...
    'nT3',size(mesh.connect3,1), ...
    'nT6',size(mesh.coord,1), ...
    'tipOld',oldTip,'tipNew',actualTip, ...
    'oldN_T6',baseline.nT6,'checkpointSaved',true);
fprintf('  tip-refined T6 nodes=%d, measured tip edge=%.10e m\n', ...
    Out.nT6,actualTip);
fprintf('  SOLVED FIELD SAFELY SAVED: %s\n',checkpointPath);
fprintf('  Phase 1 complete. Phase 2 reads checkpoint only.\n');
end

function h=local_tip_edge_median(P,T,tip)
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
if isempty(e)
    error('step38:TipTopology','No tip-adjacent T3 edge found.');
end
h=median(e);
end

function required(S,name,label)
if ~isstruct(S)||~isfield(S,name)||isempty(S.(name))
    error('step38:MissingInput','Missing %s.%s.',label,name);
end
end
