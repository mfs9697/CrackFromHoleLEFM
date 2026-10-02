function O37=main_step37_tip_core_refinement_dryrun(O25,O34,varargin)
%MAIN_STEP37_TIP_CORE_REFINEMENT_DRYRUN
% Geometry-only first test for separately refining the immediate crack-tip
% region on the ALREADY SOLVED Step-34 T3 mesh.
%
% The broad EDI annuli were genuinely refined in Steps 33/34 and EDI
% became stable. Step 36 showed cutoff-dependent native COD and scatter
% within a few tip-edge lengths. This test changes the tip resolution
% while preserving the outer EDI-region TRIANGLE CONNECTIVITY exactly
% beyond a protected radius, along with all original physical vertices,
% distinct crack-face IDs, and the physical crack/hole geometry.
%
% Reuses the tested conforming T3 edge-bisection helper on a small disk.
% The default is ONE pass. In contrast to Step 34, a decrease of the
% measured tip-edge size is now the MAIN acceptance gate.
%
% NO FEM solve, EDI integration, extrapolation, or large mesh copies from
% the Step-32 background fields. If the dry run passes, inspect its
% actual mesh and approximate resulting problem size before authorizing
% a separately developed memory-conscious tip-refinement FEM comparison.
%
% Example:
%   O37=main_step37_tip_core_refinement_dryrun(O25,O34);
%   O37.gates
%   O37.meshStats

ip=inputParser;
addParameter(ip,'CoreRadius',0.0012, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'TargetEdge',8.0e-5, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'MaxPasses',1, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x<=3&&x==round(x));
addParameter(ip,'ProtectRadius',0.0022, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'MinTipReduction',0.20, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<1);
addParameter(ip,'MaxTriangleGrowthFactor',1.30, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>1);
addParameter(ip,'Plot',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'SavePrefix',fullfile(pwd,'step37_tip_core_mesh'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});
opt=ip.Results;

for key={'config','mouth'}
    need(O25,key{1},'O25');
end
for key={'mesh','crack','hTipNew','rInner','rOuterOverA0'}
    need(O34,key{1},'O34');
end
P0=O34.mesh.coord3;
T0=O34.mesh.connect3;
cr0=O34.crack;
tip=cr0.Pmid(end,:);
a0=norm(diff(cr0.Pmid,1,1));
rMax=max(O34.rOuterOverA0(:))*a0;
if abs(opt.ProtectRadius-opt.CoreRadius)<1e-12 || ...
        opt.ProtectRadius<=opt.CoreRadius || ...
        opt.ProtectRadius>=min(rMax,0.6*a0)
    error('step37:BadProtection', ...
        'Require CoreRadius < ProtectRadius < EDI largest outer radius.');
end
if norm(cr0.Pmid(1,:)-O25.mouth)>1e-11
    error('step37:DifferentMouth','O34 and O25 mouth points disagree.');
end
if opt.TargetEdge >= O34.hTipNew
    warning('step37:TargetNearOldTip', ...
        'Target edge is not appreciably smaller than existing htip.');
end

fprintf('\n============================================================\n');
fprintf('STEP 37: IMMEDIATE-TIP SPATIAL REFINEMENT, DRY RUN ONLY\n');
fprintf('============================================================\n');
fprintf('  fixed crack mouth=[%.10e %.10e] m, a0=%.7g m\n', ...
    cr0.Pmid(1,1),cr0.Pmid(1,2),a0);
fprintf('  original T3 triangles=%d; old T6 nodes=%d\n', ...
    size(T0,1),size(O34.mesh.coord,1));
fprintf('  tip disk radius=%.6g, target edge=%.6g m; ', ...
    opt.CoreRadius,opt.TargetEdge);
fprintf('protected connectivity from radius %.6g m outward\n', ...
    opt.ProtectRadius);

% Existing helper requires 0 < rInner < rOuter. Use a tiny positive
% number to include the tip disk, with zero buffer. This is a separate
% experiment and does not modify the shared annular helper.
rTiny=min(1e-10,1e-7*opt.CoreRadius);
[P1,T1,cr1,refInfo]=refine_collapsed_t3_annulus( ...
    P0,T0,cr0,rTiny,opt.CoreRadius,opt.TargetEdge, ...
    'OuterBuffer',0,'InnerBuffer',0, ...
    'MaxPasses',opt.MaxPasses,'Verbose',true);

if ~isequal(P1(1:size(P0,1),:),P0)
    error('step37:OriginalNodesMoved', ...
        'The nested mesh changed original vertex coordinates.');
end
if norm(cr1.Pmid-cr0.Pmid,'fro')>1e-12 || ...
        cr1.tipNode~=cr0.tipNode
    error('step37:CrackPathChanged', ...
        'The tip or physical crack geometry changed.');
end
if ~isempty(setdiff(intersect(cr1.upperNodes,cr1.lowerNodes), ...
        cr1.tipNode))
    error('step37:CrackFacesMerged', ...
        'Refined crack-face node sets intersect away from the tip.');
end

% Demand EXACT same protected outer-shell T3 triangles AND old node IDs.
% This is a stricter test than comparing global node count or h-median.
oldProtected=protected_triangles(P0,T0,tip,opt.ProtectRadius,rMax);
newProtected=protected_triangles(P1,T1,tip,opt.ProtectRadius,rMax);
shellIdentical=isequal(oldProtected,newProtected);
newVertices=P1(size(P0,1)+1:end,:);
rNew=hypot(newVertices(:,1)-tip(1),newVertices(:,2)-tip(2));
if isempty(rNew)
    maxNewRadius=NaN;
else
    maxNewRadius=max(rNew);
end
insideProtected=isempty(rNew)||all(rNew<opt.ProtectRadius);
remoteBoundaryClear=isempty(newVertices)|| ...
    ~any(abs(newVertices(:,1))<1e-10 | ...
         abs(newVertices(:,1)-O25.config.A)<1e-10 | ...
         abs(newVertices(:,2)-O25.config.B)<1e-10 | ...
         abs(newVertices(:,2)+O25.config.B)<1e-10);

tipOld=tip_edge_median(P0,T0,tip);
tipNew=tip_edge_median(P1,T1,tip);
if abs(tipOld/O34.hTipNew-1)>1e-9
    error('step37:TipBaselineMismatch', ...
        'Stored old hTip does not match its own source T3 mesh.');
end
tipReduction=1-tipNew/tipOld;
growth=size(T1,1)/size(T0,1);
sizeGate=tipReduction>=opt.MinTipReduction;
memoryGate=growth<=opt.MaxTriangleGrowthFactor;
faceIncrease=numel(cr1.upperNodes)>numel(cr0.upperNodes) && ...
    numel(cr1.lowerNodes)>numel(cr0.lowerNodes);

R=[tipOld,tipNew,tipReduction, ...
    size(P0,1),size(P1,1), ...
    size(T0,1),size(T1,1),growth, ...
    numel(cr0.upperNodes),numel(cr1.upperNodes), ...
    numel(cr0.lowerNodes),numel(cr1.lowerNodes), ...
    maxNewRadius];
vars={'old_tip_h','refined_tip_h','tip_size_reduction', ...
    'old_T3_vertices','new_T3_vertices', ...
    'old_T3_triangles','new_T3_triangles','T3_growth_factor', ...
    'old_upper_face_T3','new_upper_face_T3', ...
    'old_lower_face_T3','new_lower_face_T3','max_new_vertex_radius'};
meshStats=array2table(R,'VariableNames',vars);

gates=table(shellIdentical,insideProtected,remoteBoundaryClear, ...
    sizeGate,memoryGate,faceIncrease, ...
    'VariableNames',{'protected_outer_T3_identical', ...
    'all_new_vertices_inside_protected_radius', ...
    'outer_plate_boundary_unchanged','tip_reduction_pass', ...
    'triangle_growth_within_limit','both_crack_faces_refined'});
passed=all(table2array(gates));
fprintf('\nDIRECT TIP-REFINEMENT GEOMETRY DIAGNOSTICS\n');
disp(meshStats);
fprintf('\nGATES (all must be true before any further FEM solve)\n');
disp(gates);
fprintf('  overall dry-run acceptance: %d\n',passed);
if ~passed
    fprintf(['  Stop before FEM; inspect which gate failed. ', ...
        'Do not choose mesh parameters solely to influence KII.\n']);
end

files=struct('png','','fig','');
fig=[];
if opt.Plot
    fig=figure('Name','Step37: existing vs refined physical crack-tip T3', ...
        'Color','w','NumberTitle','off','Position',[95 90 1200 570]);
    tl=tiledlayout(fig,1,2, ...
        'Padding','compact','TileSpacing','compact');
    L=max(opt.ProtectRadius*1.12,opt.CoreRadius*1.5);
    ax=nexttile(tl);
    draw_mesh(ax,P0,T0,tip,cr0,O34.rInner, ...
        opt.CoreRadius,opt.ProtectRadius,L);
    title(ax,'Existing Step 34 tip');
    ax=nexttile(tl);
    draw_mesh(ax,P1,T1,tip,cr1,O34.rInner, ...
        opt.CoreRadius,opt.ProtectRadius,L);
    title(ax,'Step 37 new tip-only nested mesh');
    title(tl,sprintf(['Exact T3 comparison; old h_{tip}=%.3f mm, ', ...
        'new h_{tip}=%.3f mm'],1e3*tipOld,1e3*tipNew));
    prefix=char(opt.SavePrefix);
    [folder,~,~]=fileparts(prefix);
    if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
    files.png=[prefix '.png'];
    files.fig=[prefix '.fig'];
    exportgraphics(fig,files.png, ...
        'Resolution',220,'BackgroundColor','white');
    savefig(fig,files.fig);
    fprintf('  Mesh PNG: %s\n',files.png);
end

% Save only the new T3 mesh, NOT any FEM system, stiffness matrix, old
% full mesh duplicate or computationally heavy EDI quadrature diagnostics.
O37=struct('settings',opt,'tip',tip, ...
    'p',P1,'t',T1,'crack',cr1,'refinement',refInfo, ...
    'meshStats',meshStats,'gates',gates,'passed',passed, ...
    'protectedOldCount',size(oldProtected,1), ...
    'protectedNewCount',size(newProtected,1), ...
    'figure',fig,'files',files);
fprintf('STEP 37 dry run completed; zero FEM solves.\n');
end

function A=protected_triangles(P,T,tip,rMin,rMax)
T=T(:,1:3);
cen=(P(T(:,1),:)+P(T(:,2),:)+P(T(:,3),:))/3;
r=hypot(cen(:,1)-tip(1),cen(:,2)-tip(2));
A=sortrows(sort(T(r>=rMin&r<=rMax,:),2));
end

function h=tip_edge_median(P,T,tip)
r=hypot(P(:,1)-tip(1),P(:,2)-tip(2));
tol=max(1e-12,1e-8*max(1,max(abs(P(:)))));
ids=find(r<=min(r)+tol);
tri=T(any(ismember(T,ids),2),:);
if isempty(tri)
    error('step37:NoTipTriangles','Tip-adjacent triangles missing.');
end
p1=P(tri(:,1),:);p2=P(tri(:,2),:);p3=P(tri(:,3),:);
L=[hypot(p1(:,1)-p2(:,1),p1(:,2)-p2(:,2)); ...
    hypot(p2(:,1)-p3(:,1),p2(:,2)-p3(:,2)); ...
    hypot(p3(:,1)-p1(:,1),p3(:,2)-p1(:,2))];
L=L(isfinite(L)&L>tol);
h=median(L);
end

function draw_mesh(ax,P,T,tip,cr,rInner,rCore,rProtect,L)
% Cull elements outside the displayed neighborhood to keep interactive
% rendering affordable on large background meshes.
cen=(P(T(:,1),:)+P(T(:,2),:)+P(T(:,3),:))/3;
near=abs(cen(:,1)-tip(1))<=1.3*L & ...
    abs(cen(:,2)-tip(2))<=1.3*L;
patch(ax,'Faces',T(near,:),'Vertices',P,'FaceColor','none', ...
    'EdgeColor',[.38 .48 .57],'LineWidth',.32);
hold(ax,'on');
plot(ax,cr.Pmid(:,1),cr.Pmid(:,2),'r-','LineWidth',2);
draw_circle(ax,tip,rInner,[.65 .20 .25],':');
draw_circle(ax,tip,rCore,[.12 .50 .24],'-');
draw_circle(ax,tip,rProtect,[.15 .35 .65],'--');
axis(ax,'equal');xlim(ax,tip(1)+[-L L]);
ylim(ax,tip(2)+[-L L]);box(ax,'on');
xlabel(ax,'x (m)');ylabel(ax,'y (m)');
end

function draw_circle(ax,tip,r,col,style)
ang=linspace(0,2*pi,200);
plot(ax,tip(1)+r*cos(ang),tip(2)+r*sin(ang), ...
    'Color',col,'LineStyle',style,'LineWidth',1.1);
end

function need(S,k,label)
if ~isstruct(S)||~isfield(S,k)||isempty(S.(k))
    error('step37:MissingInput','Required %s.%s missing.',label,k);
end
end
