function O60=main_step60_visualize_step38_mesh(varargin)
%MAIN_STEP60_VISUALIZE_STEP38_MESH
% Read-only visualization of the REAL saved Step38 asymmetric mesh.
%
% Loads ONLY mesh/crack/compact geometry metadata from the existing
% Step38 checkpoint. It never loads U, never solves FEM, never generates
% or modifies a mesh, and never evaluates an interaction integral.
%
% The figure contains:
%   (1) the complete saved Step38 T3 mesh and all free T3 boundaries;
%   (2) the actual hole/crack neighborhood in global coordinates;
%   (3) a crack-tip view in the local crack frame with the historical
%       fixed EDI inner radius and all three outer radii;
%   (4) the exact element SUPPORT of the 16-point FE-nodal-q gradient for
%       the primary r_outer/a0=0.65 domain, reproduced from the selection
%       logic in SIF_LEFM_interaction_EDI.m, but with NO U or EDI density.
%
% This is a topology/geometry diagnostic only.
%
% Usage after pulling sif-asymmetric-mesh-audit:
%   addpath(genpath(pwd));
%   O60=main_step60_visualize_step38_mesh();
%   disp(O60.summary);
%   disp(O60.supportTable);
%
% The default checkpoint is the investigator-local file created in Step38:
%   verification/step38_tip_refined_solved.mat

root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
vdir=fullfile(root,'verification');

ip=inputParser;
addParameter(ip,'CheckpointFile', ...
    fullfile(vdir,'step38_tip_refined_solved.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'PrimaryOuterRatio',0.65, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<1);
addParameter(ip,'SavePrefix', ...
    fullfile(vdir,'step60_step38_real_mesh'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'Visible','on', ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;

cp=char(opt.CheckpointFile);
prefix=char(opt.SavePrefix);
vis=char(opt.Visible);
if exist(cp,'file')~=2
    error('step60:MissingStep38Checkpoint', ...
        ['Step60 needs the REAL local Step38 checkpoint. Missing: %s\n', ...
         'Expected default: verification/step38_tip_refined_solved.mat'],cp);
end

% IMPORTANT: intentionally do NOT load U or any stiffness/force arrays.
s=load(cp,'mesh','crack','baseline','actualTip','a0');
for f={'mesh','crack','baseline','actualTip','a0'}
    need(s,f{1});
end
mesh=s.mesh;
crack=s.crack;
baseline=s.baseline;
actualTip=s.actualTip;
a0=s.a0;
clear s

for f={'coord3','connect3','coord','connect'}
    need(mesh,f{1});
end
for f={'Pmid','upperNodes','lowerNodes','tipNode'}
    need(crack,f{1});
end
for f={'rInner','rOuterOverA0','hTip'}
    need(baseline,f{1});
end

P3=mesh.coord3;
T3=mesh.connect3(:,1:3);
P6=mesh.coord;
T6=mesh.connect;
if size(P3,2)~=2 || size(T3,2)~=3 || ...
        size(P6,2)~=2 || size(T6,2)~=6 || ...
        size(T3,1)~=size(T6,1)
    error('step60:BadMesh','Expected corresponding T3/T6 Step38 meshes.');
end
if any(~isfinite(P3(:))) || any(~isfinite(P6(:))) || ...
        any(T3(:)<1) || any(T6(:)<1)
    error('step60:InvalidMesh','Stored Step38 mesh contains invalid entries.');
end

% Strong provenance checks for the documented solved Step38 checkpoint.
expectedMouth=[0.19998887220,-0.020817033905];
if abs(a0-0.008)>1e-10 || ...
        size(T3,1)~=39441 || size(P3,1)~=20164 || ...
        size(P6,1)~=79769 || ...
        abs(actualTip-5.4024650785e-5)>1e-12 || ...
        norm(crack.Pmid(1,:)-expectedMouth)>2e-10 || ...
        abs(norm(diff(crack.Pmid,1,1))-a0)>1e-10 || ...
        abs(baseline.rInner-0.0008)>1e-12
    error('step60:CheckpointProvenance', ...
        ['Local file does not match the documented Step38 checkpoint ', ...
         '(8-mm crack, 39,441 T3, 79,769 T6, fixed 0.8-mm inner radius).']);
end

rRat=baseline.rOuterOverA0(:).';
if numel(rRat)~=3 || max(abs(sort(rRat)-[.50 .65 .80]))>1e-12
    error('step60:UnexpectedEDIRadii', ...
        'Expected the documented Step38 outer ratios [0.50 0.65 0.80].');
end
[~,iPrimary]=min(abs(rRat-opt.PrimaryOuterRatio));
if abs(rRat(iPrimary)-opt.PrimaryOuterRatio)>1e-12
    error('step60:PrimaryRadiusAbsent', ...
        'PrimaryOuterRatio must match one of the stored Step38 domains.');
end

tip=crack.Pmid(end,:);
mouth=crack.Pmid(1,:);
e1=(tip-mouth)/a0;
e2=[-e1(2),e1(1)];
R=[e1(:),e2(:)];
X3=(P3-tip)*R;
X6=(P6-tip)*R;

% Free T3 boundaries by NODE ID. This retains the two distinct coincident
% crack faces separately and requires no reconstructed PDE geometry.
freeEdges=free_t3_edges(T3);

% Reproduce the element-selection part of the FE-nodal-q EDI for all
% three historical domains. No displacement field or interaction density
% is touched.
support=cell(numel(rRat),1);
supportRows=nan(numel(rRat),8);
for k=1:numel(rRat)
    ro=rRat(k)*a0;
    ids=fe_nodal_q_support(P6,T6,tip,e1,e2,baseline.rInner,ro);
    support{k}=ids(:);

    nodes=unique(T6(ids,:));
    rr=hypot(X6(nodes,1),X6(nodes,2));
    cen=(X3(T3(ids,1),:)+X3(T3(ids,2),:)+X3(T3(ids,3),:))/3;
    rc=hypot(cen(:,1),cen(:,2));
    dMouth=min(hypot(P6(nodes,1)-mouth(1),P6(nodes,2)-mouth(2)));
    supportRows(k,:)=[rRat(k),1e3*baseline.rInner,1e3*ro, ...
        numel(ids),1e3*min(rr),1e3*max(rr), ...
        1e3*max(rc),1e3*dMouth];
end
supportTable=array2table(supportRows,'VariableNames',{ ...
    'r_outer_over_a0','r_inner_mm','r_outer_mm', ...
    'n_q_support_elements','min_support_node_r_mm', ...
    'max_support_node_r_mm','max_support_centroid_r_mm', ...
    'min_support_node_distance_to_mouth_mm'});

primaryIDs=support{iPrimary};
if isempty(primaryIDs)
    error('step60:NoPrimarySupport','Primary FE-nodal-q support is empty.');
end

% Basic geometry/topology summary.
up=unique(crack.upperNodes(:));
lo=unique(crack.lowerNodes(:));
shared=intersect(up,lo);
faceGridMismatch=paired_face_coordinate_mismatch(P3,up,lo,mouth,e1);
summary=table( ...
    size(P3,1),size(T3,1),size(P6,1), ...
    1e3*a0,1e3*actualTip, ...
    numel(up),numel(lo),numel(shared), ...
    faceGridMismatch*1e3, ...
    size(freeEdges,1), ...
    'VariableNames',{ ...
      'nT3Nodes','nT3Triangles','nT6Nodes', ...
      'a0_mm','tipMedianEdge_mm', ...
      'upperCrackFaceT3Nodes','lowerCrackFaceT3Nodes', ...
      'sharedCrackFaceT3NodeIDs','pairedFaceCoordMismatch_mm', ...
      'nFreeT3Edges'});

fprintf('\n============================================================\n');
fprintf('STEP 60: REAL STEP38 ASYMMETRIC MESH — VISUALIZATION ONLY\n');
fprintf('============================================================\n');
fprintf('  Checkpoint: %s\n',cp);
fprintf('  Stored mesh: T3 nodes=%d, T3 triangles=%d, T6 nodes=%d.\n', ...
    size(P3,1),size(T3,1),size(P6,1));
fprintf('  Crack mouth=[%.12g %.12g] m; a0=%.6g m.\n', ...
    mouth(1),mouth(2),a0);
fprintf('  Tip=[%.12g %.12g] m; measured hTip=%.12g m.\n', ...
    tip(1),tip(2),actualTip);
fprintf('  NO U loaded; NO mesh generation, solve, stiffness, COD or EDI.\n');
fprintf('\nFE-NODAL-q ELEMENT SUPPORT (selection logic only)\n');
disp(supportTable);

% ---------- Figure 1: four-panel overview of the REAL mesh ----------
fig=figure('Name','Step60: real Step38 asymmetric mesh', ...
    'Color','w','NumberTitle','off','Visible',vis, ...
    'Position',[70 55 1450 900]);
tl=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');

% Panel 1: entire stored mesh.
ax1=nexttile(tl);
draw_global_mesh(ax1,P3,T3,freeEdges,crack);
title(ax1,'A. Complete saved Step38 T3 mesh');
xlabel(ax1,'x (m)');ylabel(ax1,'y (m)');

% Panel 2: hole + crack neighborhood using actual stored topology.
ax2=nexttile(tl);
draw_global_mesh(ax2,P3,T3,freeEdges,crack);
xlim(ax2,[mouth(1)-0.040, tip(1)+0.016]);
ylim(ax2,[mouth(2)-0.040, mouth(2)+0.040]);
title(ax2,'B. Actual retained-hole / crack neighborhood');
xlabel(ax2,'x (m)');ylabel(ax2,'y (m)');

% Panel 3: crack-local mesh and historical EDI circles.
ax3=nexttile(tl);
rMax=max(rRat)*a0;
L=1.12*rMax;
draw_local_mesh(ax3,X3,T3,L);
hold(ax3,'on');
draw_circle_mm(ax3,baseline.rInner,':','r_{in}=0.8 mm');
styles={'--','-','-.'};
for k=1:numel(rRat)
    draw_circle_mm(ax3,rRat(k)*a0,styles{k}, ...
        sprintf('r_o=%.1f mm',1e3*rRat(k)*a0));
end
plot_crack_faces_local(ax3,X3,up,lo,crack.tipNode);
plot(ax3,0,0,'ko','MarkerFaceColor','k','MarkerSize',5, ...
    'DisplayName','tip');
plot(ax3,-1e3*a0,0,'ks','MarkerFaceColor','w','MarkerSize',6, ...
    'DisplayName','mouth');
axis(ax3,'equal');
xlim(ax3,1e3*[-L L]);ylim(ax3,1e3*[-L L]);
box(ax3,'on');grid(ax3,'on');
xlabel(ax3,'x_1 from tip (mm)');ylabel(ax3,'x_2 (mm)');
title(ax3,'C. Real tip topology and fixed historical EDI annuli');
legend(ax3,'Location','bestoutside');

% Panel 4: exact FE-nodal-q element support for primary 0.65a0 domain.
ax4=nexttile(tl);
draw_local_mesh(ax4,X3,T3,L);
hold(ax4,'on');
patch(ax4,'Faces',T3(primaryIDs,:), ...
    'Vertices',1e3*X3,'FaceColor',[0.88 0.78 0.45], ...
    'FaceAlpha',0.42,'EdgeColor','none', ...
    'DisplayName','nonzero FE-nodal-q gradient support');
draw_circle_mm(ax4,baseline.rInner,':','r_{in}');
draw_circle_mm(ax4,rRat(iPrimary)*a0,'-','r_o=0.65a_0');
plot(ax4,0,0,'ko','MarkerFaceColor','k','MarkerSize',5, ...
    'DisplayName','tip');
plot(ax4,-1e3*a0,0,'ks','MarkerFaceColor','w','MarkerSize',6, ...
    'DisplayName','mouth');
axis(ax4,'equal');
xlim(ax4,1e3*[-L L]);ylim(ax4,1e3*[-L L]);
box(ax4,'on');grid(ax4,'on');
xlabel(ax4,'x_1 from tip (mm)');ylabel(ax4,'x_2 (mm)');
title(ax4,sprintf(['D. Actual 16-GP FE-nodal-q support, ', ...
    'r_o/a_0=%.2f (%d T6 elements)'], ...
    rRat(iPrimary),numel(primaryIDs)));
legend(ax4,'Location','bestoutside');

title(tl,{ ...
    'Step 60 — actual saved Step38 asymmetric mesh topology', ...
    'Read-only checkpoint visualization; no reconstruction and no FEM/EDI solve'});

[folder,~,~]=fileparts(prefix);
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
files=struct();
files.overviewPNG=[prefix '_overview.png'];
exportgraphics(fig,files.overviewPNG, ...
    'Resolution',240,'BackgroundColor','white');

% ---------- Figure 2: larger primary support close-up ----------
fig2=figure('Name','Step60: Step38 primary q support', ...
    'Color','w','NumberTitle','off','Visible',vis, ...
    'Position',[110 75 1050 850]);
ax=axes(fig2);
draw_local_mesh(ax,X3,T3,L);
hold(ax,'on');
patch(ax,'Faces',T3(primaryIDs,:), ...
    'Vertices',1e3*X3,'FaceColor',[0.88 0.78 0.45], ...
    'FaceAlpha',0.46,'EdgeColor','none', ...
    'DisplayName','FE-nodal-q gradient support');
draw_circle_mm(ax,baseline.rInner,':','r_{in}=0.8 mm');
draw_circle_mm(ax,rRat(iPrimary)*a0,'-', ...
    sprintf('nominal r_o=%.1f mm',1e3*rRat(iPrimary)*a0));
plot_crack_faces_local(ax,X3,up,lo,crack.tipNode);
plot(ax,0,0,'ko','MarkerFaceColor','k','MarkerSize',6, ...
    'DisplayName','crack tip');
plot(ax,-1e3*a0,0,'ks','MarkerFaceColor','w','MarkerSize',7, ...
    'DisplayName','crack mouth');
axis(ax,'equal');xlim(ax,1e3*[-L L]);ylim(ax,1e3*[-L L]);
box(ax,'on');grid(ax,'on');
xlabel(ax,'local x_1 from crack tip (mm)');
ylabel(ax,'local x_2 (mm)');
title(ax,{ ...
    'Real Step38 mesh: primary FE-nodal-q support', ...
    sprintf(['Nominal annulus 0.8–%.1f mm; support includes ', ...
             'straddling T6 elements'],1e3*rRat(iPrimary)*a0)});
legend(ax,'Location','bestoutside');
files.primarySupportPNG=[prefix '_primary_q_support.png'];
exportgraphics(fig2,files.primarySupportPNG, ...
    'Resolution',260,'BackgroundColor','white');

O60=struct( ...
    'checkpointPath',cp, ...
    'summary',summary, ...
    'supportTable',supportTable, ...
    'primaryOuterRatio',rRat(iPrimary), ...
    'primarySupportElementIDs',primaryIDs, ...
    'allSupportElementIDs',{support}, ...
    'rInner',baseline.rInner, ...
    'rOuterOverA0',rRat, ...
    'tip',tip,'mouth',mouth,'a0',a0, ...
    'files',files, ...
    'noFEM',true,'noMeshGeneration',true,'noEDI',true);

smallFile=[prefix '_small_data.mat'];
O60.files.smallData=smallFile;
save(smallFile,'O60');

fprintf('\nSTEP 60 OUTPUT FILES\n');
fprintf('  Overview PNG: %s\n',files.overviewPNG);
fprintf('  Primary support PNG: %s\n',files.primarySupportPNG);
fprintf('  Compact topology report: %s\n',smallFile);
fprintf('STEP 60 complete: real Step38 checkpoint visualized; zero solves.\n');
end

% ========================================================================
function freeEdges=free_t3_edges(T)
E=[T(:,[1 2]);T(:,[2 3]);T(:,[3 1])];
Es=sort(E,2);
[Eu,~,ic]=unique(Es,'rows');
cnt=accumarray(ic,1);
freeEdges=Eu(cnt==1,:);
end

function mismatch=paired_face_coordinate_mismatch(P,up,lo,mouth,e1)
% Pair collapsed crack-face T3 coordinates by physical crack coordinate,
% never by node number. Upper/lower IDs are intentionally distinct.
su=(P(up,:)-mouth)*e1(:);
sl=(P(lo,:)-mouth)*e1(:);
[su,iu]=sort(su);[sl,il]=sort(sl);
if numel(su)~=numel(sl)
    mismatch=Inf;
    return
end
mismatch=max(hypot(P(up(iu),1)-P(lo(il),1), ...
                   P(up(iu),2)-P(lo(il),2)));
end

function ids=fe_nodal_q_support(coord,connect,tip,e1,e2,rInner,rOuter)
% Reproduce ONLY the element-participation logic of
% SIF_LEFM_interaction_EDI for WeightFunction='fe_nodal',
% QuadratureRule=16. No displacement U, auxiliary fields, density or
% integration is evaluated.
Rgl=[e1(:),e2(:)];
Rloc=Rgl.';
xTip=tip(:);
Xloc=(Rloc*(coord-xTip.').').';
r=hypot(Xloc(:,1),Xloc(:,2));
qNode=ones(size(coord,1),1);
qNode(r>=rOuter)=0;
mid=(r>rInner)&(r<rOuter);
qNode(mid)=(rOuter-r(mid))/(rOuter-rInner);

[~,xip]=rule16_points();
used=false(size(connect,1),1);
for e=1:size(connect,1)
    nodes=connect(e,:);
    qel=qNode(nodes);
    % Only an EXACTLY constant nodal q can be skipped without changing
    % the production EDI participation test. Do not use an approximate
    % range tolerance here: tiny nodal differences can be amplified by
    % spatial differentiation on small elements.
    if all(qel==qel(1))
        continue
    end
    X=coord(nodes,:);
    Xc=mean(X(1:3,:),1).';
    rc=norm(Rloc*(Xc-xTip));
    if rc>rOuter+element_radius(X)
        continue
    end
    for igp=1:size(xip,2)
        [DetJ,dNdx]=t6_shape_grad(xip(:,igp),X);
        if DetJ<=0
            error('step60:BadT6Element', ...
                'Inverted or degenerate stored T6 element %d.',e);
        end
        qgradGl=dNdx*qel;
        qgrad=Rloc*qgradGl;
        if norm(qgrad)>1e-14
            used(e)=true;
            break
        end
    end
end
ids=find(used);
end

function rad=element_radius(X)
xc=mean(X(1:3,:),1);
rad=max(sqrt(sum((X-xc).^2,2)));
end

function [nip,xip]=rule16_points()
% Same 16-point Dunavant locations as SIF_LEFM_interaction_EDI.m.
xip=zeros(2,16);
xip(:,1)=[1/3;1/3];
a=0.170569307751760;b=0.658861384496480;
xip(:,2:4)=[a a b;a b a];
a=0.050547228317031;b=0.898905543365938;
xip(:,5:7)=[a a b;a b a];
a=0.459292588292723;b=0.081414823414554;
xip(:,8:10)=[a a b;a b a];
a=0.263112829634638;
b=0.728492392955404;
d=0.008394777409958;
xip(:,11:16)=[a a b b d d;b d a d a b];
nip=16;
end

function [DetJ,dNdx]=t6_shape_grad(xi,X)
L1=xi(1);L2=xi(2);L3=1-L1-L2;
dN_dL1=[ ...
    4*L1-1;0;-(4*L3-1);4*L2;-4*L2;4*(L3-L1)];
dN_dL2=[ ...
    0;4*L2-1;-(4*L3-1);4*L1;4*(L3-L2);-4*L1];
dNdxi=[dN_dL1.';dN_dL2.'];
J=dNdxi*X;
DetJ=det(J);
dNdx=J\dNdxi;
end

function draw_global_mesh(ax,P,T,freeEdges,crack)
patch(ax,'Faces',T,'Vertices',P,'FaceColor','none', ...
    'EdgeColor',[.72 .75 .78],'LineWidth',.12);
hold(ax,'on');
plot_edges(ax,P,freeEdges,[.10 .10 .10],.65);
plot(ax,crack.Pmid(:,1),crack.Pmid(:,2), ...
    'LineWidth',1.8,'DisplayName','crack');
plot(ax,crack.Pmid(end,1),crack.Pmid(end,2),'ko', ...
    'MarkerFaceColor','k','MarkerSize',4,'DisplayName','tip');
axis(ax,'equal');box(ax,'on');
end

function draw_local_mesh(ax,X,T,L)
cen=(X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3;
near=abs(cen(:,1))<=1.20*L & abs(cen(:,2))<=1.20*L;
patch(ax,'Faces',T(near,:),'Vertices',1e3*X, ...
    'FaceColor','none','EdgeColor',[.60 .64 .68],'LineWidth',.22);
hold(ax,'on');
end

function plot_crack_faces_local(ax,X,up,lo,tipID)
% The two face coordinate sets are coincident after collapse. Different
% marker types make their distinct topology visible without moving them.
u=setdiff(up,tipID,'stable');
l=setdiff(lo,tipID,'stable');
plot(ax,1e3*X(u,1),1e3*X(u,2),'o','MarkerSize',3.5, ...
    'LineStyle','none','DisplayName','upper crack-face IDs');
plot(ax,1e3*X(l,1),1e3*X(l,2),'+','MarkerSize',4, ...
    'LineStyle','none','DisplayName','lower crack-face IDs');
end

function draw_circle_mm(ax,r,style,label)
th=linspace(0,2*pi,360);
plot(ax,1e3*r*cos(th),1e3*r*sin(th), ...
    'LineStyle',style,'LineWidth',1.1,'DisplayName',label);
end

function plot_edges(ax,P,E,col,lw)
if isempty(E),return,end
x=[P(E(:,1),1),P(E(:,2),1),nan(size(E,1),1)].';
y=[P(E(:,1),2),P(E(:,2),2),nan(size(E,1),1)].';
plot(ax,x(:),y(:),'Color',col,'LineWidth',lw, ...
    'HandleVisibility','off');
end

function need(s,f)
if ~isstruct(s)||~isfield(s,f)||isempty(s.(f))
    error('step60:MissingField','Required field %s is missing.',f);
end
end
