function Files=plot_step34_mesh(O25,O33,O34mesh,varargin)
%PLOT_STEP34_MESH Visualize the EXACT nested T3 mesh used by Step 34.
%
% Uses ACTUAL stored T3 coordinates/connectivity; never regenerates or
% approximates the refined mesh. The first view is plate/hole/crack
% context; baseline and nested annulus views use identical axis limits.
% The final panel zooms in around the tip and the common EDI inner circle.
%
% Usage (no FEM solve; works with existing Step-34 dry-run result):
%   Files=plot_step34_mesh(O25,O33,O34mesh);
%   Files=plot_step34_mesh(O25,O33,O34mesh, ...
%       'SavePrefix','figures/step34_annulus_mesh');
%
% By default produces in the MATLAB current folder:
%   step34_mesh_comparison.png     high-resolution overview and comparison
%   step34_mesh_comparison.fig     editable MATLAB figure
%   step34_mesh_data.mat          original/refined T3 coordinates and
%                                 topology, crack and EDI radii
%
% Files.fig is the MATLAB figure handle and Files.*Path the saved paths.

p=inputParser;
addParameter(p,'SavePrefix',fullfile(pwd,'step34_mesh'), ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'Resolution',260, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=100&&x<=600);
addParameter(p,'Save',true, ...
    @(x)(islogical(x)||isnumeric(x))&&isscalar(x));
addParameter(p,'Visible',true, ...
    @(x)(islogical(x)||isnumeric(x))&&isscalar(x));
parse(p,varargin{:});
O=p.Results;

if ~isfield(O33,'mesh')||~isfield(O33,'crack') || ...
        ~isfield(O34mesh,'p') || ~isfield(O34mesh,'t') || ...
        ~isfield(O34mesh,'crack')
    error('step34plot:MissingMesh', ...
        'Pass the solved Step-33 output and Step-34 dry-run mesh output.');
end
if ~isfield(O25,'config')
    error('step34plot:MissingConfig','Expected O25.config.');
end
P0=O33.mesh.coord3;
T0=O33.mesh.connect3;
P1=O34mesh.p;
T1=O34mesh.t;
cr=O34mesh.crack;
if size(P0,2)~=2||size(P1,2)~=2 || ...
        size(T0,2)~=3||size(T1,2)~=3
    error('step34plot:BadT3','Both input meshes must be 2D T3 meshes.');
end
if size(P1,1)<size(P0,1)|| ...
        max(abs(P1(1:size(P0,1),:)-P0),[],'all')>1e-12
    error('step34plot:NotNested', ...
        'Step-34 refined mesh does not preserve original T3 vertex IDs.');
end
if ~isfield(O33,'rOuterOverA0')||~isfield(O33,'rInner')
    error('step34plot:MissingEDI','O33 must contain EDI integration radii.');
end
a0=norm(cr.Pmid(end,:)-cr.Pmid(1,:));
ro=sort(O33.rOuterOverA0(:).')*a0;
ri=O33.rInner;
tip=cr.Pmid(end,:);
mouth=cr.Pmid(1,:);
if any(abs(cr.Pmid(:)-O33.crack.Pmid(:))>1e-11)
    error('step34plot:MovedCrack', ...
        'Refined and original crack geometries differ.');
end

C=O25.config;
cx=C.hole.center(1);
cy=C.hole.center(2);
R=C.hole.r;
A=C.A;
B=C.B;
if logical(O.Visible)
    visibility='on';
else
    visibility='off';
end

fig=figure('Name','Step 34 | exact nested crack-tip mesh', ...
    'NumberTitle','off','Color','w', ...
    'Units','pixels','Position',[80,70,1390,1000], ...
    'Visible',visibility);
gridLayout=tiledlayout(fig,2,2, ...
    'TileSpacing','compact','Padding','compact');
title(gridLayout, ...
    'Step 34: actual nested T3 refinement around the crack tip', ...
    'FontWeight','bold','FontSize',14);

% Panel 1: actual physical setup, without covering a full plate with
% microscopic element edges that would be illegible at page scale.
ax1=nexttile(gridLayout,1);
hold(ax1,'on');
plot(ax1,[0 A A 0 0],[-B -B B B -B],'-', ...
    'Color',[.25 .33 .40],'LineWidth',1.1);
tt=linspace(0,2*pi,500);
plot(ax1,cx+R*cos(tt),cy+R*sin(tt),'-', ...
    'Color',[.30 .30 .30],'LineWidth',1.2);
plot(ax1,cr.Pmid(:,1),cr.Pmid(:,2),'-', ...
    'Color',[.75 .23 .13],'LineWidth',2.5);
plot(ax1,mouth(1),mouth(2),'o', ...
    'Color',[.75 .23 .13],'MarkerFaceColor',[.75 .23 .13], ...
    'MarkerSize',5);
plot(ax1,tip(1),tip(2),'o','Color',[.15 .30 .73], ...
    'MarkerFaceColor',[.15 .30 .73],'MarkerSize',5);
circleLine(ax1,tip,max(ro),[.23 .48 .78],1.2,'--');
axis(ax1,'equal');
xlim(ax1,[-.012 A+.012]);
ylim(ax1,[-B-.012 B+.012]);
title(ax1,'Geometry: plate, hole, crack and EDI region');
xlabel(ax1,'x (m)');ylabel(ax1,'y (m)');
grid(ax1,'on');box(ax1,'on');
text(ax1,tip(1)+.012,tip(2)-.024, ...
    'Highlighted: largest EDI domain', ...
    'Color',[.23 .48 .78],'FontSize',9);

% Identical zoom bounds let element sizes be compared visually.
annulusWindow=1.13*(max(ro)+max(0,O34mesh.settings.OuterBuffer));
xLimits=[tip(1)-annulusWindow tip(1)+annulusWindow];
yLimits=[tip(2)-annulusWindow tip(2)+annulusWindow];

ax2=nexttile(gridLayout,2);
drawT3(ax2,P0,T0,[.57 .62 .67],.26);
drawCrackAndEDI(ax2,cr,tip,ri,ro);
setMeshAxes(ax2,xLimits,yLimits);
title(ax2,sprintf('Before: Step 33 (%s T3)', ...
    formatInteger(size(T0,1))));
xlabel(ax2,'x (m)');ylabel(ax2,'y (m)');

ax3=nexttile(gridLayout,3);
drawT3(ax3,P1,T1,[.32 .42 .54],.27);
drawCrackAndEDI(ax3,cr,tip,ri,ro);
setMeshAxes(ax3,xLimits,yLimits);
title(ax3,sprintf('After: Step 34 (%s T3)', ...
    formatInteger(size(T1,1))));
xlabel(ax3,'x (m)');ylabel(ax3,'y (m)');

ax4=nexttile(gridLayout,4);
drawT3(ax4,P1,T1,[.29 .39 .52],.35);
drawCrackAndEDI(ax4,cr,tip,ri,ro);
tipWindow=max(1.8*ri,.00155);
setMeshAxes(ax4, ...
    [tip(1)-tipWindow tip(1)+tipWindow], ...
    [tip(2)-tipWindow tip(2)+tipWindow]);
title(ax4,'Refined tip close-up (inner EDI boundary dashed)');
xlabel(ax4,'x (m)');ylabel(ax4,'y (m)');

prefix=char(O.SavePrefix);
Files=struct('fig',fig,'pngPath','','figPath','','matPath','');
if logical(O.Save)
    [folder,~,~]=fileparts(prefix);
    if ~isempty(folder)&&exist(folder,'dir')~=7
        mkdir(folder);
    end
    Files.pngPath=[prefix '_comparison.png'];
    Files.figPath=[prefix '_comparison.fig'];
    Files.matPath=[prefix '_data.mat'];
    exportgraphics(fig,Files.pngPath, ...
        'Resolution',O.Resolution,'BackgroundColor','white');
    savefig(fig,Files.figPath);

    meshBefore=struct('p',P0,'t',T0,'crack',O33.crack);
    meshAfter=struct('p',P1,'t',T1,'crack',cr);
    EDIRadii=struct('rInner',ri,'rOuter',ro);
    save(Files.matPath,'meshBefore','meshAfter','EDIRadii','-v7.3');
    fprintf('\nSTEP 34 ACTUAL T3 MESH VISUALIZATION SAVED:\n');
    fprintf('  PNG (high resolution): %s\n',Files.pngPath);
    fprintf('  MATLAB editable figure: %s\n',Files.figPath);
    fprintf('  Original/refined T3 mesh data: %s\n',Files.matPath);
end
end

function drawT3(ax,P,T,edgeColor,lineWidth)
patch(ax,'Faces',T,'Vertices',P, ...
    'FaceColor','none','EdgeColor',edgeColor, ...
    'LineWidth',lineWidth,'EdgeAlpha',0.9);
hold(ax,'on');
end

function drawCrackAndEDI(ax,cr,tip,ri,ro)
plot(ax,cr.Pmid(:,1),cr.Pmid(:,2),'-', ...
    'Color',[.82 .25 .12],'LineWidth',2.0);
plot(ax,tip(1),tip(2),'o', ...
    'Color',[.82 .25 .12], ...
    'MarkerFaceColor',[.82 .25 .12],'MarkerSize',4);
circleLine(ax,tip,ri,[.80 .25 .15],1.1,'--');
for k=1:numel(ro)
    if k==numel(ro)
        circleLine(ax,tip,ro(k),[.18 .47 .72],1.2,'-');
    else
        circleLine(ax,tip,ro(k),[.30 .55 .75],0.9,':');
    end
end
end

function circleLine(ax,xy,r,col,wd,style)
tt=linspace(0,2*pi,240);
plot(ax,xy(1)+r*cos(tt),xy(2)+r*sin(tt), ...
    'Color',col,'LineWidth',wd,'LineStyle',style);
end

function setMeshAxes(ax,xlimits,ylimits)
axis(ax,'equal');
xlim(ax,xlimits);
ylim(ax,ylimits);
box(ax,'on');
set(ax,'FontSize',9,'Layer','top');
end

function s=formatInteger(n)
% Compatible with MATLAB versions without locale-dependent formatters.
s=sprintf('%d',n);
end
