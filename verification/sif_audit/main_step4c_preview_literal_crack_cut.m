function Out=main_step4c_preview_literal_crack_cut(varargin)
%MAIN_STEP4C_PREVIEW_LITERAL_CRACK_CUT Geometry gate for the approved lattice.
%   Out=main_step4c_preview_literal_crack_cut() opens three geometry figures.
%   Out=main_step4c_preview_literal_crack_cut('Visible','off', ...
%       'OutputDir','path/to/results') also writes PNG, MAT and JSON results.
%   No displacement field, loading, FEM solve or SIF extraction is performed.

ip=inputParser;
addParameter(ip,'Visible','on',@(x)any(strcmpi(string(x),["on","off"])));
addParameter(ip,'OutputDir','',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'OverlayParent',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:}); opts=ip.Results;

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));
[mesh,info,parent]=build_literal_ring_crack_cut_mesh();
assert(info.auditT3.passed&&info.auditT6.passed, ...
    'literalCutPreview:GeometryGate','The T3 and T6 geometry gates must pass.');

splitChildren=ismember(info.cut.parentElementID(:),info.cut.splitParentIDs);
nearCut=ismember(info.cut.parentElementID(:),info.cut.cutNeighborhoodParentIDs);
[area,quality,minAngle]=triangle_metrics(mesh.connect3,mesh.coord3);
summary=struct('geometryGatePassed',true, ...
    'parentNtheta',info.parent.Ntheta,'parentNr',info.parent.Nr, ...
    'parentR0',info.parent.r0,'parentR1',info.parent.r1, ...
    'parentQ',info.parent.qActual, ...
    'parentVertices',size(parent.coord3,1), ...
    'parentTriangles',size(parent.connect3,1), ...
    'splitParentTriangles',numel(info.cut.splitParentIDs), ...
    'cutNeighborhoodParentTriangles',numel(info.cut.cutNeighborhoodParentIDs), ...
    'unchangedTrianglesOutsideNeighborhood', ...
        size(parent.connect3,1)-numel(info.cut.cutNeighborhoodParentIDs), ...
    'insertedEdgeVertices',numel(info.cut.insertedEdgeNodeIDs), ...
    't3Vertices',size(mesh.coord3,1),'t3Triangles',size(mesh.connect3,1), ...
    't6Nodes',size(mesh.coord,1),'t6Triangles',size(mesh.connect,1), ...
    't3CrackNodePairs',numel(info.cut.crackUpperIDs), ...
    't6CrackNodePairs',numel(mesh.crackUpperT6IDs), ...
    'minimumSignedArea',min(area), ...
    'minimumNearCutQuality',min(quality(nearCut)), ...
    'minimumNearCutAngleDegrees',min(minAngle(nearCut)), ...
    'auditT3',info.auditT3,'auditT6',info.auditT6);

fprintf('\nSTEP 4C: LITERAL LATTICE CRACK-CUT GEOMETRY GATE\n');
fprintf('Parent N=%d, Nr=%d, r0=%.6g, r1=%.6g, q=%.8f\n', ...
    summary.parentNtheta,summary.parentNr,summary.parentR0, ...
    summary.parentR1,summary.parentQ);
fprintf('Parent triangles split: %d / %d\n', ...
    summary.splitParentTriangles,summary.parentTriangles);
fprintf('Unchanged triangles outside cut neighborhood: %d\n', ...
    summary.unchangedTrianglesOutsideNeighborhood);
fprintf('T3 vertices/triangles: %d / %d; T6 nodes/triangles: %d / %d\n', ...
    summary.t3Vertices,summary.t3Triangles,summary.t6Nodes,summary.t6Triangles);
fprintf('Distinct coincident crack-node pairs: T3=%d, T6=%d\n', ...
    summary.t3CrackNodePairs,summary.t6CrackNodePairs);
fprintf('Minimum signed T3 area: %.8g\n',summary.minimumSignedArea);
fprintf('Near-cut quality min: %.6f; angle min: %.6f deg\n', ...
    summary.minimumNearCutQuality,summary.minimumNearCutAngleDegrees);
fprintf('T3 geometry gate: PASS; T6 geometry gate: PASS\n');

figures=gobjects(3,1);
figures(1)=new_figure('Step 4C: full literal crack-cut mesh',[1100,1000],opts);
ax=axes(figures(1));
draw_mesh(ax,mesh,parent,splitChildren,opts.OverlayParent,false);
axis(ax,[-1,1,-1,1]*info.parent.r1*1.04);
title(ax,sprintf('Literal ring lattice with negative-x crack | N=64, N_r=46\n%d parent triangles split; geometry gate PASS', ...
    summary.splitParentTriangles));

figures(2)=new_figure('Step 4C: negative-x seam zoom',[1400,650],opts);
ax=axes(figures(2));
draw_mesh(ax,mesh,parent,splitChildren,opts.OverlayParent,true);
xlim(ax,[-0.08,-0.004]); ylim(ax,[-0.014,0.014]);
title(ax,{'Negative-x crack seam: local triangle splitting only', ...
    'Upper/lower node markers coincide exactly; different node IDs'});

figures(3)=new_figure('Step 4C: inner-boundary crack entry',[1100,850],opts);
ax=axes(figures(3));
draw_mesh(ax,mesh,parent,splitChildren,opts.OverlayParent,true);
r0=info.parent.r0;
xlim(ax,[-1.65,-0.70]*r0); ylim(ax,[-0.43,0.43]*r0);
title(ax,{'Crack enters the inner circular boundary', ...
    'Original inner-ring vertices and boundary chords are retained'});
upper=mesh.crackUpperT6IDs(:); lower=mesh.crackLowerT6IDs(:);
[~,entryIndex]=max(mesh.coord(upper,1));
entryUpper=upper(entryIndex);
% Match by coordinates rather than relying on face-array ordering.
entryXY=mesh.coord(entryUpper,:);
[~,lowerIndex]=min(sum((mesh.coord(lower,:)-entryXY).^2,2));
entryLower=lower(lowerIndex);
entryText=sprintf('Inner entry: (%.6g, 0)\nUpper ID %d; lower ID %d\nIdentical coordinates, separate topology', ...
    entryXY(1),entryUpper,entryLower);
text(ax,0.03,0.94,entryText,'Units','normalized', ...
    'VerticalAlignment','top','BackgroundColor','w','Margin',5,'FontSize',10);

outputDir=char(opts.OutputDir);
if ~isempty(outputDir)
    if ~isfolder(outputDir), mkdir(outputDir); end
    names={'literal_crack_cut_full.png','literal_crack_cut_seam.png', ...
        'literal_crack_cut_inner_entry.png'};
    for k=1:numel(figures)
        exportgraphics(figures(k),fullfile(outputDir,names{k}),'Resolution',200);
    end
    save(fullfile(outputDir,'literal_crack_cut_geometry.mat'), ...
        'mesh','info','parent','summary');
    jsonFile=fullfile(outputDir,'literal_crack_cut_audit.json');
    fid=fopen(jsonFile,'w');
    assert(fid>=0,'literalCutPreview:OutputFile','Cannot write %s.',jsonFile);
    closeFile=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));
    clear closeFile;
    fprintf('Geometry artifacts written to: %s\n',outputDir);
end
Out=struct('mesh',mesh,'info',info,'parent',parent, ...
    'summary',summary,'figures',figures,'outputDir',outputDir);
end

function f=new_figure(name,dimensions,opts)
f=figure('Name',name,'NumberTitle','off','Color','w', ...
    'Visible',char(opts.Visible),'Position',[80,80,dimensions]);
end

function draw_mesh(ax,mesh,parent,splitChildren,overlayParent,showNodes)
hold(ax,'on');
X=mesh.coord3; T=mesh.connect3;
if overlayParent
    triplot(parent.connect3,parent.coord3(:,1),parent.coord3(:,2), ...
        'Parent',ax,'Color',[0.83,0.83,0.83],'LineWidth',0.45);
end
triplot(T,X(:,1),X(:,2),'Parent',ax,'Color',[0.36,0.39,0.43],'LineWidth',0.50);
upperColor=[0.15,0.45,0.77]; lowerColor=[0.85,0.32,0.15];
centroidY=mean(reshape(X(T,2),size(T)),2);
draw_children(ax,T(splitChildren&centroidY>=0,:),X,upperColor);
draw_children(ax,T(splitChildren&centroidY<0,:),X,lowerColor);
hUpper=patch(ax,nan,nan,upperColor,'FaceAlpha',0.32,'EdgeColor',upperColor);
hLower=patch(ax,nan,nan,lowerColor,'FaceAlpha',0.32,'EdgeColor',lowerColor);
handles=[hUpper,hLower]; labels={'Split children: upper','Split children: lower'};
if showNodes
    up=mesh.crackUpperT6IDs(:); lo=mesh.crackLowerT6IDs(:);
    hUp=plot(ax,mesh.coord(up,1),mesh.coord(up,2),'o', ...
        'Color',upperColor,'MarkerSize',5.5,'LineWidth',0.8,'LineStyle','none');
    hLo=plot(ax,mesh.coord(lo,1),mesh.coord(lo,2),'x', ...
        'Color',lowerColor,'MarkerSize',4.5,'LineWidth',0.8,'LineStyle','none');
    handles=[handles,hUp,hLo];
    labels=[labels,{'Upper T6 face IDs','Lower T6 face IDs'}];
end
axis(ax,'equal'); grid(ax,'on'); box(ax,'on');
xlabel(ax,'x_1'); ylabel(ax,'x_2');
legend(ax,handles,labels,'Location','southoutside','Orientation','horizontal');
set(ax,'FontSize',11,'Layer','top');
end

function draw_children(ax,T,X,color)
if ~isempty(T)
    patch(ax,'Faces',T,'Vertices',X,'FaceColor',color, ...
        'FaceAlpha',0.32,'EdgeColor',color,'LineWidth',0.8);
end
end

function [area,quality,minAngle]=triangle_metrics(T,X)
p1=X(T(:,1),:); p2=X(T(:,2),:); p3=X(T(:,3),:);
v1=p2-p1; v2=p3-p1;
area=0.5*(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1));
a=sqrt(sum((p2-p3).^2,2));
b=sqrt(sum((p3-p1).^2,2));
c=sqrt(sum((p1-p2).^2,2));
quality=4*sqrt(3)*area./(a.^2+b.^2+c.^2);
cosAngles=[(b.^2+c.^2-a.^2)./(2*b.*c), ...
    (c.^2+a.^2-b.^2)./(2*c.*a),(a.^2+b.^2-c.^2)./(2*a.*b)];
minAngle=min(acosd(max(-1,min(1,cosAngles))),[],2);
end
