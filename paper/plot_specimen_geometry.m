function S=plot_specimen_geometry(varargin)
% Symbolic, publication-style vector schematic of the specimen and loading.
% Drawing coordinates are illustrative, not FEM geometry or a solved trajectory.
% All figure labels use the manuscript's shared publication font size.
ip=inputParser;
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});
paper=fileparts(mfilename('fullpath'));
out=char(ip.Results.OutputDir);
if isempty(out),out=fullfile(paper,'figures','geometry_loading');end
if exist(out,'dir')~=7,mkdir(out);end

P=publication_style();
W=10; H=5.7; cx=3.15; cy=2.85; r=1.08;
pageWidth=.78*P.textWidth_cm;
pageHeight=9.2;
fig=figure('Color','w','Visible','off','Units','centimeters', ...
    'Position',[2 2 pageWidth pageHeight]);
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
ax=axes(fig);hold(ax,'on');

% Plate, circular hole, and symbolic construction lines.
plot(ax,[0 W W 0 0],[0 0 H H 0],'k-','LineWidth',P.lineWidth_pt);
circleAngle=linspace(0,2*pi,361);
plot(ax,cx+r*cos(circleAngle),cy+r*sin(circleAngle),'k-', ...
    'LineWidth',P.lineWidth_pt);
plot(ax,[0 cx],[cy cy],'k--','LineWidth',.7);
plot(ax,[cx cx],[0 H],'k--','LineWidth',.7);
plot(ax,cx,cy,'ko','MarkerFaceColor','k','MarkerSize',P.markerSize_pt);

% The crack leaves P0 tangentially and bends smoothly downward.
% This symbolic curve is independent of all accepted numerical crack states.
P0=[cx+r,cy];
plot(ax,P0(1),P0(2),'ko','MarkerFaceColor','k', ...
    'MarkerSize',P.markerSize_pt);
u=linspace(0,1,220);
xCr=P0(1)+(W-.85-P0(1))*u;
smoothStep=3*u.^2-2*u.^3; % zero slope at both ends
yCr=P0(2)+.37*sin(pi*u).^2-.92*smoothStep;
plot(ax,xCr,yCr,'k-','LineWidth',1.1);
label(P0(1)+.38,P0(2)+.38,'$P_0$');

% Tensile tractions are drawn outward and labeled by unsigned sigma.
for xArrow=linspace(.18,W-.18,15)
    arrow([xArrow H],[xArrow H+.68],.08);
    arrow([xArrow 0],[xArrow -.68],.08);
end
label(W/2,H+1.05,'$\sigma$');
label(W/2,-.92,'$\sigma$');

% Place both "Free" labels inside the plate, clear of the hole-center
% construction line and of the outer right-side B dimension.
label(.47,H/2+.62,'Free');
label(W-.47,H/2+.62,'Free');

% Overall plate height B: its dimension is the outermost right annotation.
plot(ax,[W+.08 W+1.10],[0 0],'k-','LineWidth',.6);
plot(ax,[W+.08 W+1.10],[H H],'k-','LineWidth',.6);
double_arrow([W+1.00 0],[W+1.00 H]);
label(W+1.38,H/2,'$B$');

% Hole-center ordinate y_c is measured from the bottom edge and dimensioned
% outside the left boundary, away from the coordinate axes and "Free" label.
plot(ax,[-.80 0],[0 0],'k-','LineWidth',.6);
plot(ax,[-.80 0],[cy cy],'k-','LineWidth',.6);
double_arrow([-.63 0],[-.63 cy]);
label(-.90,cy/2,'$y_c$');

% Bottom annotations occupy distinct levels: loading arrows, x_c, then 2A.
% The x_c dimension uses the left edge and the hole-center projection.
plot(ax,[0 0],[-.84 -2.17],'k-','LineWidth',.6);
plot(ax,[cx cx],[-.84 -1.51],'k-','LineWidth',.6);
double_arrow([0 -1.36],[cx -1.36]);
label(cx/2,-1.15,'$x_c$');
plot(ax,[W W],[-.84 -2.17],'k-','LineWidth',.6);
double_arrow([0 -2.06],[W -2.06]);
label(W/2,-2.31,'$2A$');

% Radius arrow and nonintersecting R label. The offset is PERPENDICULAR to
% the radial shaft; do not return to a label position on the arrow.
radialUnit=[.76 sqrt(1-.76^2)];
radiusEnd=[cx cy]+r*radialUnit;
arrow([cx cy],radiusEnd,.08);
radiusMid=([cx cy]+radiusEnd)/2;
labelNormal=[radialUnit(2) -radialUnit(1)];
rLabelPos=radiusMid+.42*labelNormal;
assert(abs(dot(rLabelPos-radiusMid,labelNormal))>.35, ...
    'specimenFig:RadiusLabelClearance');
label(rLabelPos(1),rLabelPos(2),'$R$');
label(cx-.08,cy-.40,'$(x_c,y_c)$');

% Lower-left physical origin: short Cartesian arrows share the plate edges
% but do not intersect either center-coordinate dimension.
plot(ax,0,0,'ko','MarkerFaceColor','k','MarkerSize',2.5);
arrow([0 0],[1.11 0],.08);
arrow([0 0],[0 1.09],.08);
label(-.16,-.20,'$O$');
label(1.28,.15,'$x$');
label(-.19,1.22,'$y$');

axis(ax,'equal');
xlim(ax,[-1.32 W+1.77]);
ylim(ax,[-2.60 H+1.39]);
axis(ax,'off');
set(ax,'Position',[.02 .02 .96 .96]);
for obj=findall(ax,'Type','text').'
    set(obj,'Interpreter','latex','FontSize',P.label_pt);
end

set(fig,'PaperUnits','centimeters','PaperSize',[pageWidth pageHeight], ...
    'PaperPosition',[0 0 pageWidth pageHeight], ...
    'PaperPositionMode','manual','Renderer','painters');
pdf=fullfile(out,'specimen_geometry.pdf');
png=fullfile(out,'specimen_geometry.png');
print(fig,pdf,'-dpdf','-painters');
print(fig,png,'-dpng','-r220');
S=struct('pdf',pdf,'vector',true,'illustrativeCrack',true, ...
    'numericalParametersShown',false,'origin','lower-left O', ...
    'coordinateTransform','y_comp = y_sketch - B/2', ...
    'font_pt',P.label_pt,'width_cm',pageWidth);

    function label(x,y,s)
        text(ax,x,y,s,'Interpreter','latex', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',P.label_pt);
    end

    function arrow(a,b,head)
        displacement=b-a;
        direction=displacement/norm(displacement);
        normal=[-direction(2) direction(1)];
        base=b-2*head*direction;
        plot(ax,[a(1) b(1)],[a(2) b(2)],'k-','LineWidth',.75);
        patch(ax,[b(1),base(1)+head*normal(1),base(1)-head*normal(1)], ...
            [b(2),base(2)+head*normal(2),base(2)-head*normal(2)], ...
            'k','EdgeColor','none');
    end

    function double_arrow(a,b)
        arrow(a,b,.06);
        arrow(b,a,.06);
    end
end
