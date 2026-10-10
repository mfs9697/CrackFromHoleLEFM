function S=plot_specimen_geometry(varargin)
% Symbolic vector schematic: specimen geometry, tractions and possible crack.
% This function is PLOT ONLY. All coordinates below are drawing coordinates.
%
% FontMode='euclid' (default) requires installed, licensed MathType Euclid and
% Euclid Symbol fonts. MATLAB's 'latex' interpreter cannot select arbitrary
% math fonts, so this mode uses TeX for Latin symbols and the Euclid Symbol
% glyph for sigma. Do NOT silently substitute Computer Modern or other fonts.
% FontMode='latex' is an explicit portability/proofing fallback only.
%
% Call from repository root:
%   addpath(genpath(pwd));
%   S=plot_specimen_geometry();

ip=inputParser;
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'FontMode','euclid',@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});
fontMode=lower(char(ip.Results.FontMode));
assert(ismember(fontMode,{'euclid','latex'}), ...
    'specimenFig:FontMode','FontMode must be ''euclid'' or ''latex''.');

paper=fileparts(mfilename('fullpath'));
out=char(ip.Results.OutputDir);
if isempty(out),out=fullfile(paper,'figures','geometry_loading');end

% Validate typography BEFORE touching the publication image assets.
if strcmp(fontMode,'euclid')
    installed=listfonts;
    haveEuclid=any(strcmpi(installed,'Euclid'));
    haveSymbol=any(strcmpi(installed,'Euclid Symbol'));
    if ~(haveEuclid && haveSymbol)
        error('specimenFig:EuclidFontsMissing', ...
            ['The licensed Euclid and Euclid Symbol fonts must be installed ' ...
             'for Figure 1. MATLAB''s LaTeX interpreter uses its own math ' ...
             'font and is not a Euclid substitute. For a NON-PUBLICATION ' ...
             'proof only, use plot_specimen_geometry(''FontMode'',''latex'').']);
    end
end
if exist(out,'dir')~=7,mkdir(out);end

P=publication_style();
W=10; H=5.7; cx=3.15; cy=2.85; r=1.08;
pageWidth=.78*P.textWidth_cm;
pageHeight=9.2;
% A slightly larger label size, with a smaller plate-to-letter ratio
% than the initial schematic. The export page and LaTeX width are unchanged.
labelPt=12;
fig=figure('Color','w','Visible','off','Units','centimeters', ...
    'Position',[2 2 pageWidth pageHeight]);
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
ax=axes(fig);hold(ax,'on');

% Plate, hole and auxiliary center lines.
plot(ax,[0 W W 0 0],[0 0 H H 0],'k-','LineWidth',1.35);
a=linspace(0,2*pi,361);
plot(ax,cx+r*cos(a),cy+r*sin(a),'k-','LineWidth',1.35);
plot(ax,[0 cx],[cy cy],'k--','LineWidth',.85);
plot(ax,[cx cx],[0 H],'k--','LineWidth',.85);
plot(ax,cx,cy,'ko','MarkerFaceColor','k','MarkerSize',3.5);

% Schematic path: horizontal initial tangent, smooth downward bend,
% horizontal final tangent. Never substitute a calculated trajectory here.
P0=[cx+r,cy];
plot(ax,P0(1),P0(2),'ko','MarkerFaceColor','k','MarkerSize',3.5);
u=linspace(0,1,220);
smoothStep=3*u.^2-2*u.^3;
xCr=P0(1)+(W-.85-P0(1))*u;
yCr=P0(2)+.37*sin(pi*u).^2-.92*smoothStep;
plot(ax,xCr,yCr,'k-','LineWidth',1.5);
label(P0(1)+.38,P0(2)+.42,'$P_0$');

% Fewer, longer and thicker arrows with larger arrowheads.
for xArrow=linspace(.30,W-.30,12)
    arrow([xArrow H],[xArrow H+.83],.15,1.05);
    arrow([xArrow 0],[xArrow -.83],.15,1.05);
end
stressLabel(W/2,H+1.18);
stressLabel(W/2,-1.12);

% Free-edge captions stay well INSIDE the plate, away from both dimensions.
plainLabel(.92,H/2+.70,'Free');
plainLabel(W-.92,H/2+.70,'Free');

% Overall plate height B, outside the right edge.
plot(ax,[W+.12 W+1.14],[0 0],'k-','LineWidth',.8);
plot(ax,[W+.12 W+1.14],[H H],'k-','LineWidth',.8);
doubleArrow([W+1.04 0],[W+1.04 H],.12);
label(W+1.47,H/2,'$B$');

% Center ordinate y_c, from the bottom edge.
plot(ax,[-.85 0],[0 0],'k-','LineWidth',.8);
plot(ax,[-.85 0],[cy cy],'k-','LineWidth',.8);
doubleArrow([-.69 0],[-.69 cy],.12);
label(-1.04,cy/2,'$y_c$');

% Coordinate dimensions below the traction arrows, at separate levels.
plot(ax,[0 0],[-.91 -2.35],'k-','LineWidth',.8);
plot(ax,[cx cx],[-.91 -1.62],'k-','LineWidth',.8);
doubleArrow([0 -1.48],[cx -1.48],.12);
label(cx/2,-1.25,'$x_c$');
plot(ax,[W W],[-.91 -2.35],'k-','LineWidth',.8);
doubleArrow([0 -2.24],[W -2.24],.12);
label(W/2,-2.49,'$2A$');

% Radius: its math label is offset PERPENDICULAR to the arrow, not over it.
radialUnit=[.76 sqrt(1-.76^2)];
radiusEnd=[cx cy]+r*radialUnit;
arrow([cx cy],radiusEnd,.14,.9);
radiusMid=([cx cy]+radiusEnd)/2;
normalUnit=[radialUnit(2) -radialUnit(1)];
rLabel=radiusMid+.54*normalUnit;
assert(norm(rLabel-radiusMid)>.50,'specimenFig:RadiusClearance');
label(rLabel(1),rLabel(2),'$R$');
label(cx-.05,cy-.48,'$(x_c,y_c)$');

% Global schematic axes start at the lower-left plate corner.
plot(ax,0,0,'ko','MarkerFaceColor','k','MarkerSize',3);
arrow([0 0],[1.14 0],.15,.95);
arrow([0 0],[0 1.19],.15,.95);
label(-.17,-.23,'$O$');
label(1.34,.27,'$x$');
label(-.18,1.39,'$y$');

% Preserve vector geometry and reserve physical space for the larger letters.
axis(ax,'equal');
xlim(ax,[-2.07 W+2.25]);
ylim(ax,[-2.88 H+1.60]);
axis(ax,'off');
set(ax,'Position',[.015 .015 .97 .97]);

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
    'fontMode',fontMode,'font_pt',labelPt,'width_cm',pageWidth, ...
    'mathFont',fontDescription());

    function name=fontDescription()
        if strcmp(fontMode,'euclid')
            name='Euclid / Euclid Symbol';
        else
            name='MATLAB LaTeX (fallback, not Euclid)';
        end
    end

    function t=label(x,y,mathString)
        % Preserve subscript syntax such as x_c and P_0.
        if strcmp(fontMode,'euclid')
            mathString=regexprep(mathString,'^\$|\$$','');
            t=text(ax,x,y,mathString,'Interpreter','tex', ...
                'FontName','Euclid','FontSize',labelPt, ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        else
            t=text(ax,x,y,mathString,'Interpreter','latex', ...
                'FontSize',labelPt,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle');
        end
    end

    function plainLabel(x,y,s)
        if strcmp(fontMode,'euclid')
            text(ax,x,y,s,'Interpreter','none','FontName','Euclid', ...
                'FontSize',labelPt,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle');
        else
            text(ax,x,y,s,'Interpreter','none', ...
                'FontSize',labelPt,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle');
        end
    end

    function stressLabel(x,y)
        if strcmp(fontMode,'euclid')
            % Euclid Symbol uses the classic Symbol encoding: ASCII 's'
            % represents the lowercase sigma. Inspect the actual PDF glyph
            % on the workstation; this is not a Unicode substitution.
            text(ax,x,y,'s','Interpreter','none','FontName','Euclid Symbol', ...
                'FontSize',labelPt+1,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle');
        else
            label(x,y,'$\sigma$');
        end
    end

    function arrow(p1,p2,head,lw)
        d=p2-p1;unit=d/norm(d);perp=[-unit(2),unit(1)];
        foot=p2-2*head*unit;
        plot(ax,[p1(1) p2(1)],[p1(2) p2(2)],'k-','LineWidth',lw);
        patch(ax,[p2(1) foot(1)+.83*head*perp(1) foot(1)-.83*head*perp(1)], ...
            [p2(2) foot(2)+.83*head*perp(2) foot(2)-.83*head*perp(2)], ...
            'k','EdgeColor','none');
    end

    function doubleArrow(p1,p2,head)
        arrow(p1,p2,head,.8);
        arrow(p2,p1,head,.8);
    end
end
