function S=plot_specimen_geometry(varargin)
% Symbolic vector schematic: specimen geometry, tractions and possible crack.
% This function is PLOT ONLY. All coordinates below are drawing coordinates.
%
% The original publication math style uses MATLAB's LaTeX interpreter.
% Restore that style as the default: consistent italic math, subscripts,
% and Greek sigma, without experimental TeX font-name overrides.
% Optional Euclid mode is retained for comparison only. Its special
% glyph behavior is NOT equivalent to the original approved appearance.
%
% Call from repository root:
%   addpath(genpath(pwd));
%   S=plot_specimen_geometry();              % standard pale-gray version
%   S=plot_specimen_geometry('PlateGray',.92); % optional slightly darker fill
%   winopen(S.pdf);                          % show PDF on Windows

ip=inputParser;
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'FontMode','latex',@(x)ischar(x)||isstring(x));
addParameter(ip,'PlateGray',0.94,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
parse(ip,varargin{:});
fontMode=lower(char(ip.Results.FontMode));
assert(ismember(fontMode,{'euclid','latex'}), ...
    'specimenFig:FontMode','FontMode must be ''euclid'' or ''latex''.');

paper=fileparts(mfilename('fullpath'));
out=char(ip.Results.OutputDir);
if isempty(out),out=fullfile(paper,'figures','geometry_loading');end

% Only explicitly requested Euclid mode requires local font validation.
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

% Light-gray plate with a white, open hole; all mesh-free linework
% remains black. Do not change the numerical model or the schematic geometry.
plateGray=double(ip.Results.PlateGray)*[1 1 1];
patch(ax,[0 W W 0],[0 0 H H],plateGray, ...
    'EdgeColor','k','LineWidth',1.35);
a=linspace(0,2*pi,361);
patch(ax,cx+r*cos(a),cy+r*sin(a),'w', ...
    'EdgeColor','k','LineWidth',1.35);
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
    'mathFont',fontDescription(),'plateGray',ip.Results.PlateGray,'png',png, ...
    'mathItalicEnabled',true);
fprintf('Figure 1 regenerated (%s):\n  %s\n  %s\n',fontMode,pdf,png);

    function name=fontDescription()
        if strcmp(fontMode,'euclid')
            name='Euclid / Euclid Symbol';
        else
            name='MATLAB LaTeX (original mathematical lettering)';
        end
    end

    function t=label(x,y,mathString)
        % MATLAB TeX recognizes \fontname and \it; MATLAB LaTeX does not
        % respect FontName. Define the exact mathematical typography rather
        % than printing Euclid variables in the upright roman face.
        if strcmp(fontMode,'euclid')
            switch mathString
                case '$P_0$'
                    glyph='\fontname{Euclid}\it P_{\rm 0}';
                case '$R$'
                    glyph='\fontname{Euclid}\it R';
                case '$B$'
                    glyph='\fontname{Euclid}\it B';
                case '$x$'
                    glyph='\fontname{Euclid}\it x';
                case '$y$'
                    glyph='\fontname{Euclid}\it y';
                case '$x_c$'
                    glyph='\fontname{Euclid}\it x_{c}';
                case '$y_c$'
                    glyph='\fontname{Euclid}\it y_{c}';
                case '$2A$'
                    glyph='\fontname{Euclid}\rm 2\it A';
                case '$(x_c,y_c)$'
                    glyph='\fontname{Euclid}\rm (\it x_{c}\rm ,\it y_{c}\rm )';
                case '$O$'
                    glyph='\fontname{Euclid}\rm O';
                otherwise
                    error('specimenFig:UnknownMathLabel', ...
                        'No Euclid math style has been defined for %s.',mathString);
            end
            t=text(ax,x,y,glyph,'Interpreter','tex', ...
                'FontName','Euclid','FontSize',labelPt, ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        else
            % Original Figure 1 math: LaTeX controls all italic and symbolic text.
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
            % Upright words rendered by the same LaTeX interpreter.
            text(ax,x,y,s,'Interpreter','latex', ...
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
