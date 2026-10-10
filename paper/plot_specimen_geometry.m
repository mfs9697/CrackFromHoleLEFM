function S=plot_specimen_geometry(varargin)
% Symbolic vector schematic, redrawn from the author's approved composition.
% Illustrative coordinates are drawing units, not FEM geometry or trajectory.
ip=inputParser;addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));parse(ip,varargin{:});
paper=fileparts(mfilename('fullpath'));out=char(ip.Results.OutputDir);
if isempty(out),out=fullfile(paper,'figures','geometry_loading');end
if exist(out,'dir')~=7,mkdir(out);end
P=publication_style();W=10;H=5.7;cx=3.15;cy=2.85;r=1.08;
fig=figure('Color','w','Visible','off','Units','centimeters','Position',[2 2 .78*P.textWidth_cm 9.2]);
cleanup=onCleanup(@()close(fig));ax=axes(fig);hold(ax,'on');
plot(ax,[0 W W 0 0],[0 0 H H 0],'k-','LineWidth',1);
t=linspace(0,2*pi,361);plot(ax,cx+r*cos(t),cy+r*sin(t),'k-','LineWidth',1);
plot(ax,[0 cx],[cy cy],'k--','LineWidth',.7);plot(ax,[cx cx],[0 H],'k--','LineWidth',.7);
plot(ax,cx,cy,'ko','MarkerFaceColor','k','MarkerSize',3);
point=[cx+r,cy];plot(ax,point(1),point(2),'ko','MarkerFaceColor','k','MarkerSize',3);
% Smooth parametric illustrative crack; never sampled from a solved path.
t=linspace(0,1,120);x=point(1)+(W-.8-point(1))*t;
y=cy+.7*t.*exp(-5*t)-.62*(1-exp(-5*t)).*t;
plot(ax,x,y,'k-','LineWidth',1.1);
label(7.3,cy-1.0,'Illustrative crack');
for x=linspace(.18,W-.18,15)
    arrow([x H],[x H+.68],.08);arrow([x 0],[x -.68],.08);
end
label(W/2,H+1.03,'$\sigma$');label(W/2,-.94,'$\sigma$');
label(-.58,H/2,'Free');label(W+.48,H/2,'Free');
plot(ax,[0 0],[-.05 -1.28],'k-','LineWidth',.6);plot(ax,[W W],[-.05 -1.28],'k-','LineWidth',.6);
double_arrow([0 -1.15],[W -1.15]);label(W/2,-1.49,'$2A$');
plot(ax,[W+.08 W+.9],[0 0],'k-','LineWidth',.6);plot(ax,[W+.08 W+.9],[H H],'k-','LineWidth',.6);
double_arrow([W+.82 0],[W+.82 H]);label(W+1.09,H/2,'$B$');
double_arrow([0 .54],[cx .54]);label(cx/2,.79,'$x_c$');
double_arrow([.48 0],[.48 cy]);label(.76,cy/2,'$y_c$');
arrow([cx cy],[cx+.72*r cy+.69*r],.08);label(cx+.48*r,cy+.48*r,'$R$');
label(cx,cy-.32,'$(x_c,y_c)$');label(point(1)+.36,point(2)+.37,'$P_0$');
plot(ax,0,0,'ko','MarkerFaceColor','k','MarkerSize',2.5);label(-.16,-.15,'$O$');
arrow([0 0],[.96 0],.08);arrow([0 0],[0 .94],.08);
label(1.08,.10,'$x$');label(-.16,1.05,'$y$');
axis(ax,'equal');xlim(ax,[-1.05 W+1.45]);ylim(ax,[-1.8 H+1.35]);axis(ax,'off');
set(ax,'Position',[.02 .02 .96 .96]);
for obj=findall(ax,'Type','text').',set(obj,'Interpreter','latex','FontSize',11);end
set(fig,'PaperUnits','centimeters','PaperSize',[.78*P.textWidth_cm 9.2], ...
    'PaperPosition',[0 0 .78*P.textWidth_cm 9.2],'PaperPositionMode','manual','Renderer','painters');
pdf=fullfile(out,'specimen_geometry.pdf');png=fullfile(out,'specimen_geometry.png');
print(fig,pdf,'-dpdf','-painters');print(fig,png,'-dpng','-r220');
S=struct('pdf',pdf,'vector',true,'illustrativeCrack',true,'numericalParametersShown',false, ...
    'origin','lower-left O','coordinateTransform','y_comp = y_sketch - B/2', ...
    'font_pt',11,'width_cm',.78*P.textWidth_cm);

    function label(x,y,s)
        text(ax,x,y,s,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',11);
    end
    function arrow(a,b,head)
        plot(ax,[a(1),b(1)],[a(2),b(2)],'k-','LineWidth',.75);
        v=(b-a)/norm(b-a);n=[-v(2),v(1)];base=b-2*head*v;
        patch(ax,[b(1),base(1)+head*n(1),base(1)-head*n(1)], ...
            [b(2),base(2)+head*n(2),base(2)-head*n(2)],'k','EdgeColor','none');
    end
    function double_arrow(a,b)
        arrow(a,b,.06);arrow(b,a,.06);
    end
end
