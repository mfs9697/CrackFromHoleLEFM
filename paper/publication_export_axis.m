function info=publication_export_axis(ax,pdf,png,kind)
% Export vector graphics on a fixed physical page; data/axis limits unchanged.
P=publication_style();
if nargin<4,kind='half';end
if strcmp(kind,'full'),width=P.fullWidth_cm;height=4.4;
elseif strcmp(kind,'spatial_full'),width=P.fullWidth_cm;height=6.2;
elseif strcmp(kind,'third'),width=P.thirdWidth_cm;height=6.6;
else,width=P.halfWidth_cm;height=5.8;end
fig=ancestor(ax,'figure');limits={ax.XLim,ax.YLim};
set(fig,'Units','centimeters','Position',[2 2 width height],'Color','w');
set(ax,'Units','normalized','Position',[.20 .23 .76 .70], ...
    'FontName',P.fontName,'FontSize',P.tick_pt,'TickLabelInterpreter','latex', ...
    'LineWidth',.8,'TickDir','out');
if any(strcmp(kind,{'full','spatial_full'})),ax.Position=[.09 .24 .88 .69];end
xlim(ax,limits{1});ylim(ax,limits{2});
for obj=[ax.XLabel,ax.YLabel,ax.ZLabel,ax.Title]
    set(obj,'Interpreter','latex','FontSize',P.label_pt,'FontName',P.fontName);
end
texts=findall(ax,'Type','text');
for obj=texts(:).'
    set(obj,'Interpreter','latex','FontName',P.fontName);
    if ~ismember(obj,[ax.XLabel,ax.YLabel,ax.ZLabel,ax.Title]),obj.FontSize=P.annotation_pt;end
end
for lg=findall(fig,'Type','legend').'
    set(lg,'Interpreter','latex','FontName',P.fontName,'FontSize',P.legend_pt);
    if numel(lg.String)==4,lg.NumColumns=1;lg.Location='northwest';end
end
for line=findall(ax,'Type','line').'
    if ~strcmp(line.LineStyle,'none'),line.LineWidth=P.lineWidth_pt;end
    if ~strcmp(line.Marker,'none'),line.MarkerSize=P.markerSize_pt;end
end
try,ax.XAxis.SecondaryLabel.Interpreter='latex';ax.YAxis.SecondaryLabel.Interpreter='latex';catch,end
set(fig,'PaperUnits','centimeters','PaperSize',[width height], ...
    'PaperPosition',[0 0 width height],'PaperPositionMode','manual','Renderer','painters');
drawnow;
print(fig,pdf,'-dpdf','-painters');
if nargin>=3&&~isempty(png),print(fig,png,'-dpng','-r220');end
info=struct('pdf',pdf,'width_cm',width,'height_cm',height, ...
    'label_pt',P.label_pt,'tick_pt',P.tick_pt,'legend_pt',P.legend_pt, ...
    'interpreter','latex','fixedPageExport',true);
end
