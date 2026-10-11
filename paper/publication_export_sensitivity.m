function info=publication_export_sensitivity(ax,pdf,png,kind)
% Fixed physical geometry for Figures 7--9 only. No numerical data changed.
paper=fileparts(mfilename('fullpath'));
L=jsondecode(fileread(fullfile(paper,'figure_7_9_layout.json')));C=L.(kind);
P=publication_style();fig=ancestor(ax,'figure');limits={ax.XLim,ax.YLim};
set(fig,'Units','centimeters','Position',[2 2 C.page(:).'/10],'Color','w');
set(ax,'Units','centimeters','Position',C.axes(:).'/10, ...
    'PositionConstraint','innerposition','FontName',P.fontName, ...
    'FontSize',L.tick_pt,'TickLabelInterpreter','latex','LineWidth',.8,'TickDir','out');
xlim(ax,limits{1});ylim(ax,limits{2});
if strcmp(kind,'top')
    % Common viewport and matched physical aspect ratio keep the hole circular.
    axis(ax,'equal');xlim(ax,[135 300]);ylim(ax,[-55 15]);
    set(ax,'Position',C.axes(:).'/10);
else
    ax.DataAspectRatioMode='auto';ax.PlotBoxAspectRatioMode='auto';
    ax.XLim=[0 100];ax.XTick=[0 50 100];
end
for obj=[ax.XLabel,ax.YLabel,ax.Title]
    set(obj,'Interpreter','latex','FontName',P.fontName,'FontSize',L.label_pt);
end
for obj=findall(ax,'Type','text').'
    set(obj,'Interpreter','latex','FontName',P.fontName);
end
try,ax.XAxis.SecondaryLabel.Interpreter='latex';ax.YAxis.SecondaryLabel.Interpreter='latex';catch,end
for obj=findall(ax,'Type','line').'
    if ~strcmp(obj.LineStyle,'none'),obj.LineWidth=P.lineWidth_pt;end
    if ~strcmp(obj.Marker,'none'),obj.MarkerSize=P.markerSize_pt;end
end
legends=findall(fig,'Type','legend');
for lg=legends(:).'
    set(lg,'Interpreter','latex','FontName',P.fontName,'FontSize',L.legend_pt, ...
        'Box','off','AutoUpdate','off','Units','centimeters');
    if strcmp(kind,'top'),lg.NumColumns=2;else,lg.NumColumns=1;end
    lg.ItemTokenSize=[16 8];
end
drawnow;
if strcmp(kind,'top')
    assert(isequal(ax.XLim,[135 300])&&isequal(ax.YLim,[-55 15]), ...
        'fig789:ChangedSpatialViewport','Spatial axes limits changed unexpectedly.');
end
for lg=legends(:).'
    pos=lg.Position;strip=C.legend(:).'/10;
    assert(pos(3)<=strip(3)&&pos(4)<=strip(4),'fig789:LegendTooLarge', ...
        'Legend %.3f x %.3f cm exceeds strip %.3f x %.3f cm: %s',pos(3),pos(4),strip(3),strip(4),pdf);
    lg.Position=[strip(1),strip(2)+(strip(4)-pos(4))/2,pos(3:4)];
end
set(fig,'PaperUnits','centimeters','PaperSize',C.page(:).'/10, ...
    'PaperPosition',[0 0 C.page(:).'/10],'PaperPositionMode','manual','Renderer','painters');
print(fig,pdf,'-dpdf','-painters');
publication_pdf_tool({'--exact-box',pdf,'--class',kind});
if nargin>=3&&~isempty(png),print(fig,png,'-dpng','-r220');end
info=struct('class',kind,'page_mm',C.page,'axes_mm',C.axes,'legendStrip_mm',C.legend, ...
    'legendPolicy',L.legend_policy,'dataChanged',false,'xlim',ax.XLim,'ylim',ax.YLim);
info.legendBounds_mm=[];
for lg=legends(:).',info.legendBounds_mm(end+1,:)=lg.Position*10;end
fid=fopen(strrep(pdf,'.pdf','.layout.json'),'w');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(info));fclose(fid);
end
