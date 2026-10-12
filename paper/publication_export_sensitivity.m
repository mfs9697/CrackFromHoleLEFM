function info=publication_export_sensitivity(ax,pdf,png,kind,deferCrop)
% Fixed physical geometry for Figures 7--9 only. No numerical data changed.
paper=fileparts(mfilename('fullpath'));
L=jsondecode(fileread(fullfile(paper,'figure_7_9_layout.json')));C=L.(kind);
P=publication_style();fig=ancestor(ax,'figure');limits={ax.XLim,ax.YLim};
if nargin<5,deferCrop=false;end
set(fig,'Units','centimeters','Position',[2 2 C.page(:).'/10],'Color','w');
set(ax,'Units','centimeters','Position',C.axes(:).'/10, ...
    'PositionConstraint','innerposition','FontName',P.fontName, ...
    'FontSize',L.tick_pt,'TickLabelInterpreter','latex','LineWidth',.8,'TickDir','out','Box','on');
xlim(ax,limits{1});ylim(ax,limits{2});
if strcmp(kind,'top')
    % Common viewport and matched physical aspect ratio keep the hole circular.
    axis(ax,'equal');xlim(ax,limits{1});ylim(ax,limits{2});
    set(ax,'Position',C.axes(:).'/10);
else
    ax.DataAspectRatioMode='auto';ax.PlotBoxAspectRatioMode='auto';
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
    lg.NumColumns=1;
    lg.ItemTokenSize=[16 8];
    if contains(pdf,'tip_resolution'),lg.Title.String='M0 exterior';lg.Title.Interpreter='latex';end
    if contains(pdf,'increment')&&contains(pdf,{'trajectory','turn_density'}),lg.Title.String='Scaled M1';end
end
drawnow;
assert(isequal(ax.XLim,limits{1})&&isequal(ax.YLim,limits{2}), ...
    'fig789:ChangedSpatialViewport','Plot axis limits changed unexpectedly.');
for lg=legends(:).'
    publication_place_internal_legend(ax,lg);
end
set(fig,'PaperUnits','centimeters','PaperSize',C.page(:).'/10, ...
    'PaperPosition',[0 0 C.page(:).'/10],'PaperPositionMode','manual','Renderer','painters');
print(fig,pdf,'-dpdf','-painters');
publication_pdf_tool({'--exact-box',pdf,'--class',kind});
info=struct('class',kind,'page_mm',C.page,'axes_mm',C.axes, ...
    'legendPolicy',L.legend_policy,'dataChanged',false,'xlim',ax.XLim,'ylim',ax.YLim);
info.legendBounds_mm=[];
for lg=legends(:).',info.legendBounds_mm(end+1,:)=lg.Position*10;end
fid=fopen(strrep(pdf,'.pdf','.layout.json'),'w');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(info));fclose(fid);
if ~deferCrop,publication_pdf_tool({'--crop',pdf});info=jsondecode(fileread(strrep(pdf,'.pdf','.layout.json')));
elseif nargin>=3&&~isempty(png),print(fig,png,'-dpng','-r220');end
end
