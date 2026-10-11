function info=publication_normalize_retained_pdf(key,target)
% Wrap original vector evidence onto a shared blank axes plate. No trajectory
% history is read/reconstructed; no solver, EDI, COD or mesher is called.
paper=fileparts(mfilename('fullpath'));
profiles=jsondecode(fileread(fullfile(paper,'data','retained_figure_7_9','profiles.json')));
p=profiles.(key);fig=figure('Visible','off','Color','w');
plate=[tempname '.pdf'];cleanup=onCleanup(@()local_cleanup(fig,plate));
ax=axes(fig);box(ax,'on');ylim(ax,p.ylim(:).');
if strcmp(p.class,'top')
    xlim(ax,[135 300]);xlabel(ax,'$x$ [mm]','Interpreter','latex');
else
    xlim(ax,[0 100]);xlabel(ax,'$a$ [mm]','Interpreter','latex');
end
ylabel(ax,p.ylabel,'Interpreter','latex');
publication_export_sensitivity(ax,plate,'',p.class);
publication_pdf_tool({'--key',key,'--plate',plate,'--output',target});
info=jsondecode(fileread(strrep(target,'.pdf','.layout.json')));
end
function local_cleanup(fig,plate)
if isgraphics(fig),close(fig);end
if isfile(plate),delete(plate);end
metadata=strrep(plate,'.pdf','.layout.json');if isfile(metadata),delete(metadata);end
end
