function S=publication_m0_label_stamps(output)
% Display-only vector labels for retained M1 graphics. No numerical history.
if ~isfolder(output),mkdir(output);end
labels={'$y_{\mathrm{M1}}-y_{\mathrm{M0}}$ [$\mu$m]', ...
    '$\theta_{\mathrm{M1}}-\theta_{\mathrm{M0}}$ [mdeg]'};
names={'m1_vertical_label','m1_direction_label'};S=struct('numericalDataRead',false,'files',{{}});
for k=1:2
    fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 2 6.4]);
    cleanup=onCleanup(@()close(fig));
    ax=axes(fig,'Position',[0 0 1 1],'Visible','off','XLim',[0 2],'YLim',[0 6.4]);
    text(ax,1,3.2,labels{k},'Rotation',90,'Interpreter','latex','FontSize',11, ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Color',[.149 .149 .149]);
    set(fig,'PaperUnits','centimeters','PaperSize',[2 6.4],'PaperPosition',[0 0 2 6.4], ...
        'PaperPositionMode','manual','Renderer','painters');
    file=fullfile(output,[names{k} '.pdf']);print(fig,file,'-dpdf','-painters');
    S.files{end+1}=file;clear cleanup
end
end
