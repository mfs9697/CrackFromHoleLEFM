function S=plot_saved_increment_records(file,out)
% Redraw three panels from saved plotting records. Native trajectories are
% NOT reconstructed from the common-length samples; trajectory.pdf retained.
z=load(file,'Increment');I=z.Increment;T=I.metrics;N=I.native;
assert(I.has1mm&&I.maxCommonCrackLength_mm==92&&height(T)==23&&height(N)==161, ...
    'paperfig:IncompleteIncrementRecords','Complete exported 4/2/1-mm records required.');
if exist(out,'dir')~=7,mkdir(out);end
cols={[.85 .325 .098],[.929 .694 .125]};names={'vertical_deviation','direction_deviation','turn_density'};
for j=1:3
    fig=figure('Visible','off','Color','w');cleanup=onCleanup(@()close(fig));ax=axes(fig);hold(ax,'on');
    if j==1
        fields={'dy_2minus4_um','dy_1minus2_um'};lab='Successive $\delta y$ [$\mu$m]';
    elseif j==2
        fields={'dtheta_2minus4_mdeg','dtheta_1minus2_mdeg'};lab='Direction difference [mdeg]';
    else
        fields={};lab='$\Delta\theta/\Delta a$ [deg/mm]';
    end
    if j<3
        marks={'o','^'};labels={'$2-4$','$1-2$'};
        for k=1:2,plot(ax,T.crack_length_mm,T.(fields{k}),['-' marks{k}], ...
                'Color',cols{k},'DisplayName',labels{k});end
        yline(ax,0,':','HandleVisibility','off');
    else
        colors={[0 .447 .741],cols{1},cols{2}};marks={'o','s','^'};
        for k=1:3
            increment=[4,2,1];hit=N.increment_mm==increment(k);
            plot(ax,N.crack_length_mm(hit),N.turn_density_deg_per_mm(hit),['-' marks{k}], ...
                'Color',colors{k},'DisplayName',['$' num2str(increment(k)) '$ mm']);
        end
    end
    xlim(ax,[0,100]);xlabel(ax,'$a$ [mm]','Interpreter','latex');ylabel(ax,lab,'Interpreter','latex');
    if j==1,location='northeast';else,location='northwest';end
    legend(ax,'Location',location,'Interpreter','latex','Box','off');box(ax,'on');
    publication_export_axis(ax,fullfile(out,[names{j} '.pdf']),fullfile(out,[names{j} '.png']),'third');
    clear cleanup
end
S=struct('status','regenerated from saved plot records','source',file,'regeneratedPanels',{names}, ...
    'trajectoryStatus','retained: native vertex histories unavailable','numericalExperiments',0);
end
