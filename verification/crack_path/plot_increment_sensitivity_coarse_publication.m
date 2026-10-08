function S = plot_increment_sensitivity_coarse_publication(varargin)
%PLOT_INCREMENT_SENSITIVITY_COARSE_PUBLICATION
% Compare independent Delta a = 4 mm and 2 mm trajectories produced by
% main_increment_sensitivity_coarse on the same dimensionless coarse family.
%
% Panels:
%   (a) complete crack trajectories;
%   (b) matched-length vertical coordinate difference, 2 mm - 4 mm;
%   (c) matched-length absolute-direction difference;
%   (d) q_K = KII/KI histories on their native sampling.
%
% The two curves use the same dimensionless M1 exterior + 2h0 tip family.
% No interpolation is used for panels (b,c): comparison points are exact
% common crack lengths 4,8,... mm.

    ip=inputParser;
    addParameter(ip,'Run4Dir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Run2Dir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Export',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ExportCombined',false,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));

    run4=char(opt.Run4Dir);
    if isempty(run4)
        run4=fullfile(root,'verification','crack_path', ...
            'increment_sensitivity_coarse','da_4mm');
    elseif ~local_is_absolute_path(run4)
        run4=fullfile(root,run4);
    end

    run2=char(opt.Run2Dir);
    if isempty(run2)
        run2=fullfile(root,'verification','crack_path', ...
            'increment_sensitivity_coarse','da_2mm');
    elseif ~local_is_absolute_path(run2)
        run2=fullfile(root,run2);
    end

    outDir=char(opt.OutputDir);
    if isempty(outDir)
        outDir=fullfile(root,'paper','figures','increment_sensitivity');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    f4=fullfile(run4,'states.csv');
    f2=fullfile(run2,'states.csv');
    v4f=fullfile(run4,'vertices.csv');
    v2f=fullfile(run2,'vertices.csv');
    assert(exist(f4,'file')==2,'incfig:Missing4mm', ...
        '4-mm states not found: %s',f4);
    assert(exist(f2,'file')==2,'incfig:Missing2mm', ...
        '2-mm states not found: %s',f2);
    assert(exist(v4f,'file')==2&&exist(v2f,'file')==2, ...
        'incfig:MissingVertices','Both vertices.csv files are required.');

    T4=sortrows(readtable(f4),'crack_length_mm');
    T2=sortrows(readtable(f2),'crack_length_mm');
    V4=sortrows(readtable(v4f),'crack_length_mm');
    V2=sortrows(readtable(v2f),'crack_length_mm');

    req={'segment','tip_x_m','tip_y_m','theta_deg','KI_unit','KII_unit', ...
        'KII_over_KI','delta_theta_next_deg','theta_next_deg','crack_length_mm'};
    assert(all(ismember(req,T4.Properties.VariableNames))&& ...
           all(ismember(req,T2.Properties.VariableNames)), ...
        'incfig:Columns','Increment-study tables lack required columns.');

    maxCommon=min(max(T4.crack_length_mm),max(T2.crack_length_mm));
    T4c=T4(T4.crack_length_mm<=maxCommon+1e-9,:);
    T2c=T2(T2.crack_length_mm<=maxCommon+1e-9,:);

    % Exact common lengths. The 4-mm states are a subset of the 2-mm grid.
    a4=T4c.crack_length_mm;
    a2=T2c.crack_length_mm;
    common=a4(ismembertol(a4,a2,1e-10,'DataScale',1));
    assert(~isempty(common),'incfig:NoCommonLengths','No matched crack lengths.');
    A4=local_rows_at_lengths(T4c,common);
    A2=local_rows_at_lengths(T2c,common);

    dx_um=1e6*(A2.tip_x_m-A4.tip_x_m);
    dy_um=1e6*(A2.tip_y_m-A4.tip_y_m);
    dr_um=hypot(dx_um,dy_um);
    dtheta_mdeg=1e3*(A2.theta_deg-A4.theta_deg);
    dq=A2.KII_over_KI-A4.KII_over_KI;

    [qmax4,i4]=max(T4c.KII_over_KI);
    [qmax2,i2]=max(T2c.KII_over_KI);
    aQ4=T4c.crack_length_mm(i4);
    aQ2=T2c.crack_length_mm(i2);

    [cross4,aLS4]=local_zero_crossing(T4c);
    [cross2,aLS2]=local_zero_crossing(T2c);

    % Actual trajectories.
    V4=V4(V4.crack_length_mm<=maxCommon+1e-9,:);
    V2=V2(V2.crack_length_mm<=maxCommon+1e-9,:);
    X4=1e3*V4.x_m; Y4=1e3*V4.y_m;
    X2=1e3*V2.x_m; Y2=1e3*V2.y_m;

    holeCenter=[170,-20];
    holeRadius=30;
    ang=linspace(0,2*pi,361);
    hx=holeCenter(1)+holeRadius*cos(ang);
    hy=holeCenter(2)+holeRadius*sin(ang);

    blue=[0.0000 0.4470 0.7410];
    orange=[0.8500 0.3250 0.0980];
    gray=[0.45 0.45 0.45];

    fig=figure('Color','w','Units','centimeters','Position',[2 2 18.2 18.0]);
    tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

    ax1=nexttile(tl,[1 3]); hold(ax1,'on');
    plot(ax1,hx,hy,'-','Color',gray,'LineWidth',.75,'HandleVisibility','off');
    plot(ax1,X4,Y4,'-o','Color',blue,'LineWidth',1.0,'MarkerSize',3.1, ...
        'MarkerFaceColor','w','DisplayName','$\Delta a=4$ mm');
    plot(ax1,X2,Y2,'-s','Color',orange,'LineWidth',1.0,'MarkerSize',2.6, ...
        'MarkerIndices',1:2:numel(X2),'MarkerFaceColor','w', ...
        'DisplayName','$\Delta a=2$ mm');
    xlabel(ax1,'$x$ [mm]','Interpreter','latex');
    ylabel(ax1,'$y$ [mm]','Interpreter','latex');
    axis(ax1,'equal'); box(ax1,'on');
    legend(ax1,'Location','northwest','NumColumns',2,'Box','off', ...
        'Interpreter','latex');

    ax2=nexttile(tl); hold(ax2,'on');
    yline(ax2,0,'-','Color',[.7 .7 .7],'LineWidth',.6);
    plot(ax2,common,dy_um,'-o','Color',orange,'LineWidth',1.0,'MarkerSize',3);
    xlabel(ax2,'Crack length $a$ [mm]','Interpreter','latex');
    ylabel(ax2,'$y_{2\mathrm{mm}}-y_{4\mathrm{mm}}\;[\mu\mathrm{m}]$', ...
        'Interpreter','latex'); box(ax2,'on');

    ax3=nexttile(tl); hold(ax3,'on');
    yline(ax3,0,'-','Color',[.7 .7 .7],'LineWidth',.6);
    plot(ax3,common,dtheta_mdeg,'-o','Color',orange,'LineWidth',1.0,'MarkerSize',3);
    xlabel(ax3,'Crack length $a$ [mm]','Interpreter','latex');
    ylabel(ax3,'$\theta_{2\mathrm{mm}}-\theta_{4\mathrm{mm}}$ [mdeg]', ...
        'Interpreter','latex'); box(ax3,'on');

    ax4=nexttile(tl); hold(ax4,'on');
    yline(ax4,0,'-','Color',[.7 .7 .7],'LineWidth',.6,'HandleVisibility','off');
    plot(ax4,T4c.crack_length_mm,T4c.KII_over_KI,'-o','Color',blue, ...
        'LineWidth',1.0,'MarkerSize',3.0,'DisplayName','$\Delta a=4$ mm');
    plot(ax4,T2c.crack_length_mm,T2c.KII_over_KI,'-s','Color',orange, ...
        'LineWidth',1.0,'MarkerSize',2.4,'DisplayName','$\Delta a=2$ mm');
    xlabel(ax4,'Crack length $a$ [mm]','Interpreter','latex');
    ylabel(ax4,'$K_{II}/K_I$','Interpreter','latex');
    box(ax4,'on');
    legend(ax4,'Location','southwest','Box','off','Interpreter','latex');

    allAx=[ax1 ax2 ax3 ax4];
    for ax=allAx
        set(ax,'FontName','Times New Roman','FontSize',8.5, ...
            'LineWidth',.75,'TickDir','out');
        grid(ax,'off');
    end

    Metrics=table(common, ...
        1e3*A4.tip_x_m,1e3*A4.tip_y_m, ...
        1e3*A2.tip_x_m,1e3*A2.tip_y_m, ...
        dx_um,dy_um,dr_um,dtheta_mdeg, ...
        A4.KII_over_KI,A2.KII_over_KI,dq, ...
        'VariableNames',{'crack_length_mm', ...
        'x_4mm_mm','y_4mm_mm','x_2mm_mm','y_2mm_mm', ...
        'dx_um','dy_um','tip_separation_um','dtheta_mdeg', ...
        'q_4mm','q_2mm','dq'});

    S=struct();
    S.run4Dir=run4;
    S.run2Dir=run2;
    S.maxCommonCrackLength_mm=maxCommon;
    S.metrics=Metrics;
    S.maxAbsDy_um=max(abs(dy_um));
    S.maxTipSeparation_um=max(dr_um);
    S.maxAbsDtheta_mdeg=max(abs(dtheta_mdeg));
    S.maxAbsDq=max(abs(dq));
    S.qMaximum4mm=struct('a_mm',aQ4,'q',qmax4);
    S.qMaximum2mm=struct('a_mm',aQ2,'q',qmax2);
    S.crossing4mm=struct('exists',cross4,'aLS_mm',aLS4);
    S.crossing2mm=struct('exists',cross2,'aLS_mm',aLS2);
    if cross4&&cross2
        S.localSymmetryDifference_mm=aLS2-aLS4;
    else
        S.localSymmetryDifference_mm=NaN;
    end

    fprintf('\n============================================================\n');
    fprintf('CRACK-INCREMENT SENSITIVITY: 4 mm vs 2 mm\n');
    fprintf('============================================================\n');
    fprintf('  common crack-length range   : %.6f--%.6f mm\n', ...
        min(common),max(common));
    fprintf('  max |dy|                    : %.9g micrometers\n',S.maxAbsDy_um);
    fprintf('  max tip separation          : %.9g micrometers\n',S.maxTipSeparation_um);
    fprintf('  max |dtheta|                : %.9g millidegrees\n',S.maxAbsDtheta_mdeg);
    fprintf('  max |dq| at matched lengths : %.9g\n',S.maxAbsDq);
    fprintf('  q maximum 4 mm              : a=%.6f mm, q=%+.9g\n',aQ4,qmax4);
    fprintf('  q maximum 2 mm              : a=%.6f mm, q=%+.9g\n',aQ2,qmax2);
    if cross4
        fprintf('  local symmetry 4 mm         : %.9f mm\n',aLS4);
    end
    if cross2
        fprintf('  local symmetry 2 mm         : %.9f mm\n',aLS2);
    end
    if cross4&&cross2
        fprintf('  crossing shift 2mm-4mm      : %+.9f mm\n',aLS2-aLS4);
    end

    if opt.Export
        names={'trajectory','vertical_deviation','direction_deviation','mode_mixity'};
        axesList={ax1,ax2,ax3,ax4};
        panelPdf=cell(4,1); panelPng=cell(4,1);
        for jj=1:4
            panelPdf{jj}=fullfile(outDir,[names{jj} '.pdf']);
            panelPng{jj}=fullfile(outDir,[names{jj} '.png']);
            local_export_axis(axesList{jj},panelPdf{jj},panelPng{jj},names{jj});
        end
        csvFile=fullfile(outDir,'figure_increment_sensitivity_metrics.csv');
        writetable(Metrics,csvFile);
        S.panelPdf=panelPdf;
        S.panelPng=panelPng;
        S.metricsCsv=csvFile;

        if opt.ExportCombined
            combinedPdf=fullfile(outDir,'combined_preview.pdf');
            combinedPng=fullfile(outDir,'combined_preview.png');
            exportgraphics(fig,combinedPdf,'ContentType','vector');
            exportgraphics(fig,combinedPng,'Resolution',600);
            S.combinedPdf=combinedPdf;
            S.combinedPng=combinedPng;
        end
    end
end

function T=local_rows_at_lengths(Tin,a)
    idx=zeros(numel(a),1);
    for j=1:numel(a)
        [d,k]=min(abs(Tin.crack_length_mm-a(j)));
        assert(d<=1e-8,'incfig:MatchTolerance', ...
            'Could not match crack length %.12g mm.',a(j));
        idx(j)=k;
    end
    T=Tin(idx,:);
end

function [ok,aLS]=local_zero_crossing(T)
    a=T.crack_length_mm;
    q=T.KII_over_KI;
    ok=false; aLS=NaN;
    for j=1:numel(q)-1
        if q(j)==0
            ok=true; aLS=a(j); return
        end
        if q(j)*q(j+1)<0
            aLS=a(j)-q(j)*(a(j+1)-a(j))/(q(j+1)-q(j));
            ok=true; return
        end
    end
end

function local_export_axis(srcAx,pdfFile,pngFile,panelName)
    if strcmp(panelName,'trajectory')
        sz=[17.0 6.7]; pos=[0.09 0.17 0.87 0.77];
    else
        sz=[6.6 5.2]; pos=[0.18 0.19 0.77 0.75];
    end
    f=figure('Color','w','Visible','off','Units','centimeters', ...
        'Position',[2 2 sz(1) sz(2)]);
    c=onCleanup(@()local_close_if_valid(f)); %#ok<NASGU>
    ax=copyobj(srcAx,f);
    set(ax,'Units','normalized','Position',pos);
    if any(strcmp(panelName,{'trajectory','mode_mixity'}))
        if strcmp(panelName,'trajectory'),loc='northwest';else,loc='southwest';end
        legend(ax,'show','Location',loc,'Box','off','Interpreter','latex');
    end
    drawnow;
    exportgraphics(ax,pdfFile,'ContentType','vector');
    exportgraphics(ax,pngFile,'Resolution',600);
end

function local_close_if_valid(h)
    if ~isempty(h)&&isgraphics(h),close(h);end
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
