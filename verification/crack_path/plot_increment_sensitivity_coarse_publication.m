function S = plot_increment_sensitivity_coarse_publication(varargin)
%PLOT_INCREMENT_SENSITIVITY_COARSE_PUBLICATION
% Three-level crack-increment comparison for Delta a = 4, 2, and 1 mm.
%
% All trajectories must come from main_increment_sensitivity_coarse and use
% the same dimensionless M1-exterior + CoreScale=2 numerical family.
%
% Publication panels:
%   (a) complete independently propagated trajectories;
%   (b) successive-level vertical-coordinate differences at exact common
%       4-mm crack lengths: (2 mm - 4 mm) and (1 mm - 2 mm);
%   (c) successive-level crack-direction differences at the same exact
%       common lengths;
%   (d) native turning density, (Delta theta_next)/(Delta a), versus crack
%       length for all three increments.
%
% Raw q_K = KII/KI is retained in the returned diagnostics and exported
% metrics, including its positive maximum and linearly interpolated zero
% crossing for each increment.  No path interpolation is used in panels
% (b,c).

    ip=inputParser;
    addParameter(ip,'Run4Dir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Run2Dir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Run1Dir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Export',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ExportCombined',false,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));

    run4=local_resolve_run_dir(root,char(opt.Run4Dir),'da_4mm');
    run2=local_resolve_run_dir(root,char(opt.Run2Dir),'da_2mm');
    run1=local_resolve_run_dir(root,char(opt.Run1Dir),'da_1mm');

    outDir=char(opt.OutputDir);
    if isempty(outDir)
        outDir=fullfile(root,'paper','figures','increment_sensitivity');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    [T4,V4]=local_load_run(run4,4);
    [T2,V2]=local_load_run(run2,2);
    [T1,V1]=local_load_run(run1,1);

    maxCommon=min([max(T4.crack_length_mm), ...
                   max(T2.crack_length_mm), ...
                   max(T1.crack_length_mm)]);
    T4=T4(T4.crack_length_mm<=maxCommon+1e-9,:);
    T2=T2(T2.crack_length_mm<=maxCommon+1e-9,:);
    T1=T1(T1.crack_length_mm<=maxCommon+1e-9,:);
    V4=V4(V4.crack_length_mm<=maxCommon+1e-9,:);
    V2=V2(V2.crack_length_mm<=maxCommon+1e-9,:);
    V1=V1(V1.crack_length_mm<=maxCommon+1e-9,:);

    % Exact common lengths.  The 4-mm state grid is a subset of both finer
    % grids when all runs reach the same physical length.
    common=T4.crack_length_mm;
    common=common(ismembertol(common,T2.crack_length_mm,1e-10,'DataScale',1));
    common=common(ismembertol(common,T1.crack_length_mm,1e-10,'DataScale',1));
    assert(~isempty(common),'incfig:NoCommonLengths','No exact common crack lengths.');

    A4=local_rows_at_lengths(T4,common);
    A2=local_rows_at_lengths(T2,common);
    A1=local_rows_at_lengths(T1,common);

    % Successive-level trajectory differences.
    dx24_um=1e6*(A2.tip_x_m-A4.tip_x_m);
    dy24_um=1e6*(A2.tip_y_m-A4.tip_y_m);
    dr24_um=hypot(dx24_um,dy24_um);
    dtheta24_mdeg=1e3*(A2.theta_deg-A4.theta_deg);

    dx12_um=1e6*(A1.tip_x_m-A2.tip_x_m);
    dy12_um=1e6*(A1.tip_y_m-A2.tip_y_m);
    dr12_um=hypot(dx12_um,dy12_um);
    dtheta12_mdeg=1e3*(A1.theta_deg-A2.theta_deg);

    % Raw mode mixity at exact common lengths: diagnostic only.
    dq24=A2.KII_over_KI-A4.KII_over_KI;
    dq12=A1.KII_over_KI-A2.KII_over_KI;

    % Native turning-density histories.  delta_theta_next_deg is the MTS
    % correction applied to generate the following leg.
    kappa4=T4.delta_theta_next_deg/4;
    kappa2=T2.delta_theta_next_deg/2;
    kappa1=T1.delta_theta_next_deg/1;

    % Native q diagnostics.
    D4=local_q_diagnostics(T4,4);
    D2=local_q_diagnostics(T2,2);
    D1=local_q_diagnostics(T1,1);

    % Actual trajectories.
    X4=1e3*V4.x_m; Y4=1e3*V4.y_m;
    X2=1e3*V2.x_m; Y2=1e3*V2.y_m;
    X1=1e3*V1.x_m; Y1=1e3*V1.y_m;

    holeCenter=[170,-20];
    holeRadius=30;
    ang=linspace(0,2*pi,361);
    hx=holeCenter(1)+holeRadius*cos(ang);
    hy=holeCenter(2)+holeRadius*sin(ang);

    blue=[0.0000 0.4470 0.7410];
    orange=[0.8500 0.3250 0.0980];
    yellow=[0.9290 0.6940 0.1250];
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
    plot(ax1,X1,Y1,'-^','Color',yellow,'LineWidth',1.0,'MarkerSize',2.5, ...
        'MarkerIndices',1:4:numel(X1),'MarkerFaceColor','w', ...
        'DisplayName','$\Delta a=1$ mm');
    xlabel(ax1,'$x$ [mm]','Interpreter','latex');
    ylabel(ax1,'$y$ [mm]','Interpreter','latex');
    axis(ax1,'equal'); box(ax1,'on');
    legend(ax1,'Location','northwest','NumColumns',3,'Box','off', ...
        'Interpreter','latex');

    ax2=nexttile(tl); hold(ax2,'on');
    yline(ax2,0,'-','Color',[.7 .7 .7],'LineWidth',.6,'HandleVisibility','off');
    plot(ax2,common,dy24_um,'-o','Color',orange,'LineWidth',1.0, ...
        'MarkerSize',3,'DisplayName','$2-4$ mm');
    plot(ax2,common,dy12_um,'-^','Color',yellow,'LineWidth',1.0, ...
        'MarkerSize',3,'DisplayName','$1-2$ mm');
    xlabel(ax2,'Crack length $a$ [mm]','Interpreter','latex');
    ylabel(ax2,'Successive $\Delta y$ [$\mu$m]','Interpreter','latex');
    legend(ax2,'Location','best','Box','off','Interpreter','latex');
    box(ax2,'on');

    ax3=nexttile(tl); hold(ax3,'on');
    yline(ax3,0,'-','Color',[.7 .7 .7],'LineWidth',.6,'HandleVisibility','off');
    plot(ax3,common,dtheta24_mdeg,'-o','Color',orange,'LineWidth',1.0, ...
        'MarkerSize',3,'DisplayName','$2-4$ mm');
    plot(ax3,common,dtheta12_mdeg,'-^','Color',yellow,'LineWidth',1.0, ...
        'MarkerSize',3,'DisplayName','$1-2$ mm');
    xlabel(ax3,'Crack length $a$ [mm]','Interpreter','latex');
    ylabel(ax3,'Successive $\Delta\theta$ [mdeg]','Interpreter','latex');
    legend(ax3,'Location','best','Box','off','Interpreter','latex');
    box(ax3,'on');

    ax4=nexttile(tl); hold(ax4,'on');
    yline(ax4,0,'-','Color',[.7 .7 .7],'LineWidth',.6,'HandleVisibility','off');
    plot(ax4,T4.crack_length_mm,kappa4,'-o','Color',blue, ...
        'LineWidth',1.0,'MarkerSize',3.0,'DisplayName','$\Delta a=4$ mm');
    plot(ax4,T2.crack_length_mm,kappa2,'-s','Color',orange, ...
        'LineWidth',1.0,'MarkerSize',2.4,'DisplayName','$\Delta a=2$ mm');
    plot(ax4,T1.crack_length_mm,kappa1,'-^','Color',yellow, ...
        'LineWidth',1.0,'MarkerSize',2.2,'DisplayName','$\Delta a=1$ mm');
    xlabel(ax4,'Crack length $a$ [mm]','Interpreter','latex');
    ylabel(ax4,'$(\Delta\theta/\Delta a)$ [deg/mm]','Interpreter','latex');
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
        1e3*A1.tip_x_m,1e3*A1.tip_y_m, ...
        dx24_um,dy24_um,dr24_um,dtheta24_mdeg, ...
        dx12_um,dy12_um,dr12_um,dtheta12_mdeg, ...
        A4.KII_over_KI,A2.KII_over_KI,A1.KII_over_KI,dq24,dq12, ...
        'VariableNames',{'crack_length_mm', ...
        'x_4mm_mm','y_4mm_mm','x_2mm_mm','y_2mm_mm','x_1mm_mm','y_1mm_mm', ...
        'dx_2minus4_um','dy_2minus4_um','dr_2minus4_um','dtheta_2minus4_mdeg', ...
        'dx_1minus2_um','dy_1minus2_um','dr_1minus2_um','dtheta_1minus2_mdeg', ...
        'q_4mm','q_2mm','q_1mm','dq_2minus4','dq_1minus2'});

    Native=table();
    Native.increment_mm=[4*ones(height(T4),1);2*ones(height(T2),1);ones(height(T1),1)];
    Native.crack_length_mm=[T4.crack_length_mm;T2.crack_length_mm;T1.crack_length_mm];
    Native.q=[T4.KII_over_KI;T2.KII_over_KI;T1.KII_over_KI];
    Native.q_over_da=[T4.KII_over_KI/4;T2.KII_over_KI/2;T1.KII_over_KI];
    Native.turn_density_deg_per_mm=[kappa4;kappa2;kappa1];

    S=struct();
    S.run4Dir=run4; S.run2Dir=run2; S.run1Dir=run1;
    S.maxCommonCrackLength_mm=maxCommon;
    S.metrics=Metrics;
    S.native=Native;

    S.twoMinusFour=local_pair_metrics(dy24_um,dr24_um,dtheta24_mdeg,dq24);
    S.oneMinusTwo=local_pair_metrics(dy12_um,dr12_um,dtheta12_mdeg,dq12);

    S.level4mm=D4;
    S.level2mm=D2;
    S.level1mm=D1;

    % Simple successive-difference ratios are diagnostics only; three
    % increments do not by themselves prove asymptotic convergence.
    S.successiveRatioMaxSeparation= ...
        S.twoMinusFour.maxTipSeparation_um / max(S.oneMinusTwo.maxTipSeparation_um,eps);
    S.successiveRatioMaxDirection= ...
        S.twoMinusFour.maxAbsDtheta_mdeg / max(S.oneMinusTwo.maxAbsDtheta_mdeg,eps);

    fprintf('\n============================================================\n');
    fprintf('CRACK-INCREMENT SENSITIVITY: 4 mm / 2 mm / 1 mm\n');
    fprintf('============================================================\n');
    fprintf('  exact common range       : %.6f--%.6f mm\n',min(common),max(common));
    local_print_pair('2mm - 4mm',S.twoMinusFour);
    local_print_pair('1mm - 2mm',S.oneMinusTwo);
    local_print_level('4 mm',D4);
    local_print_level('2 mm',D2);
    local_print_level('1 mm',D1);
    fprintf('  ratio max separation (2-4)/(1-2) : %.6g\n', ...
        S.successiveRatioMaxSeparation);
    fprintf('  ratio max direction  (2-4)/(1-2) : %.6g\n', ...
        S.successiveRatioMaxDirection);

    if opt.Export
        names={'trajectory','vertical_deviation','direction_deviation','turn_density'};
        axesList={ax1,ax2,ax3,ax4};
        panelPdf=cell(4,1); panelPng=cell(4,1);
        for jj=1:4
            panelPdf{jj}=fullfile(outDir,[names{jj} '.pdf']);
            panelPng{jj}=fullfile(outDir,[names{jj} '.png']);
            local_export_axis(axesList{jj},panelPdf{jj},panelPng{jj},names{jj});
        end

        metricsCsv=fullfile(outDir,'figure_increment_sensitivity_metrics.csv');
        nativeCsv=fullfile(outDir,'figure_increment_sensitivity_native.csv');
        writetable(Metrics,metricsCsv);
        writetable(Native,nativeCsv);

        S.panelPdf=panelPdf;
        S.panelPng=panelPng;
        S.metricsCsv=metricsCsv;
        S.nativeCsv=nativeCsv;

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

function runDir=local_resolve_run_dir(root,userValue,defaultTag)
    runDir=userValue;
    if isempty(runDir)
        runDir=fullfile(root,'verification','crack_path', ...
            'increment_sensitivity_coarse',defaultTag);
    elseif ~local_is_absolute_path(runDir)
        runDir=fullfile(root,runDir);
    end
end

function [T,V]=local_load_run(runDir,da_mm)
    f=fullfile(runDir,'states.csv');
    vf=fullfile(runDir,'vertices.csv');
    assert(exist(f,'file')==2,'incfig:MissingStates', ...
        '%.6g-mm states not found: %s',da_mm,f);
    assert(exist(vf,'file')==2,'incfig:MissingVertices', ...
        '%.6g-mm vertices not found: %s',da_mm,vf);

    T=sortrows(readtable(f),'crack_length_mm');
    V=sortrows(readtable(vf),'crack_length_mm');

    req={'segment','tip_x_m','tip_y_m','theta_deg','KI_unit','KII_unit', ...
        'KII_over_KI','delta_theta_next_deg','theta_next_deg','crack_length_mm'};
    assert(all(ismember(req,T.Properties.VariableNames)), ...
        'incfig:Columns','%.6g-mm state table lacks required columns.',da_mm);
    assert(all(abs(T.crack_length_mm-da_mm*T.segment)<=1e-8), ...
        'incfig:IncrementMismatch','%.6g-mm state table has inconsistent lengths.',da_mm);
end

function P=local_pair_metrics(dy,dr,dtheta,dq)
    P=struct();
    P.maxAbsDy_um=max(abs(dy));
    P.maxTipSeparation_um=max(dr);
    P.maxAbsDtheta_mdeg=max(abs(dtheta));
    P.maxAbsDq=max(abs(dq));
end

function D=local_q_diagnostics(T,da_mm)
    [qmax,i]=max(T.KII_over_KI);
    [cross,aLS]=local_zero_crossing(T);
    D=struct();
    D.increment_mm=da_mm;
    D.qMaximum=qmax;
    D.qMaximumCrackLength_mm=T.crack_length_mm(i);
    D.qMaximumOverIncrement=qmax/da_mm;
    D.turnDensityAtQMaximum_deg_per_mm=T.delta_theta_next_deg(i)/da_mm;
    D.crossingExists=cross;
    D.localSymmetry_mm=aLS;
end

function local_print_pair(label,P)
    fprintf('\n  %s\n',label);
    fprintf('    max |dy|           : %.9g micrometers\n',P.maxAbsDy_um);
    fprintf('    max separation     : %.9g micrometers\n',P.maxTipSeparation_um);
    fprintf('    max |dtheta|       : %.9g millidegrees\n',P.maxAbsDtheta_mdeg);
    fprintf('    max |dq|           : %.9g\n',P.maxAbsDq);
end

function local_print_level(label,D)
    fprintf('\n  Delta a = %s\n',label);
    fprintf('    q maximum          : a=%.6f mm, q=%+.9g\n', ...
        D.qMaximumCrackLength_mm,D.qMaximum);
    fprintf('    q_max / Delta a    : %+.9g 1/mm\n',D.qMaximumOverIncrement);
    fprintf('    turn density there : %+.9g deg/mm\n', ...
        D.turnDensityAtQMaximum_deg_per_mm);
    if D.crossingExists
        fprintf('    local symmetry     : %.9f mm\n',D.localSymmetry_mm);
    else
        fprintf('    local symmetry     : no sign change found\n');
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
    if any(strcmp(panelName,{'trajectory','vertical_deviation', ...
            'direction_deviation','turn_density'}))
        legend(ax,'show','Location','best','Box','off','Interpreter','latex');
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
