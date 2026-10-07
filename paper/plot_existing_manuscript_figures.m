function S = plot_existing_manuscript_figures(varargin)
%PLOT_EXISTING_MANUSCRIPT_FIGURES Redraw the existing manuscript figures in MATLAB.
%
% This is a plot-only publication step. It reads the audited manuscript
% evidence in paper/data/evidence_exact.mat and performs no FE solve, mesh
% generation, EDI replay, COD refit, or MTS recomputation.
%
% Each panel is exported separately as vector PDF (plus a 600-dpi PNG
% preview). Panel letters and panel titles are intentionally omitted from
% MATLAB; LaTeX supplies them through subcaptions in paper/figures/*/figure.tex.
%
% Typical use from the repository root:
%   addpath(genpath(pwd));
%   S = plot_existing_manuscript_figures();
%
% Optional:
%   S = plot_existing_manuscript_figures('ShowFigures',true);
%
% Output:
%   paper/figures/trajectory/
%   paper/figures/intensities/
%   paper/figures/late/
%   paper/figures/cod/
%   paper/figures/quality/

    ip=inputParser;
    addParameter(ip,'EvidenceFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'OutputRoot','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'ShowFigures',false,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    paperDir=fileparts(mfilename('fullpath'));
    root=fileparts(paperDir);

    evidenceFile=char(opt.EvidenceFile);
    if isempty(evidenceFile)
        evidenceFile=fullfile(paperDir,'data','evidence_exact.mat');
    elseif ~local_is_absolute_path(evidenceFile)
        evidenceFile=fullfile(root,evidenceFile);
    end

    outRoot=char(opt.OutputRoot);
    if isempty(outRoot)
        outRoot=fullfile(paperDir,'figures');
    elseif ~local_is_absolute_path(outRoot)
        outRoot=fullfile(root,outRoot);
    end

    assert(exist(evidenceFile,'file')==2,'paperfig:MissingEvidence', ...
        'Audited evidence file not found: %s',evidenceFile);

    A=load(evidenceFile,'E','T','COD','Q');
    assert(isfield(A,'E')&&isfield(A,'T')&&isfield(A,'COD')&&isfield(A,'Q'), ...
        'paperfig:BadEvidence','Evidence MAT lacks E, T, COD, or Q.');

    E=A.E;
    T=A.T;
    COD=A.COD;
    Q=A.Q;

    assert(height(T)==23 && all(T.segment==(1:23)'), ...
        'paperfig:BadHistory','Expected accepted states P1--P23.');
    assert(height(COD)==176,'paperfig:BadCOD','Expected 176 stored COD fits.');
    assert(E.modeMixityPeakSegment==17,'paperfig:Peak','Expected q maximum at P17.');

    vis='off';
    if opt.ShowFigures,vis='on';end

    % Publication palette, matched to plot_m1_vs_reference_publication.m.
    C.blue=[0.0000 0.4470 0.7410];
    C.orange=[0.8500 0.3250 0.0980];
    C.green=[0.4660 0.6740 0.1880];
    C.purple=[0.4940 0.1840 0.5560];
    C.gray=[0.45 0.45 0.45];

    S=struct();
    S.evidenceFile=evidenceFile;
    S.outputRoot=outRoot;
    S.files={};

    %% Figure 1: trajectory and late detail
    out=fullfile(outRoot,'trajectory');
    local_mkdir(out);

    Vmm=1e3*E.vertices_m;
    accepted=Vmm(1:24,:);       % P0--P23
    extension=Vmm(24:25,:);     % qualified P23--P24; no P24 physical SIFs
    p17=Vmm(18,:);
    p23=Vmm(24,:);
    p24=Vmm(25,:);
    pLS=1e3*E.interpolatedZeroPoint_m;

    holeCenter=[170,-20];
    holeRadius=30;
    ang=linspace(0,2*pi,481);
    hx=holeCenter(1)+holeRadius*cos(ang);
    hy=holeCenter(2)+holeRadius*sin(ang);
    rightBoundary_mm=1e3*E.C.A;

    f=local_figure(vis,[17.0 6.7]); ax=axes(f); hold(ax,'on');
    plot(ax,hx,hy,'k-','LineWidth',1.0,'HandleVisibility','off');
    plot(ax,[rightBoundary_mm rightBoundary_mm],[-54 14],'k-','LineWidth',0.9,'HandleVisibility','off');
    plot(ax,accepted(:,1),accepted(:,2),'-o','Color',C.blue,'LineWidth',1.25, ...
        'MarkerSize',2.9,'MarkerFaceColor','w','DisplayName','Accepted P_0--P_{23}');
    plot(ax,extension(:,1),extension(:,2),'--x','Color',C.orange,'LineWidth',1.15, ...
        'MarkerSize',5,'DisplayName','Qualified, unsolved extension');
    plot(ax,p23(1),p23(2),'s','Color','k','MarkerFaceColor','k','MarkerSize',4, ...
        'HandleVisibility','off');
    axis(ax,'equal'); xlim(ax,[136 302]); ylim(ax,[-54 14]);
    xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');
    legend(ax,'Location','northeast','Box','off','Interpreter','tex');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'overview'); %#ok<AGROW>

    f=local_figure(vis,[17.0 4.5]); ax=axes(f); hold(ax,'on');
    plot(ax,accepted(16:end,1),accepted(16:end,2),'-o','Color',C.blue, ...
        'LineWidth',1.25,'MarkerSize',3.0,'MarkerFaceColor','w');
    plot(ax,extension(:,1),extension(:,2),'--x','Color',C.orange, ...
        'LineWidth',1.15,'MarkerSize',5);
    plot(ax,[rightBoundary_mm rightBoundary_mm],[-27.3 -21.3],'k-','LineWidth',0.9);
    plot(ax,p17(1),p17(2),'o','Color',C.green,'LineWidth',1.2,'MarkerSize',5);
    plot(ax,pLS(1),pLS(2),'d','Color',C.purple,'LineWidth',1.2,'MarkerSize',5);
    plot(ax,p23(1),p23(2),'s','Color','k','MarkerFaceColor','k','MarkerSize',4);
    text(ax,p17(1),p17(2),'  P_{17}','VerticalAlignment','bottom','FontSize',8);
    text(ax,pLS(1),pLS(2),'  q_K=0 estimate','VerticalAlignment','bottom','FontSize',8);
    text(ax,p23(1),p23(2),'  P_{23}','VerticalAlignment','top','HorizontalAlignment','right','FontSize',8);
    text(ax,p24(1),p24(2),'  P_{24}: unsolved','Color',C.orange, ...
        'VerticalAlignment','top','FontSize',8);
    axis(ax,'equal'); xlim(ax,[258 302]); ylim(ax,[-27.3 -21.3]);
    xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'late_detail'); %#ok<AGROW>

    %% Figure 2: KI, KII, absolute direction
    out=fullfile(outRoot,'intensities');
    local_mkdir(out);
    a=T.crack_length_mm;

    f=local_figure(vis,[14.0 4.6]); ax=axes(f);
    plot(ax,a,T.KI_unit,'-o','Color',C.blue,'LineWidth',1.2, ...
        'MarkerSize',3.1,'MarkerFaceColor','w');
    xlim(ax,[4 92]); xticks(ax,[4 20 36 52 68 84 92]);
    ylabel(ax,'K_I [MPa m^{1/2}]');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'KI'); %#ok<AGROW>

    f=local_figure(vis,[14.0 4.6]); ax=axes(f); hold(ax,'on');
    plot(ax,a,T.KII_unit,'-o','Color',C.blue,'LineWidth',1.2, ...
        'MarkerSize',3.1,'MarkerFaceColor','w');
    yline(ax,0,'--','Color',C.gray,'HandleVisibility','off');
    xlim(ax,[4 92]); xticks(ax,[4 20 36 52 68 84 92]);
    ylabel(ax,'K_{II} [MPa m^{1/2}]');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'KII'); %#ok<AGROW>

    f=local_figure(vis,[14.0 4.6]); ax=axes(f);
    plot(ax,a,T.theta_deg,'-o','Color',C.blue,'LineWidth',1.2, ...
        'MarkerSize',3.1,'MarkerFaceColor','w');
    xlim(ax,[4 92]); xticks(ax,[4 20 36 52 68 84 92]);
    xlabel(ax,'Crack length a [mm]'); ylabel(ax,'\theta_k [deg]');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'theta'); %#ok<AGROW>

    %% Figure 3: late mode mixity and MTS turn
    out=fullfile(outRoot,'late');
    local_mkdir(out);
    late=T.crack_length_mm>=60;
    aLate=T.crack_length_mm(late);
    aLS=E.interpolatedZeroLength_mm;

    f=local_figure(vis,[14.0 4.6]); ax=axes(f); hold(ax,'on');
    plot(ax,aLate,1e3*T.KII_over_KI(late),'-o','Color',C.blue,'LineWidth',1.2, ...
        'MarkerSize',3.1,'MarkerFaceColor','w');
    yline(ax,0,'--','Color',C.gray,'HandleVisibility','off');
    xline(ax,aLS,':','Color',C.purple,'LineWidth',1.0,'HandleVisibility','off');
    xlim(ax,[60 92]); xticks(ax,[60 68 76 84 88 92]);
    ylabel(ax,'10^3 K_{II}/K_I');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'mode_mixity'); %#ok<AGROW>

    f=local_figure(vis,[14.0 4.6]); ax=axes(f); hold(ax,'on');
    plot(ax,aLate,T.delta_theta_next_deg(late),'-o','Color',C.blue,'LineWidth',1.2, ...
        'MarkerSize',3.1,'MarkerFaceColor','w');
    yline(ax,0,'--','Color',C.gray,'HandleVisibility','off');
    xline(ax,aLS,':','Color',C.purple,'LineWidth',1.0,'HandleVisibility','off');
    xlim(ax,[60 92]); xticks(ax,[60 68 76 84 88 92]);
    xlabel(ax,'Crack length a [mm]'); ylabel(ax,'\Delta\theta_{k+1} [deg]');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'turn'); %#ok<AGROW>

    %% Figure 4: COD--EDI sensitivity, all four windows and both degrees
    out=fullfile(outRoot,'cod');
    local_mkdir(out);
    codLate=COD.crack_length_mm>=60;
    D=COD(codLate,:);
    windows=[0.04 0.20;0.04 0.30;0.08 0.30;0.12 0.30];
    labels={'[0.04,0.20] \Delta a','[0.04,0.30] \Delta a', ...
            '[0.08,0.30] \Delta a','[0.12,0.30] \Delta a'};
    cols={C.blue,C.orange,C.green,C.purple};
    marks={'o','s','^','d'};

    S.files{end+1}=local_plot_cod_panel(D,1,'turn_error_deg',1, ...
        'COD-EDI turn [deg]','turn_linear',out,vis,windows,labels,cols,marks,false); %#ok<AGROW>
    S.files{end+1}=local_plot_cod_panel(D,2,'turn_error_deg',1, ...
        'COD-EDI turn [deg]','turn_quadratic',out,vis,windows,labels,cols,marks,false); %#ok<AGROW>
    S.files{end+1}=local_plot_cod_panel(D,1,'ratio_error',1e5, ...
        '10^5 (COD-EDI mode mixity)','mixity_linear',out,vis,windows,labels,cols,marks,false); %#ok<AGROW>
    S.files{end+1}=local_plot_cod_panel(D,2,'ratio_error',1e5, ...
        '10^5 (COD-EDI mode mixity)','mixity_quadratic',out,vis,windows,labels,cols,marks,true); %#ok<AGROW>

    %% Figure S1: numerical quality
    out=fullfile(outRoot,'quality');
    local_mkdir(out);
    Qp=Q(Q.segment<=23,:);
    Tphys=T(T.segment>=2,:);

    f=local_figure(vis,[7.4 5.0]); ax=axes(f); hold(ax,'on');
    plot(ax,1e3*Qp.path_length_m,Qp.min_angle_deg,'-o','Color',C.blue, ...
        'LineWidth',1.15,'MarkerSize',3,'MarkerFaceColor','w');
    yline(ax,20,'--','Color',C.orange,'HandleVisibility','off');
    xlim(ax,[8 92]); xticks(ax,[8 36 64 92]); ylim(ax,[19 27]);
    ylabel(ax,'Minimum angle [deg]');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'min_angle'); %#ok<AGROW>

    f=local_figure(vis,[7.4 5.0]); ax=axes(f); hold(ax,'on');
    plot(ax,1e3*Qp.path_length_m,Qp.max_neighbor_ratio,'-o','Color',C.blue, ...
        'LineWidth',1.15,'MarkerSize',3,'MarkerFaceColor','w');
    yline(ax,1.8,'--','Color',C.orange,'HandleVisibility','off');
    xlim(ax,[8 92]); xticks(ax,[8 36 64 92]); ylim(ax,[1.78 1.804]);
    ylabel(ax,'Adjacent size ratio');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'size_ratio'); %#ok<AGROW>

    f=local_figure(vis,[7.4 5.0]); ax=axes(f);
    plot(ax,Tphys.crack_length_mm,Tphys.PCG_iterations,'-o','Color',C.blue, ...
        'LineWidth',1.15,'MarkerSize',3,'MarkerFaceColor','w');
    xlim(ax,[8 92]); xticks(ax,[8 36 64 92]);
    xlabel(ax,'Crack length a [mm]'); ylabel(ax,'PCG iterations');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'iterations'); %#ok<AGROW>

    f=local_figure(vis,[7.4 5.0]); ax=axes(f); hold(ax,'on');
    plot(ax,Tphys.crack_length_mm,Tphys.true_rel_residual/(5e-10),'-o', ...
        'Color',C.blue,'LineWidth',1.15,'MarkerSize',3,'MarkerFaceColor','w');
    yline(ax,1,'--','Color',C.orange,'HandleVisibility','off');
    xlim(ax,[8 92]); xticks(ax,[8 36 64 92]); ylim(ax,[0 1.1]);
    xlabel(ax,'Crack length a [mm]'); ylabel(ax,'True residual / gate');
    local_style_axis(ax);
    S.files{end+1}=local_export_axis(ax,out,'true_residual'); %#ok<AGROW>

    fprintf('\n============================================================\n');
    fprintf('MANUSCRIPT MATLAB FIGURE EXPORT COMPLETE\n');
    fprintf('============================================================\n');
    fprintf('  evidence : %s\n',evidenceFile);
    fprintf('  output   : %s\n',outRoot);
    fprintf('  panels   : %d vector PDFs (+ PNG previews)\n',numel(S.files));
    fprintf('  no FE solve / EDI replay / COD refit was performed.\n');
end

function out=local_plot_cod_panel(D,degree,varName,scale,yLab,fileName,outDir,vis, ...
        windows,labels,cols,marks,showLegend)
    f=local_figure(vis,[7.5 5.2]); ax=axes(f); hold(ax,'on');
    for j=1:size(windows,1)
        use=D.degree==degree & ...
            abs(D.lower_r_over_DeltaA-windows(j,1))<1e-12 & ...
            abs(D.upper_r_over_DeltaA-windows(j,2))<1e-12;
        X=D.crack_length_mm(use);
        Y=scale*D.(varName)(use);
        [X,ord]=sort(X);Y=Y(ord);
        plot(ax,X,Y,['-' marks{j}],'Color',cols{j},'LineWidth',1.05, ...
            'MarkerSize',3.0,'MarkerFaceColor','w','DisplayName',labels{j});
    end
    yline(ax,0,'--','Color',[0.45 0.45 0.45],'HandleVisibility','off');
    xlim(ax,[60 92]); xticks(ax,[60 76 92]);
    xlabel(ax,'Crack length a [mm]'); ylabel(ax,yLab);
    if showLegend
        legend(ax,'Location','best','Box','off','FontSize',7.2,'Interpreter','tex');
    end
    local_style_axis(ax);
    out=local_export_axis(ax,outDir,fileName);
end

function f=local_figure(vis,szcm)
    f=figure('Color','w','Visible',vis,'Units','centimeters', ...
        'Position',[2 2 szcm(1) szcm(2)]);
end

function local_style_axis(ax)
    set(ax,'FontName','Times New Roman','FontSize',8.5, ...
        'LineWidth',0.75,'TickDir','out','Box','on');
    grid(ax,'off');
end

function files=local_export_axis(ax,outDir,name)
    pdfFile=fullfile(outDir,[name '.pdf']);
    pngFile=fullfile(outDir,[name '.png']);
    exportgraphics(ax,pdfFile,'ContentType','vector');
    exportgraphics(ax,pngFile,'Resolution',600);
    files={pdfFile,pngFile};
end

function local_mkdir(p)
    if exist(p,'dir')~=7,mkdir(p);end
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
