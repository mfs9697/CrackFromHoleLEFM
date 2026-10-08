function S = plot_m1_vs_reference_publication(varargin)
%PLOT_M1_VS_REFERENCE_PUBLICATION
% Publication figure from the ACTUAL saved reference and independent-M1 runs.
%
% The figure contains:
%   (a) real reference and M1 crack trajectories;
%   (b) accumulated signed vertical tip deviation, M1-reference [micrometers];
%   (c) accumulated direction deviation, M1-reference [millidegrees];
%   (d) mode-mixity histories q=KII/KI, with the P21-P22 zero-crossing.
%
% No coordinates or SIF values are hard-coded.  All trajectory and
% mechanical data are read from path_run_state.mat files produced by the
% completed runs.
%
% Typical use from the repository root:
%   S = plot_m1_vs_reference_publication();
%
% Optional:
%   S = plot_m1_vs_reference_publication('Export',false);
%
% Output files (default):
%   paper/figures/m1_mesh_sensitivity/
%       trajectory.pdf / trajectory.png
%       vertical_deviation.pdf / vertical_deviation.png
%       direction_deviation.pdf / direction_deviation.png
%       mode_mixity.pdf / mode_mixity.png
%       figure_m1_vs_reference_publication_metrics.csv
%
% Plot titles and panel letters are intentionally NOT embedded in the
% graphics.  They belong to the LaTeX subcaptions in figure.tex.

    ip=inputParser;
    addParameter(ip,'ReferenceStateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'M1StateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Export',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ExportCombined',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'ShowMarkersEvery',2,@(x)isnumeric(x)&&isscalar(x)&& ...
        isfinite(x)&&x>=1&&x==round(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    here=fileparts(mfilename('fullpath'));
    root=fileparts(fileparts(here));

    refFile=char(opt.ReferenceStateFile);
    if isempty(refFile)
        refFile=fullfile(root,'verification','crack_path', ...
            'final_clean_run','path_run_state.mat');
    elseif ~local_is_absolute_path(refFile)
        refFile=fullfile(root,refFile);
    end

    m1File=char(opt.M1StateFile);
    if isempty(m1File)
        m1File=fullfile(root,'verification','crack_path', ...
            'm1_independent_run','trajectory','path_run_state.mat');
    elseif ~local_is_absolute_path(m1File)
        m1File=fullfile(root,m1File);
    end

    outDir=char(opt.OutputDir);
    if isempty(outDir)
        outDir=fullfile(root,'paper','figures','m1_mesh_sensitivity');
    elseif ~local_is_absolute_path(outDir)
        outDir=fullfile(root,outDir);
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    % The historical run-state MAT is convenient when present, but it is
    % not required for a portable paper checkout.  The committed paper/data
    % snapshot is the canonical fallback for the reference trajectory.
    referenceMatAvailable=(exist(refFile,'file')==2);
    assert(exist(m1File,'file')==2,'m1fig:MissingM1', ...
        'M1 state not found: %s',m1File);

    [Tr,VrAll,refCompleted,referenceSource]= ...
        local_load_reference(root,refFile,referenceMatAvailable);

    B=load(m1File,'State');
    assert(isfield(B,'State')&&isstruct(B.State),'m1fig:BadM1', ...
        'M1 MAT must contain State.');
    M=B.State;
    Tm=local_state_table(M);

    common=intersect(Tr.segment,Tm.segment,'stable');
    assert(~isempty(common),'m1fig:NoCommonStates','No common solved states.');
    Tr=Tr(ismember(Tr.segment,common),:);
    Tm=Tm(ismember(Tm.segment,common),:);
    [~,ir]=sort(Tr.segment); Tr=Tr(ir,:);
    [~,im]=sort(Tm.segment); Tm=Tm(im,:);
    assert(isequal(Tr.segment,Tm.segment),'m1fig:SegmentMismatch', ...
        'Reference and M1 common segment ordering differs.');

    kmax=min([refCompleted,M.completedPhysicalSegments,max(common)]);
    Vr=VrAll(1:kmax+1,:);
    Vm=M.vertices(1:kmax+1,:);
    assert(size(Vr,1)==size(Vm,1),'m1fig:VertexCount', ...
        'Reference and M1 accepted vertex counts differ through P%d.',kmax);

    % Crack length follows the imposed fixed increment. Infer it from the
    % first accepted segment rather than hard-coding 4 mm.
    da=norm(Vr(2,:)-Vr(1,:));
    a_mm=1e3*da*Tr.segment;

    dx_um=1e6*(Tm.tip_x_m-Tr.tip_x_m);
    dy_um=1e6*(Tm.tip_y_m-Tr.tip_y_m);
    dr_um=hypot(dx_um,dy_um);
    dtheta_mdeg=1e3*(Tm.theta_deg-Tr.theta_deg);

    qr=Tr.KII_over_KI;
    qm=Tm.KII_over_KI;

    % Event locations are derived from the actual histories.
    [~,iMaxR]=max(qr);
    [~,iMaxM]=max(qm);
    kMaxR=Tr.segment(iMaxR);
    kMaxM=Tm.segment(iMaxM);

    [crossR,aLSr_mm,xLSr_mm,yLSr_mm]=local_zero_crossing(Tr,Vr,da);
    [crossM,aLSm_mm,xLSm_mm,yLSm_mm]=local_zero_crossing(Tm,Vm,da);

    % Real accepted path coordinates in mm.
    Xr=1e3*Vr(:,1); Yr=1e3*Vr(:,2);
    Xm=1e3*Vm(:,1); Ym=1e3*Vm(:,2);

    % Fixed physical hole geometry used by these runs.
    holeCenter_mm=[170,-20];
    holeRadius_mm=30;
    ang=linspace(0,2*pi,361);
    hx=holeCenter_mm(1)+holeRadius_mm*cos(ang);
    hy=holeCenter_mm(2)+holeRadius_mm*sin(ang);

    % restrained publication palette
    blue=[0.0000 0.4470 0.7410];
    orange=[0.8500 0.3250 0.0980];
    green=[0.4660 0.6740 0.1880];
    purple=[0.4940 0.1840 0.5560];
    gray=[0.45 0.45 0.45];

    fig=figure('Color','w','Units','centimeters', ...
        'Position',[2 2 18.2 18.0]);
    tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

    % --------------------------------------------------------------
    % (a) Actual trajectory overlay: wide top panel
    % --------------------------------------------------------------
    ax1=nexttile(tl,[1 3]); hold(ax1,'on'); box(ax1,'on');
    plot(ax1,hx,hy,'k-','LineWidth',1.1,'HandleVisibility','off');

    mr=1:opt.ShowMarkersEvery:numel(Xr);
    mm=1:opt.ShowMarkersEvery:numel(Xm);
    pr=plot(ax1,Xr,Yr,'-','Color',blue,'LineWidth',1.35, ...
        'Marker','o','MarkerIndices',mr,'MarkerSize',3.3, ...
        'MarkerFaceColor','w','DisplayName','reference');
    pm=plot(ax1,Xm,Ym,'--','Color',orange,'LineWidth',1.15, ...
        'Marker','s','MarkerIndices',mm,'MarkerSize',3.0, ...
        'MarkerFaceColor','w','DisplayName','M1');

    % P17 / actual q maximum markers (computed, not assumed)
    scatter(ax1,1e3*Tr.tip_x_m(iMaxR),1e3*Tr.tip_y_m(iMaxR), ...
        45,'o','MarkerEdgeColor',purple,'LineWidth',1.3, ...
        'DisplayName',sprintf('reference q max: P%d',kMaxR));
    scatter(ax1,1e3*Tm.tip_x_m(iMaxM),1e3*Tm.tip_y_m(iMaxM), ...
        34,'o','MarkerEdgeColor',orange,'LineWidth',1.1, ...
        'HandleVisibility','off');

    if crossR
        scatter(ax1,xLSr_mm,yLSr_mm,45,'d','MarkerEdgeColor',green, ...
            'LineWidth',1.3,'DisplayName', ...
            sprintf('local symmetry: %.3f mm',aLSr_mm));
    end

    scatter(ax1,Xr(end),Yr(end),42,'s','MarkerEdgeColor',blue, ...
        'LineWidth',1.3,'DisplayName',sprintf('P%d: last accepted',kmax));

    xlabel(ax1,'$x$ [mm]','Interpreter','latex');
    ylabel(ax1,'$y$ [mm]','Interpreter','latex');
    % No title/panel label here: LaTeX supplies the subcaption.
    axis(ax1,'equal');
    xlim(ax1,[135 300]);
    ylim(ax1,[-55 15]);
    legend(ax1,'Location','northwest','NumColumns',2,'Box','off');

    % --------------------------------------------------------------
    % (b) accumulated geometric deviation
    % --------------------------------------------------------------
    ax2=nexttile(tl); hold(ax2,'on'); box(ax2,'on');
    plot(ax2,a_mm,dy_um,'-o','Color',blue,'LineWidth',1.15, ...
        'MarkerSize',3,'MarkerFaceColor','w');
    yline(ax2,0,':','Color',gray,'HandleVisibility','off');
    xlabel(ax2,'Crack length $a$ [mm]','Interpreter','latex');
    ylabel(ax2,'$y_{\mathrm{M1}}-y_{\mathrm{ref}}\;[\mu\mathrm{m}]$','Interpreter','latex');
    % No title/panel label here: LaTeX supplies the subcaption.

    % --------------------------------------------------------------
    % (c) accumulated angular deviation
    % --------------------------------------------------------------
    ax3=nexttile(tl); hold(ax3,'on'); box(ax3,'on');
    plot(ax3,a_mm,dtheta_mdeg,'-o','Color',purple,'LineWidth',1.15, ...
        'MarkerSize',3,'MarkerFaceColor','w');
    yline(ax3,0,':','Color',gray,'HandleVisibility','off');
    xlabel(ax3,'Crack length $a$ [mm]','Interpreter','latex');
    ylabel(ax3,'$\theta_{\mathrm{M1}}-\theta_{\mathrm{ref}}$ [mdeg]','Interpreter','latex');
    % No title/panel label here: LaTeX supplies the subcaption.

    % --------------------------------------------------------------
    % (d) actual mode-mixity histories
    % --------------------------------------------------------------
    ax4=nexttile(tl); hold(ax4,'on'); box(ax4,'on');
    plot(ax4,a_mm,qr,'-o','Color',blue,'LineWidth',1.15, ...
        'MarkerSize',3,'MarkerFaceColor','w','DisplayName','reference');
    plot(ax4,a_mm,qm,'--s','Color',orange,'LineWidth',1.05, ...
        'MarkerSize',3,'MarkerFaceColor','w','DisplayName','M1');
    yline(ax4,0,':','Color',gray,'HandleVisibility','off');
    if crossR
        xline(ax4,aLSr_mm,'--','Color',green,'LineWidth',0.9, ...
            'HandleVisibility','off');
    end
    xlabel(ax4,'Crack length $a$ [mm]','Interpreter','latex');
    ylabel(ax4,'$K_{II}/K_I$','Interpreter','latex');
    % No title/panel label here: LaTeX supplies the subcaption.
    legend(ax4,'Location','southwest','Box','off');

    % Consistent typography.
    axs=[ax1 ax2 ax3 ax4];
    for ax=axs
        set(ax,'FontName','Times New Roman','FontSize',8.5, ...
            'LineWidth',0.75,'TickDir','out');
        grid(ax,'off');
    end

    % Summary table saved with the figure.
    Metrics=table(Tr.segment,a_mm, ...
        1e3*Tr.tip_x_m,1e3*Tr.tip_y_m, ...
        1e3*Tm.tip_x_m,1e3*Tm.tip_y_m, ...
        dx_um,dy_um,dr_um,dtheta_mdeg,qr,qm,qm-qr, ...
        'VariableNames',{ ...
        'segment','crack_length_mm', ...
        'x_ref_mm','y_ref_mm','x_M1_mm','y_M1_mm', ...
        'dx_um','dy_um','tip_separation_um','dtheta_mdeg', ...
        'q_ref','q_M1','dq'});

    S=struct();
    S.referenceSource=referenceSource;
    S.referenceStateFile=refFile;
    S.M1StateFile=m1File;
    S.increment_mm=1e3*da;
    S.maxCommonSegment=kmax;
    S.metrics=Metrics;
    S.referenceQMaximumSegment=kMaxR;
    S.M1QMaximumSegment=kMaxM;
    S.referenceLocalSymmetry_mm=aLSr_mm;
    S.M1LocalSymmetry_mm=aLSm_mm;
    S.localSymmetryDifference_um=1e3*(aLSm_mm-aLSr_mm);
    S.maxAbsDy_um=max(abs(dy_um));
    S.maxTipSeparation_um=max(dr_um);
    S.maxAbsDtheta_mdeg=max(abs(dtheta_mdeg));
    S.maxAbsDq=max(abs(qm-qr));

    fprintf('\n============================================================\n');
    fprintf('REFERENCE vs M1 PUBLICATION FIGURE\n');
    fprintf('============================================================\n');
    fprintf('  common accepted states       : P1--P%d\n',kmax);
    fprintf('  inferred increment           : %.9f mm\n',S.increment_mm);
    fprintf('  reference source             : %s\n',referenceSource);
    fprintf('  q maximum reference / M1     : P%d / P%d\n',kMaxR,kMaxM);
    fprintf('  local symmetry reference     : %.9f mm\n',aLSr_mm);
    fprintf('  local symmetry M1            : %.9f mm\n',aLSm_mm);
    fprintf('  local-symmetry difference    : %+.6f micrometers\n', ...
        S.localSymmetryDifference_um);
    fprintf('  max |dy|                     : %.6g micrometers\n',S.maxAbsDy_um);
    fprintf('  max tip separation           : %.6g micrometers\n',S.maxTipSeparation_um);
    fprintf('  max |dtheta|                 : %.6g millidegrees\n',S.maxAbsDtheta_mdeg);
    fprintf('  max |dq|                     : %.6g\n',S.maxAbsDq);

    if opt.Export
        panelNames={'trajectory','vertical_deviation', ...
            'direction_deviation','mode_mixity'};
        panelAxes={ax1,ax2,ax3,ax4};
        panelPdf=cell(4,1);
        panelPng=cell(4,1);

        for jj=1:4
            panelPdf{jj}=fullfile(outDir,[panelNames{jj} '.pdf']);
            panelPng{jj}=fullfile(outDir,[panelNames{jj} '.png']);
            local_export_tiled_axis_copy( ...
                panelAxes{jj},panelPdf{jj},panelPng{jj},panelNames{jj});
        end

        csvFile=fullfile(outDir,'figure_m1_vs_reference_publication_metrics.csv');
        writetable(Metrics,csvFile);

        S.panelPdf=panelPdf;
        S.panelPng=panelPng;
        S.metricsCsv=csvFile;

        fprintf('  separate vector panels       : %s\n',outDir);
        for jj=1:4
            fprintf('    %-20s : %s\n',panelNames{jj},panelPdf{jj});
        end
        fprintf('  plotted-data CSV             : %s\n',csvFile);

        if opt.ExportCombined
            combinedPdf=fullfile(outDir,'combined_preview.pdf');
            combinedPng=fullfile(outDir,'combined_preview.png');
            exportgraphics(fig,combinedPdf,'ContentType','vector');
            exportgraphics(fig,combinedPng,'Resolution',600);
            S.combinedPdf=combinedPdf;
            S.combinedPng=combinedPng;
            fprintf('  combined preview             : %s\n',combinedPdf);
        end
    end
end

function local_export_tiled_axis_copy(srcAx,pdfFile,pngFile,panelName)
% Export one tiled-layout panel through an independent standalone figure.
% Some MATLAB releases invalidate a shared tiled-layout figure when
% exportgraphics is called directly on a child axes.  Cloning the axes
% isolates the publication export from the combined preview and matches the
% standalone-axis export path used by plot_existing_manuscript_figures.

    if strcmp(panelName,'trajectory')
        szcm=[17.0 6.7];
        pos=[0.09 0.17 0.87 0.77];
    else
        szcm=[6.6 5.2];
        pos=[0.18 0.19 0.77 0.75];
    end

    tmpFig=figure('Color','w','Visible','off','Units','centimeters', ...
        'Position',[2 2 szcm(1) szcm(2)]);
    cleanup=onCleanup(@()local_close_if_valid(tmpFig)); %#ok<NASGU>

    tmpAx=copyobj(srcAx,tmpFig);
    set(tmpAx,'Units','normalized','Position',pos);

    % Legends are figure/tiled-layout illustration objects rather than axes
    % children on some MATLAB releases, so recreate them from DisplayName.
    switch panelName
        case 'trajectory'
            legend(tmpAx,'show','Location','northwest', ...
                'NumColumns',2,'Box','off');
        case 'mode_mixity'
            legend(tmpAx,'show','Location','southwest','Box','off');
    end

    drawnow;
    exportgraphics(tmpAx,pdfFile,'ContentType','vector');
    exportgraphics(tmpAx,pngFile,'Resolution',600);
end

function local_close_if_valid(h)
    if ~isempty(h) && isgraphics(h)
        close(h);
    end
end


function [T,V,completed,source]=local_load_reference(root,refFile,matAvailable)
% Load the accepted h0 reference from the historical run-state MAT when
% available, otherwise from the committed publication evidence snapshot.

    if matAvailable
        A=load(refFile,'State');
        assert(isfield(A,'State')&&isstruct(A.State), ...
            'm1fig:BadReference','Reference MAT must contain State.');
        R=A.State;
        T=local_state_table(R);
        assert(isfield(R,'vertices')&&size(R.vertices,2)==2, ...
            'm1fig:ReferenceVertices','Reference State lacks vertices.');
        V=R.vertices;
        if isfield(R,'completedPhysicalSegments')
            completed=R.completedPhysicalSegments;
        else
            completed=height(T);
        end
        source=refFile;
        return
    end

    csvFile=fullfile(root,'paper','data','accepted_states.csv');
    jsonFile=fullfile(root,'paper','data','evidence.json');
    assert(exist(csvFile,'file')==2,'m1fig:MissingReferenceCSV', ...
        'Committed reference table not found: %s',csvFile);
    assert(exist(jsonFile,'file')==2,'m1fig:MissingReferenceJSON', ...
        'Committed reference evidence not found: %s',jsonFile);

    T=readtable(csvFile);
    E=jsondecode(fileread(jsonFile));
    assert(isfield(E,'vertices_m')&&size(E.vertices_m,2)==2, ...
        'm1fig:ReferenceEvidenceVertices', ...
        'paper/data/evidence.json lacks vertices_m.');
    V=double(E.vertices_m);
    completed=height(T);
    assert(completed>=1 && size(V,1)>=completed+1, ...
        'm1fig:ReferenceEvidenceLength', ...
        'Committed reference table/vertices are inconsistent.');
    source='paper/data/accepted_states.csv + evidence.json';
end

function T=local_state_table(S)
    assert(isfield(S,'rowsThroughCompleted')&&isfield(S,'rowVariableNames'), ...
        'm1fig:BadState','State lacks embedded physical-history rows.');
    names=cellstr(string(S.rowVariableNames));
    T=array2table(S.rowsThroughCompleted,'VariableNames',names);
    if isfield(S,'completedPhysicalSegments')
        T=T(1:S.completedPhysicalSegments,:);
    end
end

function [ok,aLS_mm,xLS_mm,yLS_mm]=local_zero_crossing(T,V,da)
    q=T.KII_over_KI;
    ok=false;
    aLS_mm=NaN; xLS_mm=NaN; yLS_mm=NaN;
    j=find(q(1:end-1).*q(2:end)<=0,1,'first');
    if isempty(j),return,end

    q0=q(j); q1=q(j+1);
    if q1==q0
        alpha=0.5;
    else
        alpha=-q0/(q1-q0);
    end
    alpha=max(0,min(1,alpha));

    k0=T.segment(j);
    k1=T.segment(j+1);
    a0=k0*da;
    a1=k1*da;
    aLS_mm=1e3*((1-alpha)*a0+alpha*a1);

    % Vertex row 1 is P0, so physical segment Pk is vertex k+1.
    p0=V(k0+1,:);
    p1=V(k1+1,:);
    p=(1-alpha)*p0+alpha*p1;
    xLS_mm=1e3*p(1);
    yLS_mm=1e3*p(2);
    ok=true;
end

function tf=local_is_absolute_path(p)
    p=char(p);
    tf=startsWith(p,filesep)|| ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))||startsWith(p,'\\');
end
