function Out = plot_accepted_clean_run_snapshot(varargin)
%PLOT_ACCEPTED_CLEAN_RUN_SNAPSHOT Plot the accepted P1-P23 scientific record.
%
% This lightweight fallback is intended for a cloned repository that does
% not contain the local final_clean_run MAT artifacts. It reads the compact
% repository snapshots:
%
%   verification/crack_path/reference/accepted_clean_run_P1_P23.csv
%   verification/crack_path/reference/accepted_clean_run_vertices_P0_P23.csv
%
% These CSV files reproduce the accepted clean-run scientific observables
% through P23. They are sufficient for the five scientific figures, but are
% NOT a replacement for the original per-step MAT archive. COD-vs-EDI and
% mesh/solver diagnostics still require final_clean_run.
%
% Usage:
%   P = plot_accepted_clean_run_snapshot;
%
% Options:
%   'SaveFigures' : true (default)
%   'FigureDir'   : default verification/crack_path/reference/plots
%   'Formats'     : {'eps','png'} by default
%   'Visible'     : 'on' or 'off'
%
% Output:
%   Out.stepTable
%   Out.vertices
%   Out.landmarks
%   Out.figures
%   Out.source

    ip=inputParser;
    addParameter(ip,'SaveFigures',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FigureDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Formats',{'eps','png'},@(x)ischar(x)||isstring(x)||iscell(x));
    addParameter(ip,'Visible','on', ...
        @(x)(ischar(x)||isstring(x))&&any(strcmpi(char(x),{'on','off'})));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(mfilename('fullpath'));
    refDir=fullfile(root,'verification','crack_path','reference');

    stepFile=fullfile(refDir,'accepted_clean_run_P1_P23.csv');
    vertexFile=fullfile(refDir,'accepted_clean_run_vertices_P0_P23.csv');

    if exist(stepFile,'file')~=2 || exist(vertexFile,'file')~=2
        error('snapshotplot:MissingReference', ...
            'Accepted clean-run reference CSV files are missing.');
    end

    T=readtable(stepFile);
    V=readtable(vertexFile,'TextType','string');

    required={'segment','tip_x_m','tip_y_m','theta_deg','KI_unit','KII_unit', ...
        'KII_over_KI','delta_theta_next_deg','theta_next_deg'};
    miss=required(~ismember(required,T.Properties.VariableNames));
    if ~isempty(miss)
        error('snapshotplot:StepFields','Reference step table is missing: %s', ...
            strjoin(miss,', '));
    end
    if ~all(ismember({'x_m','y_m'},V.Properties.VariableNames))
        error('snapshotplot:VertexFields','Reference vertex table is incomplete.');
    end

    C=cfg_first_segment_asymmetric();
    geometryOK=abs(C.A-0.30)<=1e-14 && abs(C.B-0.10)<=1e-14 && ...
        strcmpi(C.hole.type,'circle') && ...
        norm(C.hole.center-[0.17,-0.02])<=1e-14 && ...
        abs(C.hole.r-0.030)<=1e-14 && C.hole.npoly==480 && ...
        abs(C.a0-0.004)<=1e-14;
    if ~geometryOK
        error('snapshotplot:FrozenGeometryMismatch', ...
            'Repository configuration no longer matches the accepted geometry fingerprint.');
    end

    vertices=[V.x_m,V.y_m];
    if height(V)~=height(T)+1
        error('snapshotplot:VertexCount','Expected P0 plus one vertex per accepted step.');
    end

    % Independent consistency checks before plotting.
    if any(T.segment(:)~=(1:height(T))')
        error('snapshotplot:SegmentSequence','Accepted segment indices are not consecutive.');
    end
    if max(vecnorm(vertices(2:end,:)-[T.tip_x_m,T.tip_y_m],2,2))>2e-12
        error('snapshotplot:TipMismatch','Step-table tips do not match the vertex snapshot.');
    end
    seg=diff(vertices,1,1);
    if max(abs(vecnorm(seg,2,2)-C.a0))>2e-9
        error('snapshotplot:IncrementMismatch','Reference path is not a uniform 4-mm path.');
    end

    T.crack_length_mm=1e3*T.segment*C.a0;
    L=local_landmarks(T,vertices);

    figDir=char(opt.FigureDir);
    if isempty(strtrim(figDir)),figDir=fullfile(refDir,'plots');end
    if ~local_is_absolute_path(figDir),figDir=fullfile(root,figDir);end
    formats=local_formats(opt.Formats);
    if opt.SaveFigures && exist(figDir,'dir')~=7,mkdir(figDir);end

    F=struct();

    % 1. Trajectory with the physical geometry.
    F.trajectory=local_new_figure(opt.Visible,'Accepted clean-run trajectory');
    ax=axes(F.trajectory); hold(ax,'on'); box(ax,'on');
    local_plot_plate_and_hole(ax,C);
    plot(ax,1e3*vertices(:,1),1e3*vertices(:,2),'-o', ...
        'LineWidth',1.5,'MarkerSize',4,'DisplayName','accepted crack path');
    plot(ax,1e3*vertices(1,1),1e3*vertices(1,2),'s', ...
        'MarkerSize',7,'LineWidth',1.1,'DisplayName','crack mouth');
    local_plot_landmark_points(ax,L);
    xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');
    title(ax,'Accepted crack trajectory through P_{23}');
    axis(ax,'equal'); grid(ax,'on');
    xlim(ax,[130 302]); ylim(ax,[-56 18]);
    legend(ax,'Location','best');
    hold(ax,'off');

    % 2. Absolute local crack angle.
    F.theta=local_new_figure(opt.Visible,'Absolute crack direction');
    ax=axes(F.theta);
    plot(ax,T.crack_length_mm,T.theta_deg,'-o','LineWidth',1.5,'MarkerSize',4);
    hold(ax,'on'); yline(ax,0,'--','HandleVisibility','off');
    local_landmark_lines(ax,L);
    xlabel(ax,'Crack length [mm]'); ylabel(ax,'\theta_k [deg]');
    title(ax,'Absolute local crack direction');
    grid(ax,'on'); box(ax,'on'); hold(ax,'off');

    % 3. Incremental MTS turn.
    F.deltaTheta=local_new_figure(opt.Visible,'Incremental MTS turn');
    ax=axes(F.deltaTheta);
    plot(ax,T.crack_length_mm,T.delta_theta_next_deg,'-o', ...
        'LineWidth',1.5,'MarkerSize',4);
    hold(ax,'on'); yline(ax,0,'--','HandleVisibility','off');
    local_landmark_lines(ax,L);
    xlabel(ax,'Crack length [mm]'); ylabel(ax,'\Delta\theta_{k+1} [deg]');
    title(ax,'Incremental MTS turning angle');
    grid(ax,'on'); box(ax,'on'); hold(ax,'off');

    % 4. Mode-I SIF.
    F.KI=local_new_figure(opt.Visible,'Mode-I SIF');
    ax=axes(F.KI);
    plot(ax,T.crack_length_mm,T.KI_unit,'-o','LineWidth',1.5,'MarkerSize',4);
    hold(ax,'on'); local_landmark_lines(ax,L);
    xlabel(ax,'Crack length [mm]');
    ylabel(ax,'K_I at unit traction [MPa sqrt(m)]');
    title(ax,'Mode-I stress-intensity factor');
    grid(ax,'on'); box(ax,'on'); hold(ax,'off');

    % 5. Mode mixity.
    F.modeMixity=local_new_figure(opt.Visible,'Mode mixity');
    ax=axes(F.modeMixity);
    plot(ax,T.crack_length_mm,T.KII_over_KI,'-o','LineWidth',1.5,'MarkerSize',4);
    hold(ax,'on'); yline(ax,0,'--','HandleVisibility','off');
    local_landmark_lines(ax,L);
    plot(ax,L.aPeak_mm,L.qPeak,'o','MarkerSize',8,'LineWidth',1.5, ...
        'DisplayName',sprintf('P_{%d}: max K_{II}/K_I',L.peakSegment));
    if L.hasZeroCrossing
        plot(ax,L.aLS_mm,0,'d','MarkerSize',8,'LineWidth',1.5, ...
            'DisplayName',sprintf('local symmetry est., %.3f mm',L.aLS_mm));
    end
    xlabel(ax,'Crack length [mm]'); ylabel(ax,'K_{II}/K_I');
    title(ax,'Mode-mixity evolution');
    grid(ax,'on'); box(ax,'on'); legend(ax,'Location','best'); hold(ax,'off');

    if opt.SaveFigures
        local_export(F.trajectory,figDir,'01_trajectory_snapshot',formats);
        local_export(F.theta,figDir,'02_theta_snapshot',formats);
        local_export(F.deltaTheta,figDir,'03_delta_theta_snapshot',formats);
        local_export(F.KI,figDir,'04_KI_snapshot',formats);
        local_export(F.modeMixity,figDir,'05_mode_mixity_snapshot',formats);
    end

    Out=struct();
    Out.stepTable=T;
    Out.vertices=vertices;
    Out.landmarks=L;
    Out.figures=F;
    Out.source=struct( ...
        'mode','repository_scientific_snapshot', ...
        'stepsFile',stepFile, ...
        'verticesFile',vertexFile, ...
        'note',['Scientific P1-P23 snapshot only. Original compact MAT files ', ...
                'remain authoritative for COD, mesh, and solver diagnostics.']);
    Out.figureDir=figDir;

    fprintf('\nACCEPTED CLEAN-RUN SNAPSHOT PLOTTING COMPLETE\n');
    fprintf('  scientific states       : P1 ... P23\n');
    fprintf('  max mode mixity         : P%d at %.3f mm, KII/KI=%+.8g\n', ...
        L.peakSegment,L.aPeak_mm,L.qPeak);
    if L.hasZeroCrossing
        fprintf('  local-symmetry estimate : %.6f mm between P%d and P%d\n', ...
            L.aLS_mm,L.zeroLeftSegment,L.zeroRightSegment);
    end
    fprintf('  last accepted state     : P23 at 92.000 mm\n');
    fprintf('  COD/mesh/solver plots   : unavailable without final_clean_run MAT archive\n');
    if opt.SaveFigures
        fprintf('  figures saved in        : %s\n',figDir);
    end
end


function L=local_landmarks(T,vertices)
    q=T.KII_over_KI;
    [qPeak,iPeak]=max(q);

    L=struct();
    L.peakIndex=iPeak;
    L.peakSegment=T.segment(iPeak);
    L.aPeak_mm=T.crack_length_mm(iPeak);
    L.qPeak=qPeak;
    L.xPeak_mm=1e3*T.tip_x_m(iPeak);
    L.yPeak_mm=1e3*T.tip_y_m(iPeak);

    L.lastSegment=T.segment(end);
    L.aLast_mm=T.crack_length_mm(end);
    L.xLast_mm=1e3*T.tip_x_m(end);
    L.yLast_mm=1e3*T.tip_y_m(end);

    L.hasZeroCrossing=false;
    L.aLS_mm=NaN; L.xLS_mm=NaN; L.yLS_mm=NaN;
    L.zeroLeftSegment=NaN; L.zeroRightSegment=NaN;

    for i=iPeak:height(T)-1
        if q(i)==0
            f=0;
        elseif q(i)*q(i+1)<=0
            f=abs(q(i))/(abs(q(i))+abs(q(i+1)));
        else
            continue
        end
        L.hasZeroCrossing=true;
        L.zeroLeftSegment=T.segment(i);
        L.zeroRightSegment=T.segment(i+1);
        L.aLS_mm=T.crack_length_mm(i)+ ...
            f*(T.crack_length_mm(i+1)-T.crack_length_mm(i));

        pL=vertices(T.segment(i)+1,:);
        pR=vertices(T.segment(i+1)+1,:);
        p=pL+f*(pR-pL);
        L.xLS_mm=1e3*p(1);
        L.yLS_mm=1e3*p(2);
        break
    end
end


function fig=local_new_figure(vis,name)
    fig=figure('Visible',char(vis),'Name',name,'Color','w');
end


function local_plot_plate_and_hole(ax,C)
    A=1e3*C.A; B=1e3*C.B;
    plot(ax,[0 A A 0 0],[-B -B B B -B],'k-','LineWidth',1.0, ...
        'DisplayName','plate boundary');
    ph=linspace(0,2*pi,721);
    x=1e3*(C.hole.center(1)+C.hole.r*cos(ph));
    y=1e3*(C.hole.center(2)+C.hole.r*sin(ph));
    plot(ax,x,y,'k-','LineWidth',1.0,'DisplayName','circular hole');
end


function local_plot_landmark_points(ax,L)
    plot(ax,L.xPeak_mm,L.yPeak_mm,'o','MarkerSize',8,'LineWidth',1.5, ...
        'DisplayName',sprintf('P_{%d}: max mode mixity',L.peakSegment));
    if L.hasZeroCrossing
        plot(ax,L.xLS_mm,L.yLS_mm,'d','MarkerSize',8,'LineWidth',1.5, ...
            'DisplayName',sprintf('local symmetry est. %.2f mm',L.aLS_mm));
    end
    plot(ax,L.xLast_mm,L.yLast_mm,'s','MarkerSize',8,'LineWidth',1.5, ...
        'DisplayName','P_{23}: last accepted');
end


function local_landmark_lines(ax,L)
    xline(ax,L.aPeak_mm,'--',sprintf('P_{%d}',L.peakSegment), ...
        'LabelVerticalAlignment','bottom','HandleVisibility','off');
    if L.hasZeroCrossing
        xline(ax,L.aLS_mm,':','local symmetry', ...
            'LabelVerticalAlignment','middle','HandleVisibility','off');
    end
    xline(ax,L.aLast_mm,'-.','P_{23}', ...
        'LabelVerticalAlignment','top','HandleVisibility','off');
end


function formats=local_formats(x)
    if ischar(x)
        formats={lower(strtrim(x))};
    elseif isstring(x)
        formats=cellstr(lower(strtrim(x(:))));
    else
        formats=cellfun(@(s)lower(strtrim(char(s))),x,'UniformOutput',false);
    end
    formats=formats(~cellfun(@isempty,formats));
    allowed={'eps','png','pdf','fig'};
    bad=formats(~ismember(formats,allowed));
    if ~isempty(bad)
        error('snapshotplot:Format','Unsupported figure format(s): %s', ...
            strjoin(bad,', '));
    end
end


function local_export(fig,folder,stem,formats)
    for j=1:numel(formats)
        fmt=formats{j};
        file=fullfile(folder,[stem '.' fmt]);
        switch fmt
            case 'eps'
                print(fig,file,'-depsc','-painters');
            case 'png'
                exportgraphics(fig,file,'Resolution',300);
            case 'pdf'
                exportgraphics(fig,file,'ContentType','vector');
            case 'fig'
                savefig(fig,file);
        end
    end
end


function tf=local_is_absolute_path(p)
    p=char(p);
    if isempty(p),tf=false;return,end
    tf=startsWith(p,filesep) || startsWith(p,'\\') || ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'));
end
