function Out = plot_final_clean_run_additional_results(varargin)
%PLOT_FINAL_CLEAN_RUN_ADDITIONAL_RESULTS Additional final clean-run plots.
%
%   Out = plot_final_clean_run_additional_results()
%
% This routine is intentionally postprocessing-only. It reads the audited
% atomic path state plus compact per-step physical/qualification files from
% verification/crack_path/final_clean_run. It performs no FE solve, no mesh
% generation, and no crack-path continuation.
%
% Additional scientific figures
%   8.  KII versus crack length.
%   9.  KI and KII in aligned panels.
%   10. Late-path mode mixity and MTS turn in aligned panels.
%   11. COD sensitivity of the MTS turning prediction for all stored fits.
%   12. COD sensitivity of KII/KI for all stored fits.
%   13. Publication-style trajectory overview plus late-path detail.
%
% Name-value options
%   'RunDir'           : clean-run directory.
%   'SaveFigures'      : default true.
%   'FigureDir'        : default <RunDir>/plots.
%   'Formats'          : default {'eps','png'}.
%   'Visible'          : 'on' or 'off', default 'on'.
%   'CloseExisting'    : close figures created here first, default false.
%   'UseLatex'          : use LaTeX interpreters for all figure text, default true.
%   'LateStartSegment' : first segment in late-path zoom, default 15.
%   'PlateA'           : plate width [m], default 0.300.
%   'PlateB'           : plate half-height [m], default 0.100.
%   'HoleRadius'       : hole radius [m], default 0.030.
%
% Output
%   Out.stepTable
%   Out.codSensitivityTable
%   Out.landmarks
%   Out.figures
%   Out.validation
%
% The hole center for the trajectory figure is reconstructed from the
% audited crack mouth, the prescribed first-segment direction, and the
% supplied HoleRadius. For the accepted asymmetric clean run this gives
% [0.170,-0.020] m.

    ip=inputParser;
    addParameter(ip,'RunDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'SaveFigures',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FigureDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Formats',{'eps','png'},@(x)ischar(x)||isstring(x)||iscell(x));
    addParameter(ip,'Visible','on', ...
        @(x)(ischar(x)||isstring(x))&&any(strcmpi(char(x),{'on','off'})));
    addParameter(ip,'CloseExisting',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'UseLatex',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'LateStartSegment',15, ...
        @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x==round(x)&&x>=1);
    addParameter(ip,'PlateA',0.300,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'PlateB',0.100,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(ip,'HoleRadius',0.030,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(mfilename('fullpath'));

    runDir=char(opt.RunDir);
    if isempty(strtrim(runDir))
        runDir=fullfile(root,'verification','crack_path','final_clean_run');
    elseif ~local_is_absolute_path(runDir)
        runDir=fullfile(root,runDir);
    end
    if exist(runDir,'dir')~=7
        error('cleanplot2:MissingRunDir','Run directory not found: %s',runDir);
    end

    stateFile=fullfile(runDir,'path_run_state.mat');
    if exist(stateFile,'file')~=2
        error('cleanplot2:MissingState','Atomic path state not found: %s',stateFile);
    end
    d=load(stateFile,'State');
    if ~isfield(d,'State')||~isstruct(d.State)||~isscalar(d.State)
        error('cleanplot2:BadState','path_run_state.mat must contain scalar struct State.');
    end
    State=d.State;

    req={'vertices','thetaDeg','completedPhysicalSegments', ...
        'rowsThroughCompleted','rowVariableNames','rowHistoryComplete','fastEDI'};
    for j=1:numel(req)
        if ~isfield(State,req{j})
            error('cleanplot2:StateField','State is missing field %s.',req{j});
        end
    end
    if ~logical(State.rowHistoryComplete)
        error('cleanplot2:IncompleteHistory', ...
            'Additional plots require the complete accepted row history.');
    end

    kDone=double(State.completedPhysicalSegments);
    nGeom=size(State.vertices,1)-1;
    if kDone<2 || kDone~=round(kDone) || nGeom<kDone
        error('cleanplot2:StateShape','Invalid accepted/geometric segment counts.');
    end

    names=cellstr(string(State.rowVariableNames));
    rows=State.rowsThroughCompleted(1:kDone,:);
    if size(rows,2)~=numel(names)
        error('cleanplot2:RowShape','rowsThroughCompleted width does not match rowVariableNames.');
    end
    T=array2table(rows,'VariableNames',names);

    need={'segment','tip_x_m','tip_y_m','theta_deg','KI_unit','KII_unit', ...
        'KII_over_KI','delta_theta_next_deg','theta_next_deg', ...
        'PCG_iterations','PCG_relres','true_rel_residual','EDI_elements'};
    miss=need(~ismember(need,T.Properties.VariableNames));
    if ~isempty(miss)
        error('cleanplot2:StepFields','Step table is missing: %s',strjoin(miss,', '));
    end

    % The production increment is recovered from the accepted geometry, not
    % from a missing historical Stage-I MAT file.
    verticesAccepted=State.vertices(1:kDone+1,:);
    seg=diff(verticesAccepted,1,1);
    segLength=vecnorm(seg,2,2);
    da=median(segLength);
    if max(abs(segLength-da))>2e-12
        error('cleanplot2:Increment','Accepted path does not retain a uniform increment.');
    end
    T.crack_length_mm=1e3*T.segment*da;

    % Verify exact compact physical history and collect all stored COD fits.
    [COD,physicalCount]=local_collect_all_cod(runDir,T);

    % Characteristic states from authoritative accepted data.
    L=local_landmarks(T,State.vertices);

    % Qualified-unsolved next geometry is diagnostic only.
    [hasQualifiedUnsolved,QnextPass]=local_unsolved_qualification(runDir,nGeom,kDone);

    if opt.CloseExisting
        close(findall(groot,'Type','figure','Tag','CrackPathAdditionalPlot'));
    end

    figDir=char(opt.FigureDir);
    if isempty(strtrim(figDir)),figDir=fullfile(runDir,'plots');end
    if ~local_is_absolute_path(figDir),figDir=fullfile(root,figDir);end
    formats=local_formats(opt.Formats);
    if opt.SaveFigures && exist(figDir,'dir')~=7,mkdir(figDir);end

    lateStart=max(1,min(kDone,double(opt.LateStartSegment)));
    lateMask=T.segment>=lateStart;

    F=struct();

    % ==================================================================
    % 8. KII versus crack length
    % ==================================================================
    F.KII=local_new_figure(opt.Visible,'Mode-II SIF');
    ax=axes(F.KII);
    plot(ax,T.crack_length_mm,T.KII_unit,'-o','LineWidth',1.5,'MarkerSize',4, ...
        'DisplayName','K_{II}');
    hold(ax,'on');
    yline(ax,0,'--','HandleVisibility','off');
    local_landmark_lines(ax,L);
    plot(ax,L.aKIImax_mm,L.KIImax,'o','MarkerSize',8,'LineWidth',1.5, ...
        'DisplayName',sprintf('P_{%d}: max K_{II}',L.KIImaxSegment));
    if L.hasZeroCrossing
        plot(ax,L.aLS_mm,0,'d','MarkerSize',8,'LineWidth',1.5, ...
            'DisplayName',sprintf('local symmetry, %.3f mm',L.aLS_mm));
    end
    plot(ax,L.aLast_mm,T.KII_unit(end),'s','MarkerSize',8,'LineWidth',1.5, ...
        'DisplayName',sprintf('P_{%d}: last accepted',L.lastSegment));
    xlabel(ax,'Crack length [mm]');
    ylabel(ax,'K_{II} at unit traction [MPa sqrt(m)]');
    title(ax,'Mode-II stress-intensity factor');
    grid(ax,'on');box(ax,'on');legend(ax,'Location','best');hold(ax,'off');

    % ==================================================================
    % 9. KI and KII in aligned panels
    % ==================================================================
    F.KIKII=local_new_figure(opt.Visible,'KI and KII evolution');
    tl=tiledlayout(F.KIKII,2,1,'TileSpacing','compact','Padding','compact');

    ax1=nexttile(tl);
    plot(ax1,T.crack_length_mm,T.KI_unit,'-o','LineWidth',1.4,'MarkerSize',4);
    hold(ax1,'on');local_landmark_lines(ax1,L);
    ylabel(ax1,'K_I [MPa sqrt(m)]');
    title(ax1,'Mode-I intensity: monotone loading response');
    grid(ax1,'on');box(ax1,'on');hold(ax1,'off');

    ax2=nexttile(tl);
    plot(ax2,T.crack_length_mm,T.KII_unit,'-o','LineWidth',1.4,'MarkerSize',4);
    hold(ax2,'on');yline(ax2,0,'--','HandleVisibility','off');
    local_landmark_lines(ax2,L);
    xlabel(ax2,'Crack length [mm]');
    ylabel(ax2,'K_{II} [MPa sqrt(m)]');
    title(ax2,'Mode-II intensity: sign reversal controls redirection');
    grid(ax2,'on');box(ax2,'on');hold(ax2,'off');

    % ==================================================================
    % 10. Late-path mode mixity and MTS turn
    % ==================================================================
    F.lateMixityTurn=local_new_figure(opt.Visible,'Late-path mode mixity and MTS turn');
    tl=tiledlayout(F.lateMixityTurn,2,1,'TileSpacing','compact','Padding','compact');

    ax1=nexttile(tl);
    plot(ax1,T.crack_length_mm(lateMask),T.KII_over_KI(lateMask), ...
        '-o','LineWidth',1.5,'MarkerSize',5);
    hold(ax1,'on');yline(ax1,0,'--','HandleVisibility','off');
    local_landmark_lines(ax1,L);
    ylabel(ax1,'K_{II}/K_I');
    title(ax1,sprintf('Late-path mode mixity, P_{%d}--P_{%d}',lateStart,kDone));
    grid(ax1,'on');box(ax1,'on');hold(ax1,'off');

    ax2=nexttile(tl);
    plot(ax2,T.crack_length_mm(lateMask),T.delta_theta_next_deg(lateMask), ...
        '-o','LineWidth',1.5,'MarkerSize',5);
    hold(ax2,'on');yline(ax2,0,'--','HandleVisibility','off');
    local_landmark_lines(ax2,L);
    xlabel(ax2,'Crack length [mm]');
    ylabel(ax2,'\Delta\theta_{k+1} [deg]');
    title(ax2,'MTS turn: sign reversal follows K_{II}');
    grid(ax2,'on');box(ax2,'on');hold(ax2,'off');

    % ==================================================================
    % 11. COD sensitivity of MTS turning angle
    % ==================================================================
    degrees=unique(COD.degree,'sorted');
    windows=unique(COD{:,{'lower_r_over_DeltaA','upper_r_over_DeltaA'}}, ...
        'rows','stable');

    F.CODTurnSensitivity=local_new_figure(opt.Visible,'COD MTS sensitivity');
    tl=tiledlayout(F.CODTurnSensitivity,numel(degrees),1, ...
        'TileSpacing','compact','Padding','compact');

    for id=1:numel(degrees)
        deg=degrees(id);
        ax=nexttile(tl);
        plot(ax,T.crack_length_mm(2:end),T.delta_theta_next_deg(2:end), ...
            '-o','LineWidth',1.6,'MarkerSize',4,'DisplayName','EDI');
        hold(ax,'on');
        local_plot_cod_family(ax,COD,deg,windows,'delta_theta_COD_deg');
        yline(ax,0,'--','HandleVisibility','off');
        local_landmark_lines(ax,L);
        ylabel(ax,'\Delta\theta_{k+1} [deg]');
        title(ax,sprintf('COD MTS sensitivity, polynomial degree %d',deg));
        grid(ax,'on');box(ax,'on');legend(ax,'Location','best');
        if id==numel(degrees),xlabel(ax,'Crack length [mm]');end
        hold(ax,'off');
    end

    % ==================================================================
    % 12. COD sensitivity of mode mixity
    % ==================================================================
    F.CODMixitySensitivity=local_new_figure(opt.Visible,'COD mode-mixity sensitivity');
    tl=tiledlayout(F.CODMixitySensitivity,numel(degrees),1, ...
        'TileSpacing','compact','Padding','compact');

    for id=1:numel(degrees)
        deg=degrees(id);
        ax=nexttile(tl);
        plot(ax,T.crack_length_mm(2:end),T.KII_over_KI(2:end), ...
            '-o','LineWidth',1.6,'MarkerSize',4,'DisplayName','EDI');
        hold(ax,'on');
        local_plot_cod_family(ax,COD,deg,windows,'ratio_COD');
        yline(ax,0,'--','HandleVisibility','off');
        local_landmark_lines(ax,L);
        ylabel(ax,'K_{II}/K_I');
        title(ax,sprintf('COD mode-mixity sensitivity, polynomial degree %d',deg));
        grid(ax,'on');box(ax,'on');legend(ax,'Location','best');
        if id==numel(degrees),xlabel(ax,'Crack length [mm]');end
        hold(ax,'off');
    end

    % ==================================================================
    % 13. Publication-style trajectory overview and true-scale late detail
    %
    % The panels are stacked deliberately. A side-by-side late-path panel
    % is too narrow for an equal-scale x-y plot: the axes then expand and
    % suppress the small but meaningful curvature. Both panels below keep
    % 1 mm = 1 mm. The accepted path ends at P23; the qualified-but-unsolved
    % P23->P24 segment is shown separately as dashed diagnostic geometry.
    % ==================================================================
    F.trajectoryDetail=local_new_figure(opt.Visible, ...
        'Trajectory overview and true-scale late detail');
    set(F.trajectoryDetail,'Position',[100 100 1000 560], ...
        'PaperPositionMode','auto');
    tl=tiledlayout(F.trajectoryDetail,2,1, ...
        'TileSpacing','compact','Padding','compact');

    allVertices=State.vertices;
    p0=allVertices(1,:);
    e1=allVertices(2,:)-allVertices(1,:);
    e1=e1/norm(e1);
    holeCenter=p0-opt.HoleRadius*e1;

    % ----- overview: actual geometry at true scale
    ax1=nexttile(tl);
    local_plot_geometry(ax1,opt.PlateA,opt.PlateB,holeCenter,opt.HoleRadius);
    hold(ax1,'on');
    plot(ax1,1e3*verticesAccepted(:,1),1e3*verticesAccepted(:,2), ...
        '-o','LineWidth',1.5,'MarkerSize',4, ...
        'DisplayName',sprintf('accepted path, P_0--P_{%d}',kDone));
    local_plot_path_landmarks(ax1,L);

    if hasQualifiedUnsolved && QnextPass
        plot(ax1,1e3*allVertices(end-1:end,1),1e3*allVertices(end-1:end,2), ...
            '--','LineWidth',1.4, ...
            'DisplayName',sprintf('P_{%d}--P_{%d}: qualified, unsolved', ...
            kDone,nGeom));
        plot(ax1,1e3*allVertices(end,1),1e3*allVertices(end,2),'x', ...
            'MarkerSize',9,'LineWidth',1.5,'HandleVisibility','off');
    end

    xlabel(ax1,'x [mm]');ylabel(ax1,'y [mm]');
    title(ax1,'Crack trajectory: accepted path and next qualified geometry');
    grid(ax1,'on');box(ax1,'on');
    daspect(ax1,[1 1 1]);
    [xmin,xmax,ymin,ymax]=local_overview_limits(opt.PlateA,opt.PlateB, ...
        holeCenter,opt.HoleRadius,allVertices);
    xlim(ax1,[xmin xmax]);ylim(ax1,[ymin ymax]);
    legend(ax1,'Location','eastoutside');
    hold(ax1,'off');

    % ----- late detail: same physical x-y scale, but in a wide strip
    ax2=nexttile(tl);
    hold(ax2,'on');box(ax2,'on');

    lateVertexStart=max(1,lateStart+1); % row corresponding to P_lateStart
    lateAccepted=verticesAccepted(lateVertexStart:end,:);
    plot(ax2,1e3*lateAccepted(:,1),1e3*lateAccepted(:,2), ...
        '-o','LineWidth',1.7,'MarkerSize',5,'HandleVisibility','off');

    xline(ax2,1e3*opt.PlateA,'--','right boundary', ...
        'LabelVerticalAlignment','middle','HandleVisibility','off');
    local_plot_path_landmarks(ax2,L);

    if hasQualifiedUnsolved && QnextPass
        plot(ax2,1e3*allVertices(end-1:end,1),1e3*allVertices(end-1:end,2), ...
            '--','LineWidth',1.5,'HandleVisibility','off');
        plot(ax2,1e3*allVertices(end,1),1e3*allVertices(end,2),'x', ...
            'MarkerSize',10,'LineWidth',1.7,'HandleVisibility','off');
    end

    lateVertices=allVertices(lateVertexStart:end,:);
    xlo=min(1e3*lateVertices(:,1))-1.5;
    xhi=1e3*opt.PlateA+1.0;
    ylo=min(1e3*lateVertices(:,2))-1.0;
    yhi=max(1e3*lateVertices(:,2))+1.0;

    xlabel(ax2,'x [mm]');ylabel(ax2,'y [mm]');
    title(ax2,sprintf( ...
        'True-scale late-path detail, P_{%d}--P_{%d}; P_{%d} is unsolved', ...
        lateStart,nGeom,nGeom));
    grid(ax2,'on');
    daspect(ax2,[1 1 1]);
    xlim(ax2,[xlo xhi]);ylim(ax2,[ylo yhi]);
    hold(ax2,'off');

    if opt.UseLatex
        figNames=fieldnames(F);
        for j=1:numel(figNames)
            local_apply_latex(F.(figNames{j}));
        end
    end

    if opt.SaveFigures
        local_export(F.KII,figDir,'08_KII',formats);
        local_export(F.KIKII,figDir,'09_KI_KII_panels',formats);
        local_export(F.lateMixityTurn,figDir,'10_late_mode_mixity_MTS',formats);
        local_export(F.CODTurnSensitivity,figDir,'11_COD_turn_sensitivity',formats);
        local_export(F.CODMixitySensitivity,figDir,'12_COD_mode_mixity_sensitivity',formats);
        local_export(F.trajectoryDetail,figDir,'13_trajectory_late_detail',formats);
    end

    V=struct();
    V.physicalCompactFiles=physicalCount;
    V.expectedPhysicalCompactFiles=kDone-1;
    V.codFitRows=height(COD);
    V.expectedCODFitRows=8*(kDone-1);
    V.allPhysicalFilesPresent=physicalCount==(kDone-1);
    V.allStoredCODFitsPresent=height(COD)==8*(kDone-1);
    V.hasQualifiedUnsolvedNext=hasQualifiedUnsolved&&QnextPass;
    V.increment_m=da;
    V.holeCenter_m=holeCenter;
    V.holeRadius_m=opt.HoleRadius;
    V.UseLatex=opt.UseLatex;

    Out=struct();
    Out.stepTable=T;
    Out.codSensitivityTable=COD;
    Out.landmarks=L;
    Out.figures=F;
    Out.validation=V;
    Out.runDir=runDir;
    Out.figureDir=figDir;

    fprintf('\nADDITIONAL FINAL CLEAN-RUN PLOTTING COMPLETE\n');
    fprintf('  accepted physical states : P1 ... P%d\n',kDone);
    fprintf('  physical compact files   : %d / %d\n', ...
        physicalCount,kDone-1);
    fprintf('  stored COD fit rows      : %d / %d\n', ...
        height(COD),8*(kDone-1));
    fprintf('  max positive KII         : P%d at %.3f mm, KII=%+.12g\n', ...
        L.KIImaxSegment,L.aKIImax_mm,L.KIImax);
    fprintf('  max mode mixity          : P%d at %.3f mm, KII/KI=%+.12g\n', ...
        L.peakSegment,L.aPeak_mm,L.qPeak);
    if L.hasZeroCrossing
        fprintf('  local-symmetry estimate  : %.12f mm between P%d and P%d\n', ...
            L.aLS_mm,L.zeroLeftSegment,L.zeroRightSegment);
    end
    if hasQualifiedUnsolved&&QnextPass
        fprintf('  qualified-unsolved state : P%d shown only as diagnostic geometry\n',nGeom);
    end
    fprintf('  reconstructed hole ctr   : [%.12g, %.12g] m\n', ...
        holeCenter(1),holeCenter(2));
    if opt.SaveFigures
        fprintf('  figures saved in         : %s\n',figDir);
    end
end


% =========================================================================
function [COD,nPhysical]=local_collect_all_cod(runDir,T)
    nExpected=8*max(0,height(T)-1);
    data=nan(nExpected,13);
    nPhysical=0;
    nRow=0;

    for i=2:height(T)
        k=T.segment(i);
        f=fullfile(runDir,sprintf('step_%03d_physical_small.mat',k));
        if exist(f,'file')~=2
            error('cleanplot2:MissingPhysical','Missing compact physical file for P%d.',k);
        end
        d=load(f,'R');
        if ~isfield(d,'R')||~isstruct(d.R)||~isscalar(d.R)
            error('cleanplot2:BadPhysical','Physical compact file P%d lacks scalar R.',k);
        end
        R=d.R;
        if ~isfield(R,'pass')||~logical(R.pass)|| ...
                ~isfield(R,'fitTable')||~istable(R.fitTable)
            error('cleanplot2:BadPhysical','Physical compact file P%d is not accepted.',k);
        end
        nPhysical=nPhysical+1;

        F=R.fitTable;
        need={'lower_r_over_DeltaA','upper_r_over_DeltaA','degree','KI_COD', ...
            'KII_COD','ratio_COD','RMSE_KII','median_pointwise_ratio', ...
            'delta_theta_next_MTS_deg'};
        miss=need(~ismember(need,F.Properties.VariableNames));
        if ~isempty(miss)
            error('cleanplot2:CODFields','P%d fit table is missing: %s', ...
                k,strjoin(miss,', '));
        end
        if height(F)~=8
            error('cleanplot2:CODCount','Expected 8 stored COD fits at P%d, found %d.', ...
                k,height(F));
        end

        for j=1:height(F)
            nRow=nRow+1;
            data(nRow,:)=[ ...
                k,T.crack_length_mm(i),F.lower_r_over_DeltaA(j), ...
                F.upper_r_over_DeltaA(j),F.degree(j), ...
                F.KI_COD(j),F.KII_COD(j),F.ratio_COD(j),F.RMSE_KII(j), ...
                F.median_pointwise_ratio(j),F.delta_theta_next_MTS_deg(j), ...
                T.delta_theta_next_deg(i),T.KII_over_KI(i)];
        end
    end

    data=data(1:nRow,:);
    COD=array2table(data,'VariableNames',{ ...
        'segment','crack_length_mm','lower_r_over_DeltaA','upper_r_over_DeltaA', ...
        'degree','KI_COD','KII_COD','ratio_COD','RMSE_KII', ...
        'median_pointwise_ratio','delta_theta_COD_deg', ...
        'delta_theta_EDI_deg','ratio_EDI'});
end


function L=local_landmarks(T,vertices)
    q=T.KII_over_KI;
    a=T.crack_length_mm;

    finiteQ=isfinite(q);
    if ~any(finiteQ)
        error('cleanplot2:NoModeMixity','No finite KII/KI values are available.');
    end

    idx=find(finiteQ);
    [qPeak,ii]=max(q(idx));
    iPeak=idx(ii);

    [KIImax,iKII]=max(T.KII_unit);

    L=struct();
    L.peakIndex=iPeak;
    L.peakSegment=T.segment(iPeak);
    L.aPeak_mm=a(iPeak);
    L.qPeak=qPeak;
    L.xPeak_mm=1e3*T.tip_x_m(iPeak);
    L.yPeak_mm=1e3*T.tip_y_m(iPeak);

    L.KIImaxSegment=T.segment(iKII);
    L.aKIImax_mm=a(iKII);
    L.KIImax=KIImax;

    L.lastIndex=height(T);
    L.lastSegment=T.segment(end);
    L.aLast_mm=a(end);
    L.xLast_mm=1e3*T.tip_x_m(end);
    L.yLast_mm=1e3*T.tip_y_m(end);

    L.hasZeroCrossing=false;
    L.aLS_mm=NaN;
    L.xLS_mm=NaN;
    L.yLS_mm=NaN;
    L.zeroLeftSegment=NaN;
    L.zeroRightSegment=NaN;

    for i=max(1,iPeak):height(T)-1
        if ~isfinite(q(i))||~isfinite(q(i+1)),continue,end
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
        L.aLS_mm=a(i)+f*(a(i+1)-a(i));

        pL=vertices(T.segment(i)+1,:);
        pR=vertices(T.segment(i+1)+1,:);
        p=pL+f*(pR-pL);
        L.xLS_mm=1e3*p(1);
        L.yLS_mm=1e3*p(2);
        break
    end
end


function [hasNext,passNext]=local_unsolved_qualification(runDir,nGeom,kDone)
    hasNext=nGeom>kDone;
    passNext=false;
    if ~hasNext,return,end

    f=fullfile(runDir,sprintf('step_%03d_qualification_small.mat',nGeom));
    if exist(f,'file')~=2,return,end
    d=load(f,'Small');
    if isfield(d,'Small')&&isstruct(d.Small)&&isfield(d.Small,'pass')
        passNext=logical(d.Small.pass);
    end

    if exist(fullfile(runDir,sprintf('step_%03d_physical_small.mat',nGeom)),'file')==2 || ...
       exist(fullfile(runDir,sprintf('step_%03d_physical_solved.mat',nGeom)),'file')==2
        error('cleanplot2:UnexpectedNextPhysical', ...
            'P%d is not qualified-unsolved: a physical result/checkpoint exists.',nGeom);
    end
end


function local_plot_cod_family(ax,COD,degree,windows,varName)
    lineStyles={'-s','-^','-v','-d','-x','-+'};
    for iw=1:size(windows,1)
        mask=COD.degree==degree & ...
            abs(COD.lower_r_over_DeltaA-windows(iw,1))<=1e-12 & ...
            abs(COD.upper_r_over_DeltaA-windows(iw,2))<=1e-12;
        C=COD(mask,:);
        [~,ord]=sort(C.segment);
        C=C(ord,:);
        style=lineStyles{1+mod(iw-1,numel(lineStyles))};
        plot(ax,C.crack_length_mm,C.(varName),style, ...
            'LineWidth',1.1,'MarkerSize',4, ...
            'DisplayName',sprintf('COD [%.2f, %.2f] Delta a', ...
                windows(iw,1),windows(iw,2)));
    end
end


function local_landmark_lines(ax,L)
    xline(ax,L.aPeak_mm,'--',sprintf('P_{%d}',L.peakSegment), ...
        'LabelVerticalAlignment','bottom','HandleVisibility','off');
    if L.hasZeroCrossing
        xline(ax,L.aLS_mm,':','local symmetry', ...
            'LabelVerticalAlignment','middle','HandleVisibility','off');
    end
    xline(ax,L.aLast_mm,'-.',sprintf('P_{%d}',L.lastSegment), ...
        'LabelVerticalAlignment','top','HandleVisibility','off');
end


function local_plot_path_landmarks(ax,L)
    plot(ax,L.xPeak_mm,L.yPeak_mm,'o','MarkerSize',8,'LineWidth',1.5, ...
        'DisplayName',sprintf('P_{%d}: max mode mixity',L.peakSegment));
    if L.hasZeroCrossing
        plot(ax,L.xLS_mm,L.yLS_mm,'d','MarkerSize',8,'LineWidth',1.5, ...
            'DisplayName',sprintf('linear K_{II}/K_I=0 estimate, %.2f mm',L.aLS_mm));
    end
    plot(ax,L.xLast_mm,L.yLast_mm,'s','MarkerSize',8,'LineWidth',1.5, ...
        'DisplayName',sprintf('P_{%d}: last accepted',L.lastSegment));
end


function local_plot_geometry(ax,A,B,holeCenter,R)
    hold(ax,'on');
    A=1e3*A;B=1e3*B;
    plot(ax,[0 A A 0 0],[-B -B B B -B],'k-','LineWidth',1.0, ...
        'DisplayName','plate boundary');

    ph=linspace(0,2*pi,721);
    x=1e3*(holeCenter(1)+R*cos(ph));
    y=1e3*(holeCenter(2)+R*sin(ph));
    plot(ax,x,y,'k-','LineWidth',1.0,'DisplayName','circular hole');
end


function [xmin,xmax,ymin,ymax]=local_overview_limits(A,B,holeCenter,R,vertices)
    xx=1e3*vertices(:,1);
    yy=1e3*vertices(:,2);
    hx=1e3*[holeCenter(1)-R,holeCenter(1)+R];
    hy=1e3*[holeCenter(2)-R,holeCenter(2)+R];
    margin=5;

    xmin=max(0,min([xx;hx(:)])-margin);
    xmax=min(1e3*A+2,max([xx;1e3*A])+2);
    ymin=max(-1e3*B,min([yy;hy(:)])-margin);
    ymax=min( 1e3*B,max([yy;hy(:)])+margin);
end


function fig=local_new_figure(vis,name)
    fig=figure('Visible',char(vis),'Name',name,'Color','w', ...
        'Tag','CrackPathAdditionalPlot');
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
        error('cleanplot2:Format','Unsupported figure format(s): %s',strjoin(bad,', '));
    end
end


function local_apply_latex(fig)
    if isempty(fig) || ~isgraphics(fig),return,end

    tickObjs=findall(fig,'-property','TickLabelInterpreter');
    for k=1:numel(tickObjs)
        try
            set(tickObjs(k),'TickLabelInterpreter','latex');
        catch
        end
    end

    interpObjs=findall(fig,'-property','Interpreter');
    for k=1:numel(interpObjs)
        try
            set(interpObjs(k),'Interpreter','latex');
        catch
        end
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
