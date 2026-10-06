function Out = plot_final_clean_run_results(varargin)
%PLOT_FINAL_CLEAN_RUN_RESULTS Plot the accepted clean incremental path.
%
%   Out = plot_final_clean_run_results('FrozenState',R0)
%
% Reads the atomic path state and the compact per-step qualification/physical
% results from verification/crack_path/final_clean_run. Only accepted
% physical states are used in the scientific plots. An already-qualified but
% unsolved next candidate (for example P24 after a PCG stop) is retained only
% in the numerical-quality diagnostics.
%
% Scientific figures
%   1. Crack trajectory with plate/hole geometry.
%   2. Absolute local crack angle theta_k.
%   3. Incremental MTS turn Delta theta_{k+1}.
%   4. Mode-I stress-intensity factor K_I.
%   5. Mode mixity K_II/K_I.
%
% Verification figures
%   6. EDI-vs-COD MTS prediction using one fixed COD fit definition.
%   7. Numerical-quality diagnostics: physical clearance, minimum T3 angle,
%      maximum adjacent-size ratio, and PCG iterations.
%
% Characteristic states are identified from the data, not hard coded:
%   - positive mode-mixity maximum;
%   - first subsequent K_II/K_I zero crossing (linear interpolation);
%   - last accepted physical state.
%
% Name-value options
%   'FrozenState'     : accepted Stage-I R0 struct (preferred).
%   'StateFile'       : Stage-I MAT containing R0, used if FrozenState empty.
%   'RunDir'          : clean-run directory.
%   'CODWindow'       : [rmin rmax]/Delta a, default [0.08 0.30].
%   'CODDegree'       : polynomial degree, default 2.
%   'SaveFigures'     : default true.
%   'FigureDir'       : output directory, default <RunDir>/plots.
%   'Formats'         : e.g. {'eps','png','pdf'}, default {'eps','png'}.
%   'Visible'         : 'on' or 'off', default 'on'.
%   'CloseExisting'   : close figures created by this function first, default false.
%
% Output
%   Out.stepTable
%   Out.codTable
%   Out.meshTable
%   Out.landmarks
%   Out.figures
%   Out.acceptedVertices
%   Out.acceptedPhysicalSegments
%   Out.geometricSegments
%   Out.validation
%
% The P1 row is the frozen seed and has no physical PCG/COD result. COD and
% solver diagnostics therefore begin at P2.

    ip=inputParser;
    addParameter(ip,'FrozenState',[],@(x)isempty(x)||(isstruct(x)&&isscalar(x)));
    addParameter(ip,'StateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'RunDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'CODWindow',[0.08 0.30], ...
        @(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)>=0&&x(2)>x(1));
    addParameter(ip,'CODDegree',2, ...
        @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x==round(x)&&x>=1);
    addParameter(ip,'SaveFigures',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'FigureDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Formats',{'eps','png'}, ...
        @(x)ischar(x)||isstring(x)||iscell(x));
    addParameter(ip,'Visible','on', ...
        @(x)(ischar(x)||isstring(x))&&any(strcmpi(char(x),{'on','off'})));
    addParameter(ip,'CloseExisting',false,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    root=fileparts(mfilename('fullpath'));
    R0=local_load_frozen_state(root,opt.FrozenState,char(opt.StateFile));

    runDir=char(opt.RunDir);
    if isempty(strtrim(runDir))
        runDir=fullfile(root,'verification','crack_path','final_clean_run');
    elseif ~local_is_absolute_path(runDir)
        runDir=fullfile(root,runDir);
    end
    if exist(runDir,'dir')~=7
        error('cleanplot:MissingRunDir','Run directory not found: %s',runDir);
    end

    stateFile=fullfile(runDir,'path_run_state.mat');
    if exist(stateFile,'file')~=2
        error('cleanplot:MissingState','Atomic path state not found: %s',stateFile);
    end
    d=load(stateFile,'State');
    if ~isfield(d,'State')||~isstruct(d.State)||~isscalar(d.State)
        error('cleanplot:BadState','path_run_state.mat must contain scalar struct State.');
    end
    State=d.State;

    req={'vertices','completedPhysicalSegments','rowsThroughCompleted','rowVariableNames'};
    for j=1:numel(req)
        if ~isfield(State,req{j})
            error('cleanplot:StateField','State is missing field %s.',req{j});
        end
    end

    kDone=double(State.completedPhysicalSegments);
    nGeom=size(State.vertices,1)-1;
    if kDone<1 || kDone~=round(kDone) || nGeom<kDone
        error('cleanplot:StateShape','Invalid accepted/geometric segment counts.');
    end

    vertices=State.vertices(1:kDone+1,:);
    names=cellstr(string(State.rowVariableNames));
    rows=State.rowsThroughCompleted(1:kDone,:);
    if size(rows,2)~=numel(names)
        error('cleanplot:RowShape','rowsThroughCompleted width does not match rowVariableNames.');
    end
    T=array2table(rows,'VariableNames',names);

    required={'segment','tip_x_m','tip_y_m','theta_deg','KI_unit','KII_unit', ...
        'KII_over_KI','delta_theta_next_deg','theta_next_deg', ...
        'PCG_iterations','PCG_relres','true_rel_residual','EDI_elements'};
    miss=required(~ismember(required,T.Properties.VariableNames));
    if ~isempty(miss)
        error('cleanplot:StepFields','Step table is missing: %s',strjoin(miss,', '));
    end

    S0=R0.summary;
    if ~istable(S0)||height(S0)~=1|| ...
            ~ismember('a0_reserved_m',S0.Properties.VariableNames)
        error('cleanplot:FrozenSummary','Frozen R0.summary with a0_reserved_m is required.');
    end
    da=S0.a0_reserved_m;
    T.crack_length_mm=1e3*T.segment*da;

    % Per-step compact physical verification data.
    COD=local_collect_cod(runDir,T,opt.CODWindow,opt.CODDegree);

    % Qualification diagnostics, including a possible next unsolved candidate.
    Mesh=local_collect_mesh(runDir,nGeom,kDone,da);

    % Characteristic states from accepted data.
    L=local_landmarks(T,vertices);

    % Validation summary for the selected COD definition.
    goodCOD=isfinite(COD.delta_theta_COD_deg)&isfinite(COD.delta_theta_EDI_deg);
    if any(goodCOD)
        diffCOD=COD.delta_theta_COD_deg(goodCOD)-COD.delta_theta_EDI_deg(goodCOD);
        codRMS=sqrt(mean(diffCOD.^2));
        codMax=max(abs(diffCOD));
    else
        codRMS=NaN;
        codMax=NaN;
    end
    Validation=struct();
    Validation.CODWindow=opt.CODWindow(:).';
    Validation.CODDegree=opt.CODDegree;
    Validation.CODminusEDI_RMS_deg=codRMS;
    Validation.CODminusEDI_maxAbs_deg=codMax;
    Validation.nCODcomparisons=nnz(goodCOD);

    if opt.CloseExisting
        close(findall(groot,'Type','figure','Tag','CrackPathFinalPlot'));
    end

    figDir=char(opt.FigureDir);
    if isempty(strtrim(figDir)),figDir=fullfile(runDir,'plots');end
    if ~local_is_absolute_path(figDir),figDir=fullfile(root,figDir);end
    formats=local_formats(opt.Formats);
    if opt.SaveFigures && exist(figDir,'dir')~=7,mkdir(figDir);end

    F=struct();

    % ==================================================================
    % 1. Crack trajectory
    % ==================================================================
    F.trajectory=local_new_figure(opt.Visible,'Crack trajectory');
    ax=axes(F.trajectory); hold(ax,'on'); box(ax,'on');

    C=R0.C;
    local_plot_plate_and_holes(ax,C);
    plot(ax,1e3*vertices(:,1),1e3*vertices(:,2),'-o', ...
        'LineWidth',1.5,'MarkerSize',4,'DisplayName','accepted crack path');

    % Frozen initiation point and three characteristic states.
    plot(ax,1e3*vertices(1,1),1e3*vertices(1,2),'s', ...
        'MarkerSize',7,'LineWidth',1.1,'DisplayName','crack mouth');
    local_plot_landmark_points(ax,L,true);

    xlabel(ax,'x [mm]');
    ylabel(ax,'y [mm]');
    title(ax,sprintf('Accepted crack trajectory through P_{%d}',kDone));
    axis(ax,'equal'); grid(ax,'on');

    % Local geometry view: full hole, accepted path and right boundary.
    [xmin,xmax,ymin,ymax]=local_trajectory_limits(C,vertices);
    xlim(ax,[xmin xmax]); ylim(ax,[ymin ymax]);
    legend(ax,'Location','best');
    hold(ax,'off');

    % ==================================================================
    % 2. Absolute crack angle
    % ==================================================================
    F.theta=local_new_figure(opt.Visible,'Crack angle');
    ax=axes(F.theta);
    plot(ax,T.crack_length_mm,T.theta_deg,'-o','LineWidth',1.5,'MarkerSize',4);
    hold(ax,'on'); yline(ax,0,'--','HandleVisibility','off');
    local_landmark_lines(ax,L);
    xlabel(ax,'Crack length [mm]');
    ylabel(ax,'\theta_k [deg]');
    title(ax,'Absolute local crack direction');
    grid(ax,'on'); box(ax,'on'); hold(ax,'off');

    % ==================================================================
    % 3. Incremental MTS turn
    % ==================================================================
    F.deltaTheta=local_new_figure(opt.Visible,'Incremental MTS turn');
    ax=axes(F.deltaTheta);
    plot(ax,T.crack_length_mm,T.delta_theta_next_deg,'-o', ...
        'LineWidth',1.5,'MarkerSize',4);
    hold(ax,'on'); yline(ax,0,'--','HandleVisibility','off');
    local_landmark_lines(ax,L);
    xlabel(ax,'Crack length [mm]');
    ylabel(ax,'\Delta\theta_{k+1} [deg]');
    title(ax,'Incremental MTS turning angle');
    grid(ax,'on'); box(ax,'on'); hold(ax,'off');

    % ==================================================================
    % 4. Mode-I SIF
    % ==================================================================
    F.KI=local_new_figure(opt.Visible,'Mode-I SIF');
    ax=axes(F.KI);
    plot(ax,T.crack_length_mm,T.KI_unit,'-o','LineWidth',1.5,'MarkerSize',4);
    hold(ax,'on'); local_landmark_lines(ax,L);
    xlabel(ax,'Crack length [mm]');
    ylabel(ax,'K_I at unit traction [MPa sqrt(m)]');
    title(ax,'Mode-I stress-intensity factor');
    grid(ax,'on'); box(ax,'on'); hold(ax,'off');

    % ==================================================================
    % 5. Mode mixity
    % ==================================================================
    F.modeMixity=local_new_figure(opt.Visible,'Mode mixity');
    ax=axes(F.modeMixity);
    plot(ax,T.crack_length_mm,T.KII_over_KI,'-o','LineWidth',1.5,'MarkerSize',4);
    hold(ax,'on'); yline(ax,0,'--','HandleVisibility','off');
    local_landmark_lines(ax,L);
    plot(ax,L.aPeak_mm,L.qPeak,'o','MarkerSize',8,'LineWidth',1.5, ...
        'DisplayName',sprintf('P_{%d}: max K_{II}/K_I',L.peakSegment));
    if L.hasZeroCrossing
        plot(ax,L.aLS_mm,0,'d','MarkerSize',8,'LineWidth',1.5, ...
            'DisplayName',sprintf('local symmetry estimate, %.3f mm',L.aLS_mm));
    end
    xlabel(ax,'Crack length [mm]');
    ylabel(ax,'K_{II}/K_I');
    title(ax,'Mode-mixity evolution');
    grid(ax,'on'); box(ax,'on'); legend(ax,'Location','best'); hold(ax,'off');

    % ==================================================================
    % 6. EDI vs COD
    % ==================================================================
    F.CODvsEDI=local_new_figure(opt.Visible,'EDI vs COD');
    tl=tiledlayout(F.CODvsEDI,2,1,'TileSpacing','compact','Padding','compact');

    ax1=nexttile(tl);
    plot(ax1,COD.crack_length_mm,COD.delta_theta_EDI_deg,'-o', ...
        'LineWidth',1.4,'MarkerSize',4,'DisplayName','EDI');
    hold(ax1,'on');
    plot(ax1,COD.crack_length_mm,COD.delta_theta_COD_deg,'-s', ...
        'LineWidth',1.2,'MarkerSize',4, ...
        'DisplayName',sprintf('COD: [%.2f, %.2f] Delta a, degree %d', ...
        opt.CODWindow(1),opt.CODWindow(2),opt.CODDegree));
    yline(ax1,0,'--','HandleVisibility','off');
    local_landmark_lines(ax1,L);
    ylabel(ax1,'\Delta\theta_{k+1} [deg]');
    title(ax1,'MTS direction: EDI and COD');
    grid(ax1,'on'); box(ax1,'on'); legend(ax1,'Location','best'); hold(ax1,'off');

    ax2=nexttile(tl);
    codMinusEDI=COD.delta_theta_COD_deg-COD.delta_theta_EDI_deg;
    plot(ax2,COD.crack_length_mm,codMinusEDI,'-o','LineWidth',1.3,'MarkerSize',4);
    hold(ax2,'on'); yline(ax2,0,'--','HandleVisibility','off');
    local_landmark_lines(ax2,L);
    xlabel(ax2,'Crack length [mm]');
    ylabel(ax2,'COD - EDI [deg]');
    title(ax2,sprintf('RMS difference %.4g deg; max |difference| %.4g deg',codRMS,codMax));
    grid(ax2,'on'); box(ax2,'on'); hold(ax2,'off');

    % ==================================================================
    % 7. Numerical quality
    % ==================================================================
    F.quality=local_new_figure(opt.Visible,'Numerical quality');
    tl=tiledlayout(F.quality,2,2,'TileSpacing','compact','Padding','compact');

    acceptedQ=Mesh.accepted & Mesh.pass;
    unsolvedQ=~Mesh.accepted & Mesh.pass;

    ax1=nexttile(tl);
    plot(ax1,Mesh.crack_length_mm(acceptedQ),Mesh.physical_clearance_mm(acceptedQ), ...
        '-o','LineWidth',1.3,'MarkerSize',4,'DisplayName','accepted');
    hold(ax1,'on');
    if any(unsolvedQ)
        plot(ax1,Mesh.crack_length_mm(unsolvedQ),Mesh.physical_clearance_mm(unsolvedQ), ...
            'x','MarkerSize',8,'LineWidth',1.4,'DisplayName','qualified, unsolved');
    end
    yline(ax1,0.75*1e3*da,'--','r_{core}','HandleVisibility','off');
    local_landmark_lines(ax1,L);
    ylabel(ax1,'Physical clearance [mm]');
    title(ax1,'Tip-to-physical-boundary clearance');
    grid(ax1,'on'); box(ax1,'on'); legend(ax1,'Location','best'); hold(ax1,'off');

    ax2=nexttile(tl);
    plot(ax2,Mesh.crack_length_mm(acceptedQ),Mesh.min_angle_deg(acceptedQ), ...
        '-o','LineWidth',1.3,'MarkerSize',4);
    hold(ax2,'on');
    if any(unsolvedQ)
        plot(ax2,Mesh.crack_length_mm(unsolvedQ),Mesh.min_angle_deg(unsolvedQ), ...
            'x','MarkerSize',8,'LineWidth',1.4);
    end
    yline(ax2,20,'--','20 deg gate','HandleVisibility','off');
    local_landmark_lines(ax2,L);
    ylabel(ax2,'Minimum T3 angle [deg]');
    title(ax2,'Mesh minimum-angle gate');
    grid(ax2,'on'); box(ax2,'on'); hold(ax2,'off');

    ax3=nexttile(tl);
    plot(ax3,Mesh.crack_length_mm(acceptedQ),Mesh.max_neighbor_ratio(acceptedQ), ...
        '-o','LineWidth',1.3,'MarkerSize',4);
    hold(ax3,'on');
    if any(unsolvedQ)
        plot(ax3,Mesh.crack_length_mm(unsolvedQ),Mesh.max_neighbor_ratio(unsolvedQ), ...
            'x','MarkerSize',8,'LineWidth',1.4);
    end
    yline(ax3,1.8,'--','1.8 gate','HandleVisibility','off');
    local_landmark_lines(ax3,L);
    xlabel(ax3,'Crack length [mm]');
    ylabel(ax3,'Maximum neighbor ratio');
    title(ax3,'Adjacent element-size ratio');
    grid(ax3,'on'); box(ax3,'on'); hold(ax3,'off');

    ax4=nexttile(tl);
    validPCG=isfinite(T.PCG_iterations);
    plot(ax4,T.crack_length_mm(validPCG),T.PCG_iterations(validPCG), ...
        '-o','LineWidth',1.3,'MarkerSize',4);
    hold(ax4,'on');
    if nGeom>kDone
        xline(ax4,1e3*nGeom*da,':',sprintf('P_{%d} qualified/unsolved',nGeom), ...
            'LabelVerticalAlignment','middle','HandleVisibility','off');
    end
    local_landmark_lines(ax4,L);
    xlabel(ax4,'Crack length [mm]');
    ylabel(ax4,'PCG iterations');
    title(ax4,'Linear-solver effort for accepted states');
    grid(ax4,'on'); box(ax4,'on'); hold(ax4,'off');

    if opt.SaveFigures
        local_export(F.trajectory,figDir,'01_trajectory',formats);
        local_export(F.theta,figDir,'02_theta',formats);
        local_export(F.deltaTheta,figDir,'03_delta_theta',formats);
        local_export(F.KI,figDir,'04_KI',formats);
        local_export(F.modeMixity,figDir,'05_mode_mixity',formats);
        local_export(F.CODvsEDI,figDir,'06_COD_vs_EDI',formats);
        local_export(F.quality,figDir,'07_numerical_quality',formats);
    end

    Out=struct();
    Out.stepTable=T;
    Out.codTable=COD;
    Out.meshTable=Mesh;
    Out.landmarks=L;
    Out.figures=F;
    Out.acceptedVertices=vertices;
    Out.acceptedPhysicalSegments=kDone;
    Out.geometricSegments=nGeom;
    Out.validation=Validation;
    Out.runDir=runDir;
    Out.figureDir=figDir;

    fprintf('\nFINAL CLEAN-RUN PLOTTING COMPLETE\n');
    fprintf('  accepted physical states : P1 ... P%d\n',kDone);
    fprintf('  represented geometry     : P0 ... P%d\n',nGeom);
    fprintf('  max mode mixity          : P%d at %.3f mm, KII/KI=%+.8g\n', ...
        L.peakSegment,L.aPeak_mm,L.qPeak);
    if L.hasZeroCrossing
        fprintf('  local-symmetry estimate  : %.6f mm between P%d and P%d\n', ...
            L.aLS_mm,L.zeroLeftSegment,L.zeroRightSegment);
    else
        fprintf('  local-symmetry estimate  : no post-peak sign crossing found\n');
    end
    fprintf('  last accepted state      : P%d at %.3f mm\n',L.lastSegment,L.aLast_mm);
    fprintf('  COD-vs-EDI RMS / max     : %.6g / %.6g deg\n',codRMS,codMax);
    if nGeom>kDone
        fprintf('  next qualified geometry  : P%d is excluded from scientific curves\n',nGeom);
    end
    if opt.SaveFigures
        fprintf('  figures saved in         : %s\n',figDir);
    end
end


% =========================================================================
function R0=local_load_frozen_state(root,Rin,stateFile)
    if ~isempty(Rin)
        R0=Rin;
        return
    end
    if isempty(strtrim(stateFile))
        stateFile=fullfile(root,'verification','crack_path','stage1_starting_state.mat');
    elseif ~local_is_absolute_path(stateFile)
        stateFile=fullfile(root,stateFile);
    end
    if exist(stateFile,'file')~=2
        error('cleanplot:MissingFrozenState', ...
            'Frozen Stage-I state not found; pass ''FrozenState'',R0.');
    end
    d=load(stateFile,'R0');
    if ~isfield(d,'R0')||~isstruct(d.R0)
        error('cleanplot:BadFrozenState','StateFile must contain struct R0.');
    end
    R0=d.R0;
end


function COD=local_collect_cod(runDir,T,window,degree)
    n=height(T);
    data=nan(n,10);
    for i=1:n
        k=T.segment(i);
        data(i,1:4)=[k,T.crack_length_mm(i),T.delta_theta_next_deg(i), ...
            T.KII_over_KI(i)];
        if k<2,continue,end

        f=fullfile(runDir,sprintf('step_%03d_physical_small.mat',k));
        if exist(f,'file')~=2,continue,end
        d=load(f,'R');
        if ~isfield(d,'R')||~isstruct(d.R)||~isfield(d.R,'fitTable'),continue,end
        F=d.R.fitTable;
        if ~istable(F),continue,end

        id=abs(F.lower_r_over_DeltaA-window(1))<=1e-12 & ...
           abs(F.upper_r_over_DeltaA-window(2))<=1e-12 & ...
           F.degree==degree;
        jj=find(id,1);
        if isempty(jj),continue,end

        data(i,5:10)=[F.KI_COD(jj),F.KII_COD(jj),F.ratio_COD(jj), ...
            F.RMSE_KII(jj),F.median_pointwise_ratio(jj), ...
            F.delta_theta_next_MTS_deg(jj)];
    end

    COD=array2table(data,'VariableNames',{ ...
        'segment','crack_length_mm','delta_theta_EDI_deg','ratio_EDI', ...
        'KI_COD','KII_COD','ratio_COD','RMSE_KII', ...
        'median_pointwise_ratio','delta_theta_COD_deg'});
end


function Mesh=local_collect_mesh(runDir,nGeom,kDone,da)
    data=nan(nGeom,11);
    for k=1:nGeom
        data(k,1)=k;
        data(k,2)=1e3*k*da;
        data(k,3)=double(k<=kDone);

        if k<2,continue,end
        f=fullfile(runDir,sprintf('step_%03d_qualification_small.mat',k));
        if exist(f,'file')~=2,continue,end
        d=load(f,'Small');
        if ~isfield(d,'Small')||~isstruct(d.Small)|| ...
                ~isfield(d.Small,'summary')||~istable(d.Small.summary)|| ...
                height(d.Small.summary)~=1
            continue
        end
        S=d.Small.summary;
        vars=S.Properties.VariableNames;
        need={'physical_clearance_m','min_angle_deg','max_neighbor_ratio', ...
            'T3_elements','T6_nodes','EDI_elements','pass'};
        if ~all(ismember(need,vars)),continue,end

        data(k,4:11)=[1e3*S.physical_clearance_m,S.min_angle_deg, ...
            S.max_neighbor_ratio,S.T3_elements,S.T6_nodes,S.EDI_elements, ...
            double(S.pass),double(exist(fullfile(runDir, ...
            sprintf('step_%03d_physical_small.mat',k)),'file')==2)];
    end

    Mesh=array2table(data,'VariableNames',{ ...
        'segment','crack_length_mm','accepted', ...
        'physical_clearance_mm','min_angle_deg','max_neighbor_ratio', ...
        'T3_elements','T6_nodes','EDI_elements','pass','physical_result_present'});
    Mesh.accepted=logical(Mesh.accepted);
    Mesh.pass=logical(Mesh.pass);
    Mesh.physical_result_present=logical(Mesh.physical_result_present);
end


function L=local_landmarks(T,vertices)
    q=T.KII_over_KI;
    a=T.crack_length_mm;

    finiteQ=isfinite(q);
    if ~any(finiteQ)
        error('cleanplot:NoModeMixity','No finite KII/KI values are available.');
    end

    idx=find(finiteQ);
    [qPeak,ii]=max(q(idx));
    iPeak=idx(ii);

    L=struct();
    L.peakIndex=iPeak;
    L.peakSegment=T.segment(iPeak);
    L.aPeak_mm=a(iPeak);
    L.qPeak=qPeak;
    L.xPeak_mm=1e3*T.tip_x_m(iPeak);
    L.yPeak_mm=1e3*T.tip_y_m(iPeak);

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

        % vertices row 1 is P0, so P_k is row k+1.
        pL=vertices(T.segment(i)+1,:);
        pR=vertices(T.segment(i+1)+1,:);
        p=pL+f*(pR-pL);
        L.xLS_mm=1e3*p(1);
        L.yLS_mm=1e3*p(2);
        break
    end
end


function fig=local_new_figure(vis,name)
    fig=figure('Visible',char(vis),'Name',name,'Color','w', ...
        'Tag','CrackPathFinalPlot');
end


function local_plot_plate_and_holes(ax,C)
    A=1e3*C.A; B=1e3*C.B;
    plot(ax,[0 A A 0 0],[-B -B B B -B],'k-','LineWidth',1.0, ...
        'DisplayName','plate boundary');

    holes={};
    if isfield(C,'holes')&&~isempty(C.holes)
        holes=C.holes;
        if isstruct(holes),holes=num2cell(holes);end
    elseif isfield(C,'hole')&&~isempty(C.hole)
        holes={C.hole};
    end

    for j=1:numel(holes)
        h=holes{j};
        if ~isstruct(h)||~isfield(h,'type')||~strcmpi(h.type,'circle')
            continue
        end
        ph=linspace(0,2*pi,721);
        x=1e3*(h.center(1)+h.r*cos(ph));
        y=1e3*(h.center(2)+h.r*sin(ph));
        if j==1
            dn='circular hole';
        else
            dn=sprintf('hole %d',j);
        end
        plot(ax,x,y,'k-','LineWidth',1.0,'DisplayName',dn);
    end
end


function [xmin,xmax,ymin,ymax]=local_trajectory_limits(C,vertices)
    xx=1e3*vertices(:,1);
    yy=1e3*vertices(:,2);
    A=1e3*C.A;

    holeX=[];holeY=[];
    holes={};
    if isfield(C,'holes')&&~isempty(C.holes)
        holes=C.holes;if isstruct(holes),holes=num2cell(holes);end
    elseif isfield(C,'hole')&&~isempty(C.hole)
        holes={C.hole};
    end
    for j=1:numel(holes)
        h=holes{j};
        if isstruct(h)&&isfield(h,'type')&&strcmpi(h.type,'circle')
            holeX=[holeX,1e3*(h.center(1)-h.r),1e3*(h.center(1)+h.r)]; %#ok<AGROW>
            holeY=[holeY,1e3*(h.center(2)-h.r),1e3*(h.center(2)+h.r)]; %#ok<AGROW>
        end
    end

    if isempty(holeX),holeX=xx;holeY=yy;end
    margin=5;
    xmin=max(0,min([xx;holeX(:)])-margin);
    xmax=min(A+2,max(A,max(xx)+margin));
    ymin=max(-1e3*C.B,min([yy;holeY(:)])-margin);
    ymax=min( 1e3*C.B,max([yy;holeY(:)])+margin);
end


function local_plot_landmark_points(ax,L,withLast)
    plot(ax,L.xPeak_mm,L.yPeak_mm,'o','MarkerSize',8,'LineWidth',1.5, ...
        'DisplayName',sprintf('P_{%d}: max mode mixity',L.peakSegment));
    if L.hasZeroCrossing
        plot(ax,L.xLS_mm,L.yLS_mm,'d','MarkerSize',8,'LineWidth',1.5, ...
            'DisplayName',sprintf('local symmetry est. %.2f mm',L.aLS_mm));
    end
    if withLast
        plot(ax,L.xLast_mm,L.yLast_mm,'s','MarkerSize',8,'LineWidth',1.5, ...
            'DisplayName',sprintf('P_{%d}: last accepted',L.lastSegment));
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
        error('cleanplot:Format','Unsupported figure format(s): %s',strjoin(bad,', '));
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
