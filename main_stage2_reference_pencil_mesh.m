function R2 = main_stage2_reference_pencil_mesh(varargin)
%MAIN_STAGE2_REFERENCE_PENCIL_MESH
% First Stage-II coding step for the crack-from-hole LEFM workflow.
%
% Purpose
% -------
% Consume the FROZEN Stage-I starting-state MAT file and build exactly one
% short sharp crack using the existing appended-hole / pencil meshing chain.
% No Stage-I solve and no SIF calculation are performed here.
%
% Default target:
%   theta_1 = 0 deg
%   a0      = the frozen a0_reserved_m value (currently 4 mm)
%
% The function:
%   1) loads verification/crack_path/stage1_starting_state.mat;
%   2) reconstructs the Stage-II geometry inputs from its one-row summary;
%   3) adapts the frozen state to the legacy initiation-struct interface;
%   4) calls build_stage2_cracked_mesh_for_theta;
%   5) checks geometry/topology/collapse gates;
%   6) optionally writes a compact verification MAT file.
%
% Usage
% -----
%   R2 = main_stage2_reference_pencil_mesh();
%
%   R2 = main_stage2_reference_pencil_mesh( ...
%       'Theta1Deg', 0, ...
%       'PlotGeom', true, ...
%       'PlotMesh', true, ...
%       'PlotCollapsed', true);
%
% Name-value options
% ------------------
%   'StateFile'      frozen Stage-I MAT file
%   'Theta1Deg'      trial first-segment angle in the frozen local frame
%   'PlotGeom'       pass through to Stage-II builder
%   'PlotMesh'       pass through to Stage-II builder
%   'PlotCollapsed'  pass through to Stage-II builder
%   'SaveCompact'    save compact verification MAT (default true)
%   'OutputFile'     explicit compact MAT path; default is generated beside
%                    the Stage-I state file
%
% Output
% ------
%   R2.summary       one-row numeric summary table
%   R2.gates         boolean qualification gates
%   R2.diagnostics   tolerances and geometric diagnostics
%   R2.frozenSummary frozen Stage-I summary row used as input
%   R2.C             reconstructed Stage-II configuration
%   R2.I             adapter initiation struct
%   R2.G2,D,M,Mc     existing Stage-II builder outputs
%
% Important
% ---------
% This driver deliberately does NOT call solve_hole_only, Stage-I stress
% recovery, solve_cracked_LEFM, or any SIF routine.

    ip = inputParser;
    addParameter(ip, 'StateFile', ...
        fullfile('verification', 'crack_path', 'stage1_starting_state.mat'), ...
        @(s)ischar(s) || isstring(s));
    addParameter(ip, 'Theta1Deg', 0, ...
        @(x)isnumeric(x) && isscalar(x) && isfinite(x));
    addParameter(ip, 'PlotGeom', false, ...
        @(x)islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'PlotMesh', false, ...
        @(x)islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'PlotCollapsed', false, ...
        @(x)islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'SaveCompact', true, ...
        @(x)islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'OutputFile', '', ...
        @(s)ischar(s) || isstring(s));
    parse(ip, varargin{:});

    stateFile     = char(ip.Results.StateFile);
    theta1Deg     = ip.Results.Theta1Deg;
    plotGeom      = logical(ip.Results.PlotGeom);
    plotMesh      = logical(ip.Results.PlotMesh);
    plotCollapsed = logical(ip.Results.PlotCollapsed);
    saveCompact   = logical(ip.Results.SaveCompact);
    outputFile    = char(ip.Results.OutputFile);

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('CRACK PATH: STAGE-II REFERENCE PENCIL MESH\n');
    fprintf('============================================================\n');
    fprintf('  Frozen Stage-I state: %s\n', stateFile);
    fprintf('  theta_1 = %+ .10f deg\n', theta1Deg);
    fprintf('  NO Stage-I solve and NO SIF calculation are performed here.\n\n');

    if exist(stateFile, 'file') ~= 2
        error('main_stage2_reference_pencil_mesh:MissingStateFile', ...
            ['Frozen Stage-I state file not found:\n  %s\n', ...
             'Run main_stage1_freeze_starting_state first, or pass StateFile.'], ...
            stateFile);
    end

    S = load(stateFile);
    T = local_extract_summary_table(S);

    if height(T) ~= 1
        error('main_stage2_reference_pencil_mesh:BadSummaryHeight', ...
            'Frozen Stage-I summary must contain exactly one row; got %d.', height(T));
    end

    req = { ...
        'A_m','B_m', ...
        'hole_x_m','hole_y_m','hole_R_m','hole_npoly', ...
        'a0_reserved_m', ...
        'hmin_m','hmax_m','hgrad', ...
        'phi_star_deg','x_star_m','y_star_m', ...
        'nmat_x','nmat_y','that_x','that_y'};

    local_require_table_variables(T, req);

    row = T(1,:);

    % ------------------------------------------------------------
    % Reconstruct Stage-II configuration from the frozen summary.
    % Start from the canonical project configuration only to inherit material,
    % loading and currently tested Stage-II pencil controls. All geometry and
    % Stage-I mesh scales are overwritten by the frozen state.
    % ------------------------------------------------------------
    C = cfg_hole_initiation();

    C.A = row.A_m;
    C.B = row.B_m;

    C.hole.type   = 'circle';
    C.hole.center = [row.hole_x_m, row.hole_y_m];
    C.hole.r      = row.hole_R_m;
    C.hole.npoly  = row.hole_npoly;
    C.holes       = {C.hole};

    C.a0 = row.a0_reserved_m;

    C.mesh1.hmin  = row.hmin_m;
    C.mesh1.hmax  = row.hmax_m;
    C.mesh1.hhole = row.hmin_m;
    C.mesh1.hgrad = row.hgrad;

    C.mesh2.hmax   = row.hmax_m;
    C.mesh2.hhole  = row.hmin_m;
    C.mesh2.hcrack = row.hmin_m;
    C.mesh2.hgrad  = row.hgrad;

    % If a future frozen summary records Stage-II-specific pencil controls,
    % prefer those values. Otherwise retain the canonical tested defaults.
    if ismember('chw_m', T.Properties.VariableNames)
        C.mesh2.chw = row.chw_m;
    end
    if ismember('tip_radius_m', T.Properties.VariableNames)
        C.mesh2.tip_radius = row.tip_radius_m;
    end

    % ------------------------------------------------------------
    % Adapter to the existing Stage-II builder interface.
    % ------------------------------------------------------------
    nmat = [row.nmat_x, row.nmat_y];
    that = [row.that_x, row.that_y];

    I = struct();
    I.x_star      = [row.x_star_m, row.y_star_m];
    I.n_mat_star  = nmat;
    I.n_hole_star = -nmat;
    I.t_hat_star  = that;
    I.phi_star    = deg2rad(row.phi_star_deg);

    if ismember('sigma_tt_peak_unit', T.Properties.VariableNames)
        I.sig_tt_unit = row.sigma_tt_peak_unit;
    end
    if ismember('lambda_ini', T.Properties.VariableNames)
        I.lambda_ini = row.lambda_ini;
    end

    theta1 = deg2rad(theta1Deg);

    eExpected = cos(theta1) * nmat + sin(theta1) * that;
    eExpected = eExpected / norm(eExpected);
    xExpected = I.x_star + C.a0 * eExpected;

    fprintf('FROZEN INPUT TO STAGE II\n');
    fprintf('  x_*          = [%.12g, %.12g] m\n', I.x_star(1), I.x_star(2));
    fprintf('  n_mat        = [%.12g, %.12g]\n', nmat(1), nmat(2));
    fprintf('  t_hat        = [%.12g, %.12g]\n', that(1), that(2));
    fprintf('  a0           = %.9f mm\n', 1e3*C.a0);
    fprintf('  expected tip = [%.12g, %.12g] m\n\n', xExpected(1), xExpected(2));

    % ------------------------------------------------------------
    % Existing proven Stage-II pencil chain.
    % ------------------------------------------------------------
    [G2, D, M, Mc] = build_stage2_cracked_mesh_for_theta( ...
        C, I, theta1, ...
        'PlotGeom', plotGeom, ...
        'PlotMesh', plotMesh, ...
        'PlotCollapsed', plotCollapsed);

    % ------------------------------------------------------------
    % Qualification tolerances.
    % ------------------------------------------------------------
    scaleL  = max([1, C.A, 2*C.B, C.a0]);
    tolGeom = 1e-9 * scaleL;
    tolFrame = 1e-10;

    % ------------------------------------------------------------
    % Geometry / topology diagnostics.
    % ------------------------------------------------------------
    mouthErr = norm(G2.crack.x0 - I.x_star);
    tipErr   = norm(G2.crack.xtip - xExpected);

    crackVec = G2.crack.xtip - G2.crack.x0;
    crackLen = norm(crackVec);
    eActual  = crackVec / crackLen;

    lengthErr = abs(crackLen - C.a0);
    dirErr    = norm(eActual - eExpected);

    collapsedMouthErr = norm(Mc.crack.x0 - I.x_star);
    collapsedTipErr   = norm(Mc.crack.xtip - xExpected);

    [maxUpperLineDist, maxLowerLineDist] = ...
        local_face_distance_to_segment(Mc, I.x_star, xExpected);

    sharedFaceNodes = intersect(Mc.crack.upperNodes(:), Mc.crack.lowerNodes(:));
    sharedOnlyAtTip = numel(sharedFaceNodes) == 1 && ...
        sharedFaceNodes(1) == Mc.crack.tipNode;

    [noInversion, minAreaRatio, minAbsArea2] = ...
        local_check_collapse_element_orientation(M.p, Mc.p, Mc.t);

    crackInsideMaterial = local_segment_in_material( ...
        I.x_star, xExpected, C.hole.center, C.hole.r, C.A, C.B, tolGeom);

    localFrameUnit = ...
        abs(norm(nmat) - 1) <= tolFrame && ...
        abs(norm(that) - 1) <= tolFrame;

    localFrameOrthogonal = abs(dot(nmat, that)) <= tolFrame;

    stage1Pass = true;
    if ismember('stage1_pass', T.Properties.VariableNames)
        stage1Pass = logical(row.stage1_pass);
    end

    % ------------------------------------------------------------
    % Hard boolean gates only.
    % ------------------------------------------------------------
    gates = struct();
    gates.frozenStage1Pass          = stage1Pass;
    gates.localFrameUnit            = localFrameUnit;
    gates.localFrameOrthogonal      = localFrameOrthogonal;
    gates.mouthMatchesFrozenPoint   = mouthErr <= tolGeom;
    gates.segmentLengthCorrect      = lengthErr <= tolGeom;
    gates.segmentDirectionCorrect   = dirErr <= tolFrame;
    gates.tipMatchesExpected        = tipErr <= tolGeom;
    gates.crackSegmentInMaterial    = crackInsideMaterial;
    gates.faceNodeSetsPresent       = ...
        ~isempty(Mc.crack.upperNodes) && ~isempty(Mc.crack.lowerNodes);
    gates.faceTopologyDistinct      = sharedOnlyAtTip;
    gates.collapsedMouthMatches     = collapsedMouthErr <= tolGeom;
    gates.collapsedTipMatches       = collapsedTipErr <= tolGeom;
    gates.collapsedUpperOnMidline   = maxUpperLineDist <= tolGeom;
    gates.collapsedLowerOnMidline   = maxLowerLineDist <= tolGeom;
    gates.noElementInversion        = noInversion;

    gateVals = struct2cell(gates);
    pass = all(cellfun(@(x) isequal(logical(x), true), gateVals));

    % ------------------------------------------------------------
    % Human-readable summary.
    % ------------------------------------------------------------
    summary = table( ...
        theta1Deg, C.a0, ...
        I.x_star(1), I.x_star(2), ...
        xExpected(1), xExpected(2), ...
        size(M.p,1), size(M.t,1), ...
        size(Mc.p,1), size(Mc.t,1), ...
        Mc.crack.nUpper, Mc.crack.nLower, ...
        mouthErr, tipErr, lengthErr, dirErr, ...
        maxUpperLineDist, maxLowerLineDist, ...
        minAreaRatio, minAbsArea2, pass, ...
        'VariableNames', { ...
        'theta1_deg','a0_m', ...
        'mouth_x_m','mouth_y_m', ...
        'tip_x_m','tip_y_m', ...
        'T3_nodes_before_collapse','T3_elements_before_collapse', ...
        'T3_nodes_collapsed','T3_elements_collapsed', ...
        'upper_face_nodes','lower_face_nodes', ...
        'mouth_error_m','tip_error_m','length_error_m','direction_error', ...
        'upper_midline_error_m','lower_midline_error_m', ...
        'min_area_ratio','min_abs_twice_area_m2','stage2_geometry_pass'});

    diagnostics = struct();
    diagnostics.tolGeom = tolGeom;
    diagnostics.tolFrame = tolFrame;
    diagnostics.expectedDirection = eExpected;
    diagnostics.actualDirection = eActual;
    diagnostics.expectedTip = xExpected;
    diagnostics.sharedFaceNodes = sharedFaceNodes;
    diagnostics.minAreaRatio = minAreaRatio;
    diagnostics.minAbsArea2 = minAbsArea2;

    fprintf('STAGE-II REFERENCE PENCIL RESULT\n');
    fprintf('  mouth         = [%.12g, %.12g] m\n', G2.crack.x0(1), G2.crack.x0(2));
    fprintf('  tip           = [%.12g, %.12g] m\n', G2.crack.xtip(1), G2.crack.xtip(2));
    fprintf('  segment e1    = [%.12g, %.12g]\n', eActual(1), eActual(2));
    fprintf('  T3 mesh       = %d nodes, %d elements\n', size(Mc.p,1), size(Mc.t,1));
    fprintf('  face nodes    = upper %d, lower %d\n', ...
        Mc.crack.nUpper, Mc.crack.nLower);
    fprintf('  min area ratio after collapse = %.6e\n\n', minAreaRatio);

    fprintf('STAGE-II GEOMETRY GATES\n');
    fns = fieldnames(gates);
    for k = 1:numel(fns)
        fprintf('  %-30s : %d\n', fns{k}, gates.(fns{k}));
    end

    if pass
        fprintf('\nSTAGE-II REFERENCE PENCIL GEOMETRY PASS.\n');
        fprintf('  The collapsed short-crack mesh is qualified for the next extraction step.\n');
    else
        fprintf('\nSTAGE-II REFERENCE PENCIL GEOMETRY FAIL.\n');
        fprintf('  Do NOT proceed to the physical cracked solve or SIF extraction.\n');
    end

    % ------------------------------------------------------------
    % Full in-memory result.
    % ------------------------------------------------------------
    R2 = struct();
    R2.summary = summary;
    R2.gates = gates;
    R2.diagnostics = diagnostics;
    R2.frozenSummary = row;
    R2.C = C;
    R2.I = I;
    R2.G2 = G2;
    R2.D = D;
    R2.M = M;
    R2.Mc = Mc;

    % ------------------------------------------------------------
    % Compact verification checkpoint: exclude PDE Toolbox handle objects.
    % ------------------------------------------------------------
    if saveCompact
        if isempty(outputFile)
            outDir = fileparts(stateFile);
            if isempty(outDir)
                outDir = '.';
            end
            if exist(outDir, 'dir') ~= 7
                mkdir(outDir);
            end
            tag = local_angle_tag(theta1Deg);
            outputFile = fullfile(outDir, ...
                ['stage2_reference_pencil_', tag, '.mat']);
        else
            outDir = fileparts(outputFile);
            if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
                mkdir(outDir);
            end
        end

        Rsave = struct();
        Rsave.summary = summary;
        Rsave.gates = gates;
        Rsave.diagnostics = diagnostics;
        Rsave.frozenSummary = row;
        Rsave.C = C;
        Rsave.I = I;
        Rsave.G2_crack = G2.crack;
        Rsave.G2_tip = G2.tip;
        Rsave.D_Pmid = D.Pmid;
        Rsave.p_precollapse = M.p;
        Rsave.t_precollapse = M.t;
        Rsave.p_collapsed = Mc.p;
        Rsave.t_collapsed = Mc.t;
        Rsave.crack = Mc.crack;

        save(outputFile, 'Rsave', '-v7');
        R2.outputFile = outputFile;
        fprintf('  Compact MAT: %s\n', outputFile);
    else
        R2.outputFile = '';
    end
end


% =========================================================================
% Local helpers
% =========================================================================

function T = local_extract_summary_table(S)
% Locate the frozen Stage-I one-row summary without depending on the MAT
% variable name used by the freeze driver.

    if isfield(S, 'R0') && isstruct(S.R0) && ...
            isfield(S.R0, 'summary') && istable(S.R0.summary)
        T = S.R0.summary;
        return;
    end

    if isfield(S, 'summary') && istable(S.summary)
        T = S.summary;
        return;
    end

    fn = fieldnames(S);
    for k = 1:numel(fn)
        v = S.(fn{k});
        if isstruct(v) && isfield(v, 'summary') && istable(v.summary)
            T = v.summary;
            return;
        end
    end

    error('main_stage2_reference_pencil_mesh:NoFrozenSummary', ...
        ['Could not locate the frozen Stage-I summary table in the MAT file. ', ...
         'Expected R0.summary, summary, or a struct containing .summary.']);
end


function local_require_table_variables(T, names)
    missing = names(~ismember(names, T.Properties.VariableNames));
    if ~isempty(missing)
        error('main_stage2_reference_pencil_mesh:MissingFrozenFields', ...
            'Frozen Stage-I summary is missing required variables: %s', ...
            strjoin(missing, ', '));
    end
end


function [dU, dL] = local_face_distance_to_segment(Mc, A, B)
    U = Mc.p(Mc.crack.upperNodes, :);
    L = Mc.p(Mc.crack.lowerNodes, :);

    dU = max(local_point_segment_distance(U, A, B));
    dL = max(local_point_segment_distance(L, A, B));
end


function d = local_point_segment_distance(P, A, B)
    AB = B - A;
    L2 = dot(AB, AB);

    if L2 <= 0
        d = vecnorm(P - A, 2, 2);
        return;
    end

    tau = ((P - A) * AB.') / L2;
    tau = max(0, min(1, tau));
    Q = A + tau .* AB;
    d = vecnorm(P - Q, 2, 2);
end


function [ok, minRatio, minAbsA2] = ...
        local_check_collapse_element_orientation(p0, p1, t)

    tri = t(:,1:3);

    A0 = local_twice_signed_area(p0, tri);
    A1 = local_twice_signed_area(p1, tri);

    scaleA = max(max(abs(A0)), 1);
    areaTol = 100 * eps(scaleA);

    valid0 = abs(A0) > areaTol;
    nondegenerate1 = abs(A1) > areaTol;

    sameOrientation = true(size(A0));
    sameOrientation(valid0) = A0(valid0) .* A1(valid0) > 0;

    ok = all(nondegenerate1) && all(sameOrientation);

    ratio = abs(A1) ./ max(abs(A0), areaTol);
    minRatio = min(ratio);
    minAbsA2 = min(abs(A1));
end


function A2 = local_twice_signed_area(p, tri)
    P1 = p(tri(:,1),:);
    P2 = p(tri(:,2),:);
    P3 = p(tri(:,3),:);

    A2 = (P2(:,1)-P1(:,1)).*(P3(:,2)-P1(:,2)) - ...
         (P2(:,2)-P1(:,2)).*(P3(:,1)-P1(:,1));
end


function tf = local_segment_in_material(x0, xtip, center, R, A, B, tol)
% Sample the open segment just beyond the mouth and verify that it stays
% inside the plate and outside the circular cavity.

    s = linspace(0, 1, 101).';
    X = x0 + s .* (xtip - x0);

    insideBox = ...
        X(:,1) >= -tol & X(:,1) <= A + tol & ...
        X(:,2) >= -B - tol & X(:,2) <= B + tol;

    r = vecnorm(X - center, 2, 2);
    outsideHole = r >= R - tol;

    tf = all(insideBox) && all(outsideHole);
end


function tag = local_angle_tag(thetaDeg)
    if abs(thetaDeg) < 5e-13
        tag = 'theta0';
        return;
    end

    if thetaDeg > 0
        sgn = 'p';
    else
        sgn = 'm';
    end

    mag = abs(thetaDeg);
    raw = sprintf('%.6f', mag);
    raw = regexprep(raw, '0+$', '');
    raw = regexprep(raw, '\.$', '');
    raw = strrep(raw, '.', 'p');

    tag = ['theta_', sgn, raw];
end
