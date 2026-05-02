function [KI, KII, Dbg] = SIF_LEFM_circle2_debug_swapPQ(mesh, U, V, mat, rI, varargin)
%SIF_LEFM_CIRCLE2_DEBUG_SWAPPQ
% Diagnostic variant of SIF_LEFM_circle2_debug.
%
% This version keeps the same circular contour, same P/Q point locations,
% same element lookup, and same field evaluation as SIF_LEFM_circle2_debug,
% but swaps the P/Q convention in the mode-separation formulas.
%
% Purpose:
%   To test whether the sign/parity behavior of extracted KII is caused by
%   the upper/lower-side convention in the P/Q mirrored-point separation.
%
% Usage:
%   [KI, KII, Dbg] = SIF_LEFM_circle2_debug_swapPQ(mesh, U, V, mat, rI)
%   [KI, KII, Dbg] = SIF_LEFM_circle2_debug_swapPQ(..., ...
%       'nthet', 200, 'plot', true, 'verbose', true)
%
% Inputs:
%   mesh   : struct with fields
%              .coord    [nT6 x 2]
%              .connect  [nElem x 6]
%              .coord3   [nT3 x 2]
%              .connect3 [nElem x 3]
%   U      : displacement vector [ux1; uy1; ux2; uy2; ...]
%   V      : crack polyline, last point is crack tip
%   mat    : struct with .E, .nu, .Dmat or .D
%   rI     : circular contour radius. If empty, uses 0.7*crack length.
%
% Outputs:
%   KI, KII : stress intensity factors
%   Dbg     : diagnostic struct

    %% ------------------------------------------------------------
    % 0. Parse options
    %% ------------------------------------------------------------
    ip = inputParser;
    addParameter(ip, 'nthet', 100, @(x)isnumeric(x) && isscalar(x) && x >= 10);
    addParameter(ip, 'eps_th', 1e-3, @(x)isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'plot', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'verbose', false, @(x)islogical(x) || isnumeric(x));
    parse(ip, varargin{:});

    nthet   = round(ip.Results.nthet);
    eps_th  = ip.Results.eps_th;
    doPlot  = logical(ip.Results.plot);
    verbose = logical(ip.Results.verbose);

    %% ------------------------------------------------------------
    % 1. Get mesh and material
    %% ------------------------------------------------------------
    must_have(mesh, 'coord');
    must_have(mesh, 'connect');
    must_have(mesh, 'coord3');
    must_have(mesh, 'connect3');

    must_have(mat, 'E');
    must_have(mat, 'nu');

    if isfield(mat, 'Dmat') && ~isempty(mat.Dmat)
        Dmat = mat.Dmat;
    elseif isfield(mat, 'D') && ~isempty(mat.D)
        Dmat = mat.D;
    else
        error('SIF_LEFM_circle2_debug_swapPQ:MissingDmat', ...
            'Material struct must contain .Dmat or .D.');
    end

    coord    = mesh.coord;
    connect  = mesh.connect;
    coord3   = mesh.coord3;
    connect3 = mesh.connect3;

    E  = mat.E;
    nu = mat.nu;

    if size(connect,2) ~= 6
        error('SIF_LEFM_circle2_debug_swapPQ:BadT6Connect', ...
            'mesh.connect must have 6 columns for T6 elements.');
    end

    if size(connect3,2) ~= 3
        error('SIF_LEFM_circle2_debug_swapPQ:BadT3Connect', ...
            'mesh.connect3 must have 3 columns for T3 parent elements.');
    end

    if numel(U) ~= 2*size(coord,1)
        error('SIF_LEFM_circle2_debug_swapPQ:BadDisplacementSize', ...
            'numel(U) must be 2*size(mesh.coord,1).');
    end

    %% ------------------------------------------------------------
    % 2. Crack-tip frame and radius
    %% ------------------------------------------------------------
    if size(V,1) < 2 || size(V,2) ~= 2
        error('SIF_LEFM_circle2_debug_swapPQ:BadCrackPolyline', ...
            'V must be an [nPts x 2] crack polyline with nPts >= 2.');
    end

    x_tip = V(end,:).';
    a_tot = norm(V(end,:).' - V(1,:).');

    if nargin < 5 || isempty(rI)
        rI = 0.7 * a_tot;
    end

    if ~(isnumeric(rI) && isscalar(rI) && isfinite(rI) && rI > 0)
        error('SIF_LEFM_circle2_debug_swapPQ:BadRadius', ...
            'rI must be a positive finite scalar.');
    end

    p_tip  = V(end,:).';
    p_prev = V(end-1,:).';

    e1 = p_tip - p_prev;
    Le = norm(e1);

    if Le <= eps
        error('SIF_LEFM_circle2_debug_swapPQ:DegenerateLastSegment', ...
            'The last crack segment has nearly zero length.');
    end

    e1 = e1 / Le;          % local x1, crack direction
    e2 = [-e1(2); e1(1)];  % local x2, crack normal

    R_gl  = [e1, e2];      % x_global = R_gl*x_local
    R_loc = R_gl.';        % x_local  = R_loc*x_global

    %% ------------------------------------------------------------
    % 3. Build contour points
    %% ------------------------------------------------------------
    if eps_th >= pi/2
        error('SIF_LEFM_circle2_debug_swapPQ:BadEpsTheta', ...
            'eps_th must be smaller than pi/2.');
    end

    theta = linspace(eps_th, pi - eps_th, nthet).';

    xP_loc = [ rI*cos(theta),  rI*sin(theta) ];
    xQ_loc = [ rI*cos(theta), -rI*sin(theta) ];

    xP_gl = (R_gl * xP_loc.').';
    xQ_gl = (R_gl * xQ_loc.').';

    xP_gl = xP_gl + x_tip.';
    xQ_gl = xQ_gl + x_tip.';

    %% ------------------------------------------------------------
    % 4. Locate contour points in the T3 parent mesh
    %% ------------------------------------------------------------
    TR = triangulation(connect3, coord3);

    [elemP, baryP] = pointLocation(TR, xP_gl(:,1), xP_gl(:,2));
    [elemQ, baryQ] = pointLocation(TR, xQ_gl(:,1), xQ_gl(:,2));

    outsideP = isnan(elemP);
    outsideQ = isnan(elemQ);

    if any(outsideP) || any(outsideQ)
        Dbg = make_partial_debug();
        Dbg.error = 'Some contour points lie outside the mesh.';
        Dbg.outsideP = outsideP;
        Dbg.outsideQ = outsideQ;

        error('SIF_LEFM_circle2_debug_swapPQ:ContourOutsideMesh', ...
            'Some contour points lie outside the mesh: P=%d, Q=%d.', ...
            nnz(outsideP), nnz(outsideQ));
    end

    %% ------------------------------------------------------------
    % 5. Allocate debug arrays
    %% ------------------------------------------------------------
    J1 = zeros(nthet,1);
    J2 = zeros(nthet,1);

    epsP_all   = zeros(nthet,3);
    epsQ_all   = zeros(nthet,3);
    epsI_all   = zeros(nthet,3);
    epsII_all  = zeros(nthet,3);

    sigP_all   = zeros(nthet,3);
    sigQ_all   = zeros(nthet,3);
    sigI_all   = zeros(nthet,3);
    sigII_all  = zeros(nthet,3);

    duI_all    = zeros(nthet,2);
    duII_all   = zeros(nthet,2);

    UI_all     = zeros(nthet,1);
    UII_all    = zeros(nthet,1);

    GradP_loc_all = zeros(2,2,nthet);
    GradQ_loc_all = zeros(2,2,nthet);

    detJP = nan(nthet,1);
    detJQ = nan(nthet,1);

    baryMinP = min(baryP, [], 2);
    baryMinQ = min(baryQ, [], 2);

    sameElemPQ = elemP == elemQ;

    %% ------------------------------------------------------------
    % 6. Loop over contour points
    %% ------------------------------------------------------------
    for k = 1:nthet

        % ---------- P: nominal upper side ----------
        xiP = baryP(k,1:2).';
        eP = elemP(k);

        nodesP = connect(eP,:);
        XP = coord(nodesP,:);

        [~, detJP(k), dNdxP] = BN_local(xiP, XP);

        uelP = [U(2*nodesP - 1), U(2*nodesP)];

        dux_dxP = dNdxP(1,:) * uelP(:,1);
        dux_dyP = dNdxP(2,:) * uelP(:,1);
        duy_dxP = dNdxP(1,:) * uelP(:,2);
        duy_dyP = dNdxP(2,:) * uelP(:,2);

        GradP_gl = [dux_dxP, dux_dyP;
                    duy_dxP, duy_dyP];

        GradP_loc = R_loc * GradP_gl * R_gl;
        GradP_loc_all(:,:,k) = GradP_loc;

        eps11P = GradP_loc(1,1);
        eps22P = GradP_loc(2,2);
        gam12P = GradP_loc(1,2) + GradP_loc(2,1);

        epsP = [eps11P; eps22P; gam12P];
        sigP = Dmat * epsP;

        du_dx1_P = GradP_loc(:,1);

        % ---------- Q: nominal lower side ----------
        xiQ = baryQ(k,1:2).';
        eQ = elemQ(k);

        nodesQ = connect(eQ,:);
        XQ = coord(nodesQ,:);

        [~, detJQ(k), dNdxQ] = BN_local(xiQ, XQ);

        uelQ = [U(2*nodesQ - 1), U(2*nodesQ)];

        dux_dxQ = dNdxQ(1,:) * uelQ(:,1);
        dux_dyQ = dNdxQ(2,:) * uelQ(:,1);
        duy_dxQ = dNdxQ(1,:) * uelQ(:,2);
        duy_dyQ = dNdxQ(2,:) * uelQ(:,2);

        GradQ_gl = [dux_dxQ, dux_dyQ;
                    duy_dxQ, duy_dyQ];

        GradQ_loc = R_loc * GradQ_gl * R_gl;
        GradQ_loc_all(:,:,k) = GradQ_loc;

        eps11Q = GradQ_loc(1,1);
        eps22Q = GradQ_loc(2,2);
        gam12Q = GradQ_loc(1,2) + GradQ_loc(2,1);

        epsQ = [eps11Q; eps22Q; gam12Q];
        sigQ = Dmat * epsQ;

        du_dx1_Q = GradQ_loc(:,1);

        % ---------------------------------------------------------
        % Swapped P/Q mode separation
        % ---------------------------------------------------------
        % This is the key diagnostic change.
        %
        % Original convention in SIF_LEFM_circle2_debug:
        %   epsI  = 0.5*[P1+Q1; P2+Q2; P3-Q3]
        %   epsII = 0.5*[P1-Q1; P2-Q2; P3+Q3]
        %
        % Swapped convention:
        %   formally exchange P <-> Q in the separated fields.
        %
        % Mode I is unchanged for normal components and changes sign only
        % in the engineering shear separation. Mode II changes sign in the
        % normal-difference terms. This file is diagnostic: compare parity
        % KII(-theta) and KII(+theta) with the original routine.

        epsI = 0.5 * [ ...
            epsQ(1) + epsP(1);
            epsQ(2) + epsP(2);
            epsQ(3) - epsP(3) ];

        epsII = 0.5 * [ ...
            epsQ(1) - epsP(1);
            epsQ(2) - epsP(2);
            epsQ(3) + epsP(3) ];

        sigI = 0.5 * [ ...
            sigQ(1) + sigP(1);
            sigQ(2) + sigP(2);
            sigQ(3) - sigP(3) ];

        sigII = 0.5 * [ ...
            sigQ(1) - sigP(1);
            sigQ(2) - sigP(2);
            sigQ(3) + sigP(3) ];

        duI = 0.5 * [ ...
            du_dx1_Q(1) + du_dx1_P(1);
            du_dx1_Q(2) - du_dx1_P(2) ];

        duII = 0.5 * [ ...
            du_dx1_Q(1) - du_dx1_P(1);
            du_dx1_Q(2) + du_dx1_P(2) ];

        UI  = 0.5 * (sigI.'  * epsI);
        UII = 0.5 * (sigII.' * epsII);

        % ---------------------------------------------------------
        % J-integrands
        % ---------------------------------------------------------
        n1 = cos(theta(k));
        n2 = sin(theta(k));

        tI1 = sigI(1)*n1 + sigI(3)*n2;
        tI2 = sigI(3)*n1 + sigI(2)*n2;

        J1(k) = UI*n1 - (tI1*duI(1) + tI2*duI(2));

        tII1 = sigII(1)*n1 + sigII(3)*n2;
        tII2 = sigII(3)*n1 + sigII(2)*n2;

        J2(k) = UII*n1 - (tII1*duII(1) + tII2*duII(2));

        % ---------------------------------------------------------
        % Store pointwise debug data
        % ---------------------------------------------------------
        epsP_all(k,:)  = epsP.';
        epsQ_all(k,:)  = epsQ.';
        epsI_all(k,:)  = epsI.';
        epsII_all(k,:) = epsII.';

        sigP_all(k,:)  = sigP.';
        sigQ_all(k,:)  = sigQ.';
        sigI_all(k,:)  = sigI.';
        sigII_all(k,:) = sigII.';

        duI_all(k,:)   = duI.';
        duII_all(k,:)  = duII.';

        UI_all(k)      = UI;
        UII_all(k)     = UII;
    end

    %% ------------------------------------------------------------
    % 7. Integrate and convert to SIFs
    %% ------------------------------------------------------------
    JI  = 2 * rI * trapz(theta, J1);
    JII = 2 * rI * trapz(theta, J2);

    Eeff = E/(1 - nu^2);

    KI = sqrt(abs(JI) * Eeff);

    if JII < 0
        KII = -sqrt(abs(JII) * Eeff);
    else
        KII =  sqrt(abs(JII) * Eeff);
    end

    %% ------------------------------------------------------------
    % 8. Build debug output
    %% ------------------------------------------------------------
    Dbg = struct();

    Dbg.variant = 'swapPQ';

    Dbg.KI  = KI;
    Dbg.KII = KII;

    Dbg.JI  = JI;
    Dbg.JII = JII;
    Dbg.Eeff = Eeff;

    Dbg.rI = rI;
    Dbg.a_tot = a_tot;
    Dbg.rI_over_a_tot = rI / a_tot;

    Dbg.nthet = nthet;
    Dbg.eps_th = eps_th;

    Dbg.theta = theta;
    Dbg.thetaDeg = rad2deg(theta);

    Dbg.J1 = J1;
    Dbg.J2 = J2;
    Dbg.J1_absint = 2 * rI * trapz(theta, abs(J1));
    Dbg.J2_absint = 2 * rI * trapz(theta, abs(J2));

    Dbg.JI_over_absint  = safe_divide(JI,  Dbg.J1_absint);
    Dbg.JII_over_absint = safe_divide(JII, Dbg.J2_absint);

    Dbg.x_tip = x_tip.';
    Dbg.e1 = e1.';
    Dbg.e2 = e2.';
    Dbg.R_gl = R_gl;
    Dbg.R_loc = R_loc;

    Dbg.xP_loc = xP_loc;
    Dbg.xQ_loc = xQ_loc;
    Dbg.xP_gl = xP_gl;
    Dbg.xQ_gl = xQ_gl;

    Dbg.elemP = elemP;
    Dbg.elemQ = elemQ;
    Dbg.baryP = baryP;
    Dbg.baryQ = baryQ;

    Dbg.baryMinP = baryMinP;
    Dbg.baryMinQ = baryMinQ;
    Dbg.sameElemPQ = sameElemPQ;

    Dbg.detJP = detJP;
    Dbg.detJQ = detJQ;

    Dbg.epsP = epsP_all;
    Dbg.epsQ = epsQ_all;
    Dbg.epsI = epsI_all;
    Dbg.epsII = epsII_all;

    Dbg.sigP = sigP_all;
    Dbg.sigQ = sigQ_all;
    Dbg.sigI = sigI_all;
    Dbg.sigII = sigII_all;

    Dbg.duI = duI_all;
    Dbg.duII = duII_all;

    Dbg.UI = UI_all;
    Dbg.UII = UII_all;

    Dbg.GradP_loc = GradP_loc_all;
    Dbg.GradQ_loc = GradQ_loc_all;

    Dbg.diagnostics = local_diagnostics(Dbg);

    if verbose
        print_debug_summary(Dbg);
    end

    if doPlot
        plot_debug_contour(Dbg, coord3, connect3, V);
        plot_debug_integrands(Dbg);
        plot_debug_elements(Dbg);
    end

    %% ------------------------------------------------------------
    % Nested helper: partial debug struct for early errors
    %% ------------------------------------------------------------
    function D = make_partial_debug()
        D = struct();
        D.variant = 'swapPQ';
        D.rI = rI;
        D.a_tot = a_tot;
        D.nthet = nthet;
        D.theta = theta;
        D.x_tip = x_tip.';
        D.e1 = e1.';
        D.e2 = e2.';
        D.xP_gl = xP_gl;
        D.xQ_gl = xQ_gl;
        D.elemP = elemP;
        D.elemQ = elemQ;
        D.baryP = baryP;
        D.baryQ = baryQ;
    end
end


% ========================================================================
% Diagnostics summary
% ========================================================================

function S = local_diagnostics(Dbg)

    S = struct();

    S.KII_over_KI = safe_divide(Dbg.KII, Dbg.KI);

    S.JI = Dbg.JI;
    S.JII = Dbg.JII;

    S.JI_over_absint  = Dbg.JI_over_absint;
    S.JII_over_absint = Dbg.JII_over_absint;

    S.min_baryP = min(Dbg.baryMinP);
    S.min_baryQ = min(Dbg.baryMinQ);

    S.num_near_edge_P_1e6 = nnz(Dbg.baryMinP < 1e-6);
    S.num_near_edge_Q_1e6 = nnz(Dbg.baryMinQ < 1e-6);

    S.num_near_edge_P_1e4 = nnz(Dbg.baryMinP < 1e-4);
    S.num_near_edge_Q_1e4 = nnz(Dbg.baryMinQ < 1e-4);

    S.frac_same_elem_PQ = mean(Dbg.sameElemPQ);

    S.min_detJP = min(Dbg.detJP);
    S.min_detJQ = min(Dbg.detJQ);
    S.max_detJP = max(Dbg.detJP);
    S.max_detJQ = max(Dbg.detJQ);

    S.J2_sign_changes = nnz(diff(sign_nozero(Dbg.J2)) ~= 0);
    S.J1_sign_changes = nnz(diff(sign_nozero(Dbg.J1)) ~= 0);

    S.max_abs_J1 = max(abs(Dbg.J1));
    S.max_abs_J2 = max(abs(Dbg.J2));

    S.mean_abs_J1 = mean(abs(Dbg.J1));
    S.mean_abs_J2 = mean(abs(Dbg.J2));
end


% ========================================================================
% Console summary
% ========================================================================

function print_debug_summary(Dbg)

    S = Dbg.diagnostics;

    fprintf('\nSIF_LEFM_circle2_debug_swapPQ summary:\n');
    fprintf('  variant            = %s\n', Dbg.variant);
    fprintf('  rI                 = %.8e\n', Dbg.rI);
    fprintf('  rI/a_tot           = %.8e\n', Dbg.rI_over_a_tot);
    fprintf('  KI                 = %.8e\n', Dbg.KI);
    fprintf('  KII                = %.8e\n', Dbg.KII);
    fprintf('  KII/KI             = %.8e\n', S.KII_over_KI);
    fprintf('  JI                 = %.8e\n', Dbg.JI);
    fprintf('  JII                = %.8e\n', Dbg.JII);
    fprintf('  JI / int|J1|       = %.8e\n', S.JI_over_absint);
    fprintf('  JII / int|J2|      = %.8e\n', S.JII_over_absint);
    fprintf('  min bary P/Q       = %.3e / %.3e\n', S.min_baryP, S.min_baryQ);
    fprintf('  near-edge P/Q <1e-6= %d / %d\n', ...
        S.num_near_edge_P_1e6, S.num_near_edge_Q_1e6);
    fprintf('  near-edge P/Q <1e-4= %d / %d\n', ...
        S.num_near_edge_P_1e4, S.num_near_edge_Q_1e4);
    fprintf('  frac same elem P,Q = %.4f\n', S.frac_same_elem_PQ);
    fprintf('  J1 sign changes    = %d\n', S.J1_sign_changes);
    fprintf('  J2 sign changes    = %d\n', S.J2_sign_changes);
end


% ========================================================================
% Plot contour over parent T3 mesh
% ========================================================================

function plot_debug_contour(Dbg, coord3, connect3, V)

    figure('Name', 'SIF debug swapPQ: contour points', 'Color', 'w');
    clf;
    hold on;
    axis equal;
    box on;
    grid on;

    triplot(connect3, coord3(:,1), coord3(:,2), 'Color', [0.80 0.80 0.80]);

    plot(V(:,1), V(:,2), 'k-o', 'LineWidth', 1.5, ...
        'MarkerSize', 4, 'DisplayName', 'crack polyline');

    plot(Dbg.xP_gl(:,1), Dbg.xP_gl(:,2), 'bo-', ...
        'LineWidth', 1.0, 'MarkerSize', 3, 'DisplayName', 'P contour');

    plot(Dbg.xQ_gl(:,1), Dbg.xQ_gl(:,2), 'ro-', ...
        'LineWidth', 1.0, 'MarkerSize', 3, 'DisplayName', 'Q contour');

    plot(Dbg.x_tip(1), Dbg.x_tip(2), 'kp', ...
        'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', 'tip');

    xlabel('x');
    ylabel('y');
    title(sprintf('swapPQ circular contour, r_I = %.4e', Dbg.rI));
    legend('Location', 'best');
end


% ========================================================================
% Plot J-integrands
% ========================================================================

function plot_debug_integrands(Dbg)

    figure('Name', 'SIF debug swapPQ: J integrands', 'Color', 'w');
    clf;

    tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    hold on;
    box on;
    grid on;
    plot(Dbg.thetaDeg, Dbg.J1, 'o-', 'LineWidth', 1.2);
    yline(0, 'k--');
    xlabel('\theta on contour, deg');
    ylabel('J_1 integrand');
    title(sprintf('swapPQ mode-I integrand, J_I = %.4e', Dbg.JI));

    nexttile;
    hold on;
    box on;
    grid on;
    plot(Dbg.thetaDeg, Dbg.J2, 'o-', 'LineWidth', 1.2);
    yline(0, 'k--');
    xlabel('\theta on contour, deg');
    ylabel('J_2 integrand');
    title(sprintf('swapPQ mode-II integrand, J_{II} = %.4e', Dbg.JII));
end


% ========================================================================
% Plot element/barycentric diagnostics
% ========================================================================

function plot_debug_elements(Dbg)

    figure('Name', 'SIF debug swapPQ: element and barycentric diagnostics', ...
        'Color', 'w');
    clf;

    tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    hold on;
    box on;
    grid on;
    plot(Dbg.thetaDeg, Dbg.elemP, 'b.-', 'DisplayName', 'elemP');
    plot(Dbg.thetaDeg, Dbg.elemQ, 'r.-', 'DisplayName', 'elemQ');
    xlabel('\theta on contour, deg');
    ylabel('element id');
    title('Element IDs along contour');
    legend('Location', 'best');

    nexttile;
    hold on;
    box on;
    grid on;
    plot(Dbg.thetaDeg, Dbg.baryMinP, 'b.-', 'DisplayName', 'P');
    plot(Dbg.thetaDeg, Dbg.baryMinQ, 'r.-', 'DisplayName', 'Q');
    yline(1e-6, 'k--', '1e-6');
    yline(1e-4, 'k:',  '1e-4');
    xlabel('\theta on contour, deg');
    ylabel('min barycentric coordinate');
    title('Near-edge indicator');
    legend('Location', 'best');

    nexttile;
    hold on;
    box on;
    grid on;
    plot(Dbg.thetaDeg, Dbg.sameElemPQ, 'ko-');
    xlabel('\theta on contour, deg');
    ylabel('same element P/Q');
    ylim([-0.1, 1.1]);
    title('Whether mirrored P and Q points lie in the same T3 element');
end


% ========================================================================
% Small helpers
% ========================================================================

function must_have(S, field)
    if ~isstruct(S) || ~isfield(S, field) || isempty(S.(field))
        error('SIF_LEFM_circle2_debug_swapPQ:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end


function y = safe_divide(a, b)
    if abs(b) < eps
        y = nan;
    else
        y = a / b;
    end
end


function s = sign_nozero(x)
    s = sign(x);
    for i = 1:numel(s)
        if s(i) == 0
            if i == 1
                j = find(s ~= 0, 1, 'first');
                if isempty(j)
                    s(i) = 0;
                else
                    s(i) = s(j);
                end
            else
                s(i) = s(i-1);
            end
        end
    end
end