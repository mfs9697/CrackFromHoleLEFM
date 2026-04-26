function out = main_hole_initiation_stage12_geom()
%MAIN_HOLE_INITIATION_STAGE12_GEOM
% Current main driver for the two-stage hole-initiation study.
%
% What it does at the current development stage:
%   Stage I:
%     1) build the plate-with-hole mesh
%     2) solve the elastic problem at unit load
%     3) sample sigma_tt along the hole boundary
%     4) determine the initiation point and initiation load
%
%   Stage II:
%     5) define a short crack at the initiation point for a chosen angle theta
%     6) build the appended-hole geometry description
%     7) mesh the appended-hole geometry as a T3 mesh with local tip refinement
%
% What it does NOT yet do:
%     - collapse the two pencil faces to the crack midline
%     - upgrade the collapsed mesh to quadratic
%     - solve the Stage-II cracked problem
%     - compute KI, KII, or the J-integral on the final cracked mesh
%
% Output:
%   out struct with fields:
%       .C, .G, .S1, .B, .I, .G2, .D, .M

    clc;

    %% ============================================================
    % 0. Configuration
    %% ============================================================
    C = cfg_hole_initiation();

    %% ============================================================
    % 1. Stage I: geometry
    %% ============================================================
    fprintf('\n=== Stage I: hole-only geometry ===\n');
    G = geom_hole_only(C);

    %% ============================================================
    % 2. Stage I: solve at unit load
    %% ============================================================
    fprintf('\n=== Stage I: elastic solution at unit load ===\n');
    S1 = solve_hole_only(C, G, 'lambda', 1.0);

    %% ============================================================
    % 3. Stage I: boundary stress sampling
    %% ============================================================
    fprintf('\n=== Stage I: sample hole-boundary stress ===\n');
    B = sample_hole_boundary_stress(C, G, S1);

    %% ============================================================
    % 4. Stage I: initiation point
    %% ============================================================
    fprintf('\n=== Stage I: initiation point ===\n');
    I = find_hole_initiation_point(C, B);

    fprintf('Initiation point:\n');
    fprintf('  phi_*          = %.6f rad = %.3f deg\n', I.phi_star, I.phi_star*180/pi);
    fprintf('  x_*            = [%.8f, %.8f]\n', I.x_star(1), I.x_star(2));
    fprintf('  sigma_tt,max   = %.8e (unit load)\n', I.sig_tt_pos_unit);
    fprintf('  lambda_ini     = %.8e\n', I.lambda_ini);
    fprintf('  applied_ini    = %.8e\n', I.sig_applied_ini);

    %% ============================================================
    % 5. Stage-I diagnostic plots
    %% ============================================================
    plot_stage1_boundary_stress(B, I);
    plot_stage1_local_frame(G, I);

    %% ============================================================
    % 6. Stage II: choose a trial short-crack angle
    %% ============================================================
    % Convention:
    %   e_dir = cos(theta)*n_in + sin(theta)*t_hat
    % so theta = 0 means radial inward growth.
    theta = 0.0;

    fprintf('\n=== Stage II: short-crack geometry for theta = %.3f deg ===\n', ...
        theta*180/pi);

    G2 = geom_hole_shortcrack(C, I, theta);

    fprintf('Short crack:\n');
    fprintf('  x0            = [%.8f, %.8f]\n', G2.crack.x0(1),   G2.crack.x0(2));
    fprintf('  xtip          = [%.8f, %.8f]\n', G2.crack.xtip(1), G2.crack.xtip(2));
    fprintf('  a0            = %.8e\n', G2.crack.a0);
    fprintf('  theta         = %.6f rad = %.3f deg\n', ...
        G2.crack.theta, G2.crack.thetaDeg);

    %% ============================================================
    % 7. Stage II: appended-hole geometry description
    %% ============================================================
    fprintf('\n=== Stage II: build appended-hole geometry description ===\n');

    D = build_domain_hole_pencil_polyline( ...
        G2.crack.polyline, ...
        C.A, C.B, C.holes, C.mesh2.chw, ...
        'corner_tol', 1e-10, ...
        'epsMode', 'arclength', ...
        'nArc', 160, ...
        'orientation', 'cw');

    fprintf('Appended-hole geometry description:\n');
    fprintf('  number of holes    = %d\n', numel(D.holeLoops));

    if isfield(D, 'topology') && isfield(D.topology, 'appendedHoleArea') ...
            && ~isnan(D.topology.appendedHoleArea)
        fprintf('  appended-hole area = %.8e\n', D.topology.appendedHoleArea);
    else
        fprintf('  appended-hole area = [not stored]\n');
    end

    fprintf('  start on hole?     = %d\n', D.topology.flags.startOnHole);
    fprintf('  tip inside plate?  = %d\n', D.topology.flags.tipInsidePlate);

    %% ============================================================
    % 8. Plot appended-hole geometry description
    %% ============================================================
    plot_stage2_domain_description(D);

    %% ============================================================
    % 9. Stage II: T3 mesh with explicit local tip refinement
    %% ============================================================
    fprintf('\n=== Stage II: mesh appended-hole geometry ===\n');

    % Current geometry labels for the sharp pencil:
    %   tip vertex : V3
    %   pencil edges: E3 and E161
    %
    % These were read from the PDE geometry plot. If the appended-hole
    % geometry changes, these labels may need updating.
    %v_tip = 3;
    %e_tip = [3 161];
    
    M = mesh_hole_pencil_domain(D, ...
        'Hmin', C.mesh1.hmin, ...
        'Hmax', C.mesh1.hmax, ...
        'Hgrad', C.mesh1.hgrad, ...
        'Hedge', {e_tip, 0.8*C.mesh1.hmin}, ...
        'Hvertex', {[v_tip], 0.1*C.mesh1.hmin}, ...
        'PlotGeom', true, ...
        'PlotMesh', true);

    geomIDs = M.region.geomIDs;
    v_tip = geomIDs.v_tip;
    e_tip = geomIDs.e_tip;

    fprintf('Stage-II mesh summary:\n');
    fprintf('  number of nodes    = %d\n', size(M.p,1));
    fprintf('  number of elements = %d\n', size(M.t,1));

    %% ============================================================
    % 10. Return all current objects
    %% ============================================================
    out = struct();
    out.C  = C;
    out.G  = G;
    out.S1 = S1;
    out.B  = B;
    out.I  = I;
    out.G2 = G2;
    out.D  = D;
    out.M  = M;

    Mc = collapse_pencil_faces_to_midline(M, D, ...
        'EdgeIDs', [3 161], ...
        'TipVertexID', 3);

    out.Mc = Mc;

    plot_collapsed_pencil_mesh(Mc, 'ShowOriginalFaces', true);

    fprintf('\n=== Current driver finished successfully ===\n');
end


% ========================================================================
% local plotting helpers
% ========================================================================

function plot_stage1_boundary_stress(B, I)
% Plot sigma_tt along the hole boundary and mark the selected initiation point.

    figure('Name', 'Stage I: boundary stress', 'Color', 'w'); clf
    plot(B.phi*180/pi, B.sig_tt_eff, 'LineWidth', 1.2); hold on
    plot(I.phi_star*180/pi, I.sig_tt_unit, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5)

    xlabel('\phi, deg')
    ylabel('\sigma_{tt}')
    grid on
    box on
    xlim([0, 360])

    title('Stage I: tangential stress on the hole boundary')
end


function plot_stage1_local_frame(G, I)
% Plot the hole geometry and the selected local frame at the initiation point.

    figure('Name', 'Stage I: initiation point', 'Color', 'w'); clf
    hold on; axis equal; box on

    if isfield(G, 't') && ~isempty(G.t) && isfield(G, 'p') && ~isempty(G.p)
        triplot(G.t(:,1:3), G.p(:,1), G.p(:,2), 'Color', [0.80 0.80 0.80]);
    end

    if isfield(G, 'outerPoly') && ~isempty(G.outerPoly)
        P = G.outerPoly;
        plot([P(:,1); P(1,1)], [P(:,2); P(1,2)], 'k-', 'LineWidth', 1.0);
    end

    if isfield(G, 'holeLoops') && ~isempty(G.holeLoops)
        for k = 1:numel(G.holeLoops)
            H = G.holeLoops{k};
            plot([H(:,1); H(1,1)], [H(:,2); H(1,2)], 'r-', 'LineWidth', 1.2);
        end
    end

    x0    = I.x_star;
    nmat  = I.n_mat_star;    % into the solid
    nhole = I.n_hole_star;   % into the hole
    that  = I.t_hat_star;

    scl = 0.008;

    plot(x0(1), x0(2), 'ko', 'MarkerSize', 6, 'LineWidth', 1.2);
    quiver(x0(1), x0(2), scl*nmat(1),  scl*nmat(2),  0, ...
        'Color', [0 0.5 0], 'LineWidth', 1.4, 'MaxHeadSize', 0.8);
    quiver(x0(1), x0(2), scl*nhole(1), scl*nhole(2), 0, ...
        'Color', [0.2 0.2 0.2], 'LineWidth', 1.4, 'MaxHeadSize', 0.8);
    quiver(x0(1), x0(2), scl*that(1),  scl*that(2),  0, ...
        'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, 'MaxHeadSize', 0.8);

    xlabel('x')
    ylabel('y')
    title('Stage I: initiation point and local frame')
end


function plot_stage2_domain_description(D)
% Plot rectangle, hole loop(s), optional legacy channel polygon, and crack midline.

    figure('Name', 'Stage II: domain description', 'Color', 'w'); clf
    hold on; axis equal; box on

    % outer rectangle
    P = D.outerPoly;
    plot([P(:,1); P(1,1)], [P(:,2); P(1,2)], 'k-', 'LineWidth', 1.2);

    % holes / appended-hole loop
    for k = 1:numel(D.holeLoops)
        H = D.holeLoops{k};
        plot([H(:,1); H(1,1)], [H(:,2); H(1,2)], 'r-', 'LineWidth', 1.3);
    end

    % legacy channel polygon if present
    if isfield(D, 'channelPoly') && ~isempty(D.channelPoly)
        Cc = D.channelPoly;
        plot([Cc(:,1); Cc(1,1)], [Cc(:,2); Cc(1,2)], 'b-', 'LineWidth', 1.5);
    end

    % crack midline
    Pm = D.Pmid;
    plot(Pm(:,1), Pm(:,2), 'go-', 'LineWidth', 1.5, 'MarkerSize', 5);

    xlabel('x')
    ylabel('y')
    title('Stage II: appended-hole geometry description')
end