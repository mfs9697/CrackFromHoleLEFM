function files = stage1_export_iteration_outputs(outDir, iter, C, G, S1, B, I)
%STAGE1_EXPORT_ITERATION_OUTPUTS Save mesh and stress-field diagnostics.
%
%   files = stage1_export_iteration_outputs(outDir, iter, C, G, S1, B, I)
%
% Outputs written per iteration:
%   stage1_iter_XX_hole_contour_stress.png
%   stage1_iter_XX_hole_contour_stress.csv
%   stage1_iter_XX_mesh.png
%   stage1_iter_XX_stress_syy.png
%   stage1_iter_XX_stress_vonmises.png
%   stage1_iter_XX_fields.mat

    if nargin < 7
        error('stage1_export_iteration_outputs:NotEnoughInputs', ...
            'Expected outDir, iter, C, G, S1, B, and I.');
    end

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    tag = sprintf('stage1_iter_%02d', iter);

    files = struct();
    files.holeContourPng = fullfile(outDir, [tag, '_hole_contour_stress.png']);
    files.holeContourCsv = fullfile(outDir, [tag, '_hole_contour_stress.csv']);
    files.meshPng = fullfile(outDir, [tag, '_mesh.png']);
    files.syyPng = fullfile(outDir, [tag, '_stress_syy.png']);
    files.vonMisesPng = fullfile(outDir, [tag, '_stress_vonmises.png']);
    files.mat = fullfile(outDir, [tag, '_fields.mat']);

    save_iteration_data(files.mat, C, G, S1, B, I);
    write_hole_contour_csv(files.holeContourCsv, B);

    fig = figure('Name', sprintf('Stage I iteration %d hole contour stress', iter), ...
        'Color', 'w', 'Visible', 'off', 'Position', [100 100 1200 850]);
    plot_hole_contour_stress(G, B, I);
    exportgraphics(fig, files.holeContourPng, 'Resolution', 200);
    close(fig);

    fig = figure('Name', sprintf('Stage I iteration %d mesh', iter), ...
        'Color', 'w', 'Visible', 'off');
    plot_iteration_mesh(G, I);
    exportgraphics(fig, files.meshPng, 'Resolution', 200);
    close(fig);

    fig = figure('Name', sprintf('Stage I iteration %d sigma yy', iter), ...
        'Color', 'w', 'Visible', 'off');
    plot_stress_field(S1, S1.stress(:,2), '\sigma_{yy}', G, I);
    exportgraphics(fig, files.syyPng, 'Resolution', 200);
    close(fig);

    fig = figure('Name', sprintf('Stage I iteration %d von Mises', iter), ...
        'Color', 'w', 'Visible', 'off');
    plot_stress_field(S1, von_mises_stress(S1.stress), '\sigma_{vM}', G, I);
    exportgraphics(fig, files.vonMisesPng, 'Resolution', 200);
    close(fig);
end


function write_hole_contour_csv(filename, B)
%WRITE_HOLE_CONTOUR_CSV Save sampled hole-contour stresses.

    T = table( ...
        rad2deg(B.phi(:)), ...
        B.x(:,1), B.x(:,2), ...
        B.xq(:,1), B.xq(:,2), ...
        B.sig_tt(:), B.sig_nn(:), B.sig_nt(:), B.sig_tt_eff(:), ...
        'VariableNames', { ...
            'phi_deg', ...
            'x_boundary', 'y_boundary', ...
            'x_query', 'y_query', ...
            'sig_tt', 'sig_nn', 'sig_nt', 'sig_tt_eff'});

    writetable(T, filename);
end


function save_iteration_data(filename, C, G, S1, B, I)
%SAVE_ITERATION_DATA Save compact diagnostics without stiffness matrices.

    meshT3 = struct();
    meshT3.p = G.p;
    meshT3.t = G.t;
    meshT3.outerPoly = G.outerPoly;
    meshT3.holeLoops = G.holeLoops;
    meshT3.edgeSets = G.edgeSets;

    meshT6 = S1.mesh;
    stress = S1.stress;
    displacement = S1.U;
    coord_def = S1.coord_def;
    load = S1.load;
    bc = S1.bc;

    meta = struct();
    if isfield(G, 'meta')
        meta = G.meta;
    end

    save(filename, ...
        'C', 'meshT3', 'meshT6', 'stress', 'displacement', 'coord_def', ...
        'B', 'I', 'load', 'bc', 'meta');
end


function plot_hole_contour_stress(G, B, I)
%PLOT_HOLE_CONTOUR_STRESS Plot contour stress and hole mesh points.

    [meshPhi, meshPts] = hole_mesh_angles(G, B);
    sigAtMesh = periodic_interp(B.phi(:), B.sig_tt_eff(:), meshPhi(:));

    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile([1 2]);
    hold on;
    box on;
    grid on;

    plot(rad2deg(B.phi), B.sig_tt_eff, 'LineWidth', 1.2, ...
        'DisplayName', '\sigma_{\theta\theta}^{eff}');

    if ~isempty(meshPhi)
        plot(rad2deg(meshPhi), sigAtMesh, 'k.', ...
            'MarkerSize', 9, 'DisplayName', 'hole mesh nodes');
    end

    plot(rad2deg(I.phi_star), I.sig_tt_unit, ...
        'ro', 'MarkerSize', 6, 'LineWidth', 1.4, ...
        'DisplayName', 'selected peak');

    plot_refinement_windows(G);

    xlabel('\phi, deg');
    ylabel('\sigma_{\theta\theta}^{eff}');
    title(sprintf('Hole-contour stress (%d boundary mesh nodes)', numel(meshPhi)));
    xlim([0, 360]);
    legend('Location', 'best');

    nexttile;
    hold on;
    axis equal;
    box on;

    scatter(B.x(:,1), B.x(:,2), 12, B.sig_tt_eff, 'filled');
    colorbar;
    colormap(parula);

    if ~isempty(meshPts)
        plot(meshPts(:,1), meshPts(:,2), 'ko', ...
            'MarkerSize', 3.0, 'LineWidth', 0.7, ...
            'DisplayName', 'hole mesh nodes');
    end

    plot(I.x_star(1), I.x_star(2), 'rp', ...
        'MarkerSize', 9, 'LineWidth', 1.4, ...
        'MarkerFaceColor', 'y', 'DisplayName', 'selected peak');

    xlabel('x');
    ylabel('y');
    title('Hole contour nodes');

    nexttile;
    plot_peak_mesh_zoom(G, B, I, meshPts);
end


function plot_iteration_mesh(G, I)
%PLOT_ITERATION_MESH Plot the Stage-I T3 mesh and selected peak point.

    hold on;
    axis equal;
    box on;

    triplot(G.t(:,1:3), G.p(:,1), G.p(:,2), 'Color', [0.72 0.72 0.72]);

    if isfield(G, 'outerPoly') && ~isempty(G.outerPoly)
        plot_closed_loop(G.outerPoly, 'k-', 1.1);
    end

    if isfield(G, 'holeLoops') && ~isempty(G.holeLoops)
        for k = 1:numel(G.holeLoops)
            plot_closed_loop(G.holeLoops{k}, 'r-', 1.3);
        end
    end

    plot(I.x_star(1), I.x_star(2), 'ko', ...
        'MarkerSize', 5, 'LineWidth', 1.3, 'MarkerFaceColor', 'y');

    xlabel('x');
    ylabel('y');
    title(sprintf('Stage I mesh, peak \\phi = %.2f deg', rad2deg(I.phi_star)));
end


function plot_stress_field(S1, nodalValue, labelText, G, I)
%PLOT_STRESS_FIELD Plot a nodal stress scalar on the quadratic mesh.

    coord = S1.mesh.coord;
    tri = split_t6_to_t3(S1.mesh.connect);

    patch('Faces', tri, ...
        'Vertices', coord, ...
        'FaceVertexCData', nodalValue(:), ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'none');

    hold on;
    axis equal;
    box on;
    colorbar;
    colormap(parula);

    if nargin >= 4 && isfield(G, 'outerPoly') && ~isempty(G.outerPoly)
        plot_closed_loop(G.outerPoly, 'k-', 1.0);
    end
    if nargin >= 4 && isfield(G, 'holeLoops') && ~isempty(G.holeLoops)
        for k = 1:numel(G.holeLoops)
            plot_closed_loop(G.holeLoops{k}, 'k-', 1.0);
        end
    end
    if nargin >= 5 && isfield(I, 'x_star')
        plot(I.x_star(1), I.x_star(2), 'ko', ...
            'MarkerSize', 5, 'LineWidth', 1.3, 'MarkerFaceColor', 'y');
    end

    xlabel('x');
    ylabel('y');
    title(['Stage I stress field: ', labelText]);
end


function [phi, pts] = hole_mesh_angles(G, B)
%HOLE_MESH_ANGLES Return sorted polar angles of T3 mesh nodes on the hole.

    phi = [];
    pts = [];

    if ~isfield(G, 'p') || isempty(G.p)
        return;
    end

    ids = hole_boundary_node_ids_from_loops(G);

    if isempty(ids) && isfield(G, 'edgeSets') && isfield(G.edgeSets, 'hole')
        ids = unique(G.edgeSets.hole(:));
    end

    ids = ids(ids >= 1 & ids <= size(G.p,1));
    if isempty(ids)
        return;
    end

    c = B.hole.center(:).';
    pts = G.p(ids,:);
    phi = mod(atan2(pts(:,2) - c(2), pts(:,1) - c(1)), 2*pi);

    [phi, order] = sort(phi);
    pts = pts(order,:);
end


function ids = hole_boundary_node_ids_from_loops(G)
%HOLE_BOUNDARY_NODE_IDS_FROM_LOOPS Find nodes on polygonal hole edges.

    ids = [];

    if ~isfield(G, 'holeLoops') || isempty(G.holeLoops)
        return;
    end

    p = G.p;

    Hscale = 1.0;
    if isfield(G, 'meta')
        Hscale = max([Hscale, abs(G.meta.A), abs(G.meta.B)]);
        if isfield(G.meta, 'Hhole') && ~isempty(G.meta.Hhole)
            Hscale = max(Hscale, G.meta.Hhole);
        end
    end
    tol = 1e-8 * Hscale;

    keep = false(size(p,1), 1);

    for ih = 1:numel(G.holeLoops)
        H = G.holeLoops{ih};
        if isempty(H)
            continue;
        end
        if norm(H(end,:) - H(1,:), inf) ~= 0
            Hc = [H; H(1,:)];
        else
            Hc = H;
        end

        for k = 1:size(Hc,1)-1
            d = point_segment_distance(p, Hc(k,:), Hc(k+1,:));
            keep = keep | (d <= tol);
        end
    end

    ids = find(keep);
end


function d = point_segment_distance(P, A, B)
%POINT_SEGMENT_DISTANCE Distance from point rows P to segment AB.

    AB = B - A;
    L2 = max(dot(AB, AB), 1e-30);

    t = ((P - A) * AB.') / L2;
    t = max(0, min(1, t));

    Q = A + t .* AB;
    d = vecnorm(P - Q, 2, 2);
end


function plot_refinement_windows(G)
%PLOT_REFINEMENT_WINDOWS Overlay active angular refinement windows.

    if ~isfield(G, 'meta') || ~isfield(G.meta, 'refinement')
        return;
    end
    R = G.meta.refinement;
    if ~isfield(R, 'active') || ~R.active
        return;
    end
    if ~isfield(R, 'refineHoleAngles') || isempty(R.refineHoleAngles) || ...
            ~isfield(R, 'refineAngleHalfWidth') || isempty(R.refineAngleHalfWidth)
        return;
    end

    angles = mod(R.refineHoleAngles(:), 2*pi);
    widths = R.refineAngleHalfWidth(:);
    if isscalar(widths)
        widths = repmat(widths, size(angles));
    end

    yl = ylim;
    for k = 1:numel(angles)
        cdeg = rad2deg(angles(k));
        wdeg = rad2deg(widths(k));
        xline(cdeg, ':', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off');
        xline(mod(cdeg - wdeg, 360), '--', 'Color', [0.65 0.65 0.65], ...
            'HandleVisibility', 'off');
        xline(mod(cdeg + wdeg, 360), '--', 'Color', [0.65 0.65 0.65], ...
            'HandleVisibility', 'off');
    end
    ylim(yl);
end


function plot_peak_mesh_zoom(G, B, I, meshPts)
%PLOT_PEAK_MESH_ZOOM Show mesh refinement near the selected concentrator.

    hold on;
    axis equal;
    box on;

    triplot(G.t(:,1:3), G.p(:,1), G.p(:,2), 'Color', [0.70 0.70 0.70]);

    if isfield(G, 'holeLoops') && ~isempty(G.holeLoops)
        for k = 1:numel(G.holeLoops)
            plot_closed_loop(G.holeLoops{k}, 'r-', 1.2);
        end
    end

    if ~isempty(meshPts)
        plot(meshPts(:,1), meshPts(:,2), 'ko', ...
            'MarkerSize', 3.0, 'LineWidth', 0.7);
    end

    plot(I.x_star(1), I.x_star(2), 'rp', ...
        'MarkerSize', 9, 'LineWidth', 1.4, 'MarkerFaceColor', 'y');

    R = B.hole.r;
    span = 1.35 * R;
    xlim(I.x_star(1) + span*[-1 1]);
    ylim(I.x_star(2) + span*[-1 1]);

    xlabel('x');
    ylabel('y');
    title('Peak mesh zoom');
end


function yq = periodic_interp(phi, y, phiq)
%PERIODIC_INTERP Linear interpolation of 2*pi-periodic contour data.

    if isempty(phiq)
        yq = [];
        return;
    end

    phi = mod(phi(:), 2*pi);
    y = y(:);
    phiq = mod(phiq(:), 2*pi);

    [phi, order] = sort(phi);
    y = y(order);

    phiExt = [phi; phi(1) + 2*pi];
    yExt = [y; y(1)];

    yq = interp1(phiExt, yExt, phiq, 'linear');
end


function tri = split_t6_to_t3(connect)
%SPLIT_T6_TO_T3 Split each T6 element into four T3 sub-triangles.

    if size(connect,2) == 3
        tri = connect(:,1:3);
        return;
    end

    if size(connect,2) < 6
        error('stage1_export_iteration_outputs:BadConnectivity', ...
            'Expected T3 or T6 connectivity.');
    end

    tri = [
        connect(:,[1 4 6]);
        connect(:,[4 2 5]);
        connect(:,[6 5 3]);
        connect(:,[4 5 6])];
end


function vm = von_mises_stress(stress)
%VON_MISES_STRESS Plane-stress-style scalar diagnostic from [sxx syy sxy].

    sxx = stress(:,1);
    syy = stress(:,2);
    sxy = stress(:,3);

    vm = sqrt(max(sxx.^2 - sxx.*syy + syy.^2 + 3*sxy.^2, 0));
end


function plot_closed_loop(P, style, lineWidth)
%PLOT_CLOSED_LOOP Plot polygon loop P with repeated first point.

    plot([P(:,1); P(1,1)], [P(:,2); P(1,2)], style, 'LineWidth', lineWidth);
end
