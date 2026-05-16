function files = stage1_export_iteration_outputs(outDir, iter, C, G, S1, B, I)
%STAGE1_EXPORT_ITERATION_OUTPUTS Save mesh and stress-field diagnostics.
%
%   files = stage1_export_iteration_outputs(outDir, iter, C, G, S1, B, I)
%
% Outputs written per iteration:
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
    files.meshPng = fullfile(outDir, [tag, '_mesh.png']);
    files.syyPng = fullfile(outDir, [tag, '_stress_syy.png']);
    files.vonMisesPng = fullfile(outDir, [tag, '_stress_vonmises.png']);
    files.mat = fullfile(outDir, [tag, '_fields.mat']);

    save_iteration_data(files.mat, C, G, S1, B, I);

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
