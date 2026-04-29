function [G2, D, M, Mc] = build_stage2_cracked_mesh_for_theta(C, I, theta, varargin)
%BUILD_STAGE2_CRACKED_MESH_FOR_THETA
% Build and collapse the appended-hole short-crack mesh for a given angle.
%
% theta is measured from the material-side normal to the hole boundary.

    ip = inputParser;
    addParameter(ip, 'PlotGeom', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'PlotMesh', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'PlotCollapsed', false, @(x)islogical(x) || isnumeric(x));
    parse(ip, varargin{:});

    plotGeom      = logical(ip.Results.PlotGeom);
    plotMesh      = logical(ip.Results.PlotMesh);
    plotCollapsed = logical(ip.Results.PlotCollapsed);

    G2 = geom_hole_shortcrack(C, I, theta);

    D = build_domain_hole_pencil_polyline( ...
        G2.crack.polyline, ...
        C.A, C.B, C.holes, C.mesh2.chw, ...
        'corner_tol', 1e-10, ...
        'epsMode', 'arclength', ...
        'nArc', 160, ...
        'orientation', 'cw');

    if plotGeom
        local_plot_stage2_domain_description(D);
    end

    M = mesh_hole_pencil_domain(D, ...
        'Hmin', C.mesh1.hmin, ...
        'Hmax', C.mesh1.hmax, ...
        'Hgrad', C.mesh1.hgrad, ...
        'PlotGeom', plotGeom, ...
        'PlotMesh', plotMesh);

    geomIDs = M.region.geomIDs;

    Mc = collapse_pencil_faces_to_midline(M, D, ...
        'EdgeIDs', geomIDs.e_tip, ...
        'TipVertexID', geomIDs.v_tip);

    if plotCollapsed
        plot_collapsed_pencil_mesh(Mc, 'ShowOriginalFaces', true);
    end
end


function local_plot_stage2_domain_description(D)

    figure('Name', 'Stage II: appended-hole geometry description', ...
        'Color', 'w');
    clf;
    hold on;
    axis equal;
    box on;

    P = D.outerPoly;
    plot([P(:,1); P(1,1)], [P(:,2); P(1,2)], ...
        'k-', 'LineWidth', 1.2);

    for k = 1:numel(D.holeLoops)
        H = D.holeLoops{k};
        plot([H(:,1); H(1,1)], [H(:,2); H(1,2)], ...
            'r-', 'LineWidth', 1.3);
    end

    if isfield(D, 'channelPoly') && ~isempty(D.channelPoly)
        Cc = D.channelPoly;
        plot([Cc(:,1); Cc(1,1)], [Cc(:,2); Cc(1,2)], ...
            'b-', 'LineWidth', 1.5);
    end

    Pm = D.Pmid;
    plot(Pm(:,1), Pm(:,2), ...
        'go-', 'LineWidth', 1.5, 'MarkerSize', 5);

    xlabel('x');
    ylabel('y');
    title('Stage II: appended-hole geometry description');
end