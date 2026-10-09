function [G2, D, M, Mc] = build_stage2_cracked_mesh_for_theta(C, I, theta, varargin)
%BUILD_STAGE2_CRACKED_MESH_FOR_THETA
% Build and collapse the appended-hole short-crack mesh for a given angle.
%
% theta is measured from the material-side normal to the hole boundary.
%
% Optional NArc (default 160): number of retained circular-hole arc points
% in the temporary appended-hole polygon.
%
% Optional A0Override: explicit short-crack length for the temporary Stage-II
% carrier. When empty, geom_hole_shortcrack retains its historical C.a0
% default. This keeps existing callers unchanged while allowing controlled
% crack-increment studies to override the reserved Stage-II increment without
% mutating the accepted Stage-I configuration struct.

    ip = inputParser;
    addParameter(ip, 'PlotGeom', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'PlotMesh', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'PlotCollapsed', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'NArc', 160, @(x)isnumeric(x) && isscalar(x) && isfinite(x) && x>=16 && x==round(x));
    addParameter(ip, 'A0Override', [], @(x)isempty(x) || ...
        (isnumeric(x) && isscalar(x) && isfinite(x) && x>0));
    parse(ip, varargin{:});

    plotGeom      = logical(ip.Results.PlotGeom);
    plotMesh      = logical(ip.Results.PlotMesh);
    plotCollapsed = logical(ip.Results.PlotCollapsed);

    if isempty(ip.Results.A0Override)
        G2 = geom_hole_shortcrack(C, I, theta);
    else
        G2 = geom_hole_shortcrack(C, I, theta, ...
            'a0', ip.Results.A0Override);
    end

    D = build_domain_hole_pencil_polyline( ...
        G2.crack.polyline, ...
        C.A, C.B, C.holes, C.mesh2.chw, ...
        'corner_tol', 1e-10, ...
        'epsMode', 'arclength', ...
        'nArc', ip.Results.NArc, ...
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