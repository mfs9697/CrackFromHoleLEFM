function G = geom_hole_only(C)
%GEOM_HOLE_ONLY  Build geometry and mesh for a rectangular plate with hole(s).
%
%   G = geom_hole_only(C)
%
% Required fields in C:
%   C.A              plate length in x-direction, domain x in [0,A]
%   C.B              plate half-height, domain y in [-B,B]
%
% Hole specification:
%   Either
%       C.holes = { struct('type','circle','center',[xc yc],'r',R,'npoly',N), ... }
%   or
%       C.hole  = struct('type','circle','center',[xc yc],'r',R,'npoly',N)
%
% Optional mesh controls:
%   C.mesh1.hmax
%   C.mesh1.hhole
%   C.mesh1.hgrad
%   C.mesh1.refineBox   (currently stored only; not actively used)
%
% Optional plotting:
%   C.plot.show_mesh1
%
% Output:
%   G struct with fields:
%       .p           node coordinates [np x 2]
%       .t           element connectivity [nt x 3]
%       .geom        PDE model
%       .meshobj     PDE mesh object
%       .outerPoly   outer boundary polygon
%       .holeLoops   cell array of hole polygons
%       .hole        first hole spec if exactly one hole exists, else []
%       .edgeSets    boundary node sets
%       .meta        auxiliary metadata
%
% Notes
%   - This is a Stage-I geometry builder only: plate + hole(s), no crack.
%   - It reuses the hole infrastructure from the CrackPath project:
%         holes_to_loops
%         build_pde_geometry_with_holes
%   - Mesh generation currently uses a global Hmax/Hgrad only.
%     The near-hole target size hhole is stored in metadata for future
%     refinement logic, but PDE Toolbox's generateMesh interface here does
%     not enforce a local sizing field directly.

    % -------------------- required input --------------------
    must(C, 'A');
    must(C, 'B');

    A = C.A;
    B = C.B;

    % -------------------- hole input normalization --------------------
    holes = normalize_holes(C);

    if isempty(holes)
        warning('geom_hole_only:NoHoles', ...
            'No holes were provided. Building a plain rectangular plate.');
    end

    % -------------------- mesh options --------------------
    mesh1 = getf(C, 'mesh1', struct());

    Hmin  = getf(mesh1, 'hmin', 0.001);
    Hmax  = getf(mesh1, 'hmax', 0.10);
    Hhole = getf(mesh1, 'hhole', []);
    Hgrad = getf(mesh1, 'hgrad', 1.35);

    % currently only stored for future use
    refineBox = getf(mesh1, 'refineBox', []);

    % -------------------- plotting options --------------------
    plotStruct = getf(C, 'plot', struct());
    showMesh   = getf(plotStruct, 'show_mesh1', false);

    % -------------------- outer boundary --------------------
    % CCW rectangle:
    %   (0,-B) -> (A,-B) -> (A,B) -> (0,B)
    outerPoly = [ ...
        0, -B;
        A, -B;
        A,  B;
        0,  B ];

    % -------------------- hole loops --------------------
    holeLoops = holes_to_loops(holes);

    % -------------------- PDE geometry --------------------
    [mdl, dl, bt, gd, ns, sf] = build_pde_geometry_with_holes(outerPoly, holeLoops);

    % -------------------- mesh generation --------------------
    msh = generateMesh(mdl, ...
        'Hmax', Hmin, ...
        'Hmax', Hmax, ...
        'Hgrad', max(1.01, Hgrad), ...
        'GeometricOrder', 'linear');

    p = msh.Nodes.';      % [np x 2]
    t = msh.Elements.';   % [nt x 3]

    % -------------------- boundary node sets --------------------
    edgeSets = build_edge_sets(p, A, B, holes);

    % -------------------- optional plots --------------------
    if showMesh
        figure; clf; hold on; axis equal; box on

        triplot(t, p(:,1), p(:,2), 'Color', [0.75 0.75 0.75]);

        % outer boundary
        plot([outerPoly(:,1); outerPoly(1,1)], ...
             [outerPoly(:,2); outerPoly(1,2)], ...
             'k-', 'LineWidth', 1.2);

        % hole boundaries
        for ih = 1:numel(holeLoops)
            H = holeLoops{ih};
            plot([H(:,1); H(1,1)], ...
                 [H(:,2); H(1,2)], ...
                 'r-', 'LineWidth', 1.2);
        end

        % mark boundary-node sets if useful
        if ~isempty(edgeSets.hole)
            plot(p(edgeSets.hole,1), p(edgeSets.hole,2), 'ro', ...
                'MarkerSize', 4, 'LineWidth', 1.0);
        end

        title('Stage I geometry: plate with hole(s)');
        xlabel('x');
        ylabel('y');
        xlim([0, A]);
        ylim([-B, B]);
    end

    % -------------------- outputs --------------------
    G = struct();

    G.p         = p;
    G.t         = t;
    G.geom      = mdl;
    G.meshobj   = msh;

    G.outerPoly = outerPoly;
    G.holeLoops = holeLoops;

    if numel(holes) == 1
        G.hole = holes{1};
    else
        G.hole = [];
    end

    G.edgeSets  = edgeSets;

    G.meta = struct();
    G.meta.A         = A;
    G.meta.B         = B;
    G.meta.holes     = holes;
    G.meta.mesh1     = mesh1;
    G.meta.Hmin      = Hmax;
    G.meta.Hmax      = Hmax;
    G.meta.Hhole     = Hhole;
    G.meta.Hgrad     = Hgrad;
    G.meta.refineBox = refineBox;

    G.meta.decsg = struct( ...
        'dl', dl, ...
        'bt', bt, ...
        'gd', gd, ...
        'ns', ns, ...
        'sf', sf);

end


% =========================================================================
% helpers
% =========================================================================

function holes = normalize_holes(C)
%NORMALIZE_HOLES Normalize hole input to a cell array of structs.

    if isfield(C, 'holes') && ~isempty(C.holes)
        holes = C.holes;
    elseif isfield(C, 'hole') && ~isempty(C.hole)
        holes = {C.hole};
    else
        holes = {};
    end

    if isstruct(holes)
        holes = num2cell(holes);
    end

    if ~iscell(holes)
        error('geom_hole_only:BadHoleInput', ...
            'C.holes must be a cell array or struct array; or provide C.hole.');
    end
end


function edgeSets = build_edge_sets(p, A, B, holes)
%BUILD_EDGE_SETS Build node sets on the outer boundary and hole boundary.
%
% Output fields:
%   .left
%   .right
%   .top
%   .bottom
%   .hole
%   .corners.left_bottom
%   .corners.right_bottom
%   .corners.left_top
%   .corners.right_top

    x = p(:,1);
    y = p(:,2);

    Lx = max(A, 1.0);
    Ly = max(2*B, 1.0);
    tol = 1e-8 * max(Lx, Ly);

    % Outer rectangle
    edgeSets = struct();

    edgeSets.left   = find(abs(x - 0) < tol);
    edgeSets.right  = find(abs(x - A) < tol);
    edgeSets.bottom = find(abs(y + B) < tol);
    edgeSets.top    = find(abs(y - B) < tol);

    % Corner anchor candidates
    edgeSets.corners = struct();
    edgeSets.corners.left_bottom  = nearest_node(p, [0, -B]);
    edgeSets.corners.right_bottom = nearest_node(p, [A, -B]);
    edgeSets.corners.left_top     = nearest_node(p, [0,  B]);
    edgeSets.corners.right_top    = nearest_node(p, [A,  B]);

    % Hole boundary nodes
    edgeSets.hole = [];

    if isempty(holes)
        return;
    end

    hole_ids = false(size(p,1),1);

    for ih = 1:numel(holes)
        hk = holes{ih};

        if ~isstruct(hk) || ~isfield(hk, 'type')
            error('geom_hole_only:BadHoleSpec', ...
                'Each hole must be a struct with field "type".');
        end

        switch lower(strtrim(hk.type))
            case 'circle'
                must_field(hk, 'center', 'geom_hole_only:BadHoleSpec');
                must_field(hk, 'r',      'geom_hole_only:BadHoleSpec');

                c = hk.center(:).';
                r = hk.r;

                rr = hypot(x - c(1), y - c(2));
                hole_ids = hole_ids | (abs(rr - r) < 5*tol);

            case 'polygon'
                % Conservative fallback: nodes near any polygon edge
                must_field(hk, 'vertices', 'geom_hole_only:BadHoleSpec');
                V = hk.vertices;
                ids = nodes_near_polygon_edges(p, V, 5*tol);
                hole_ids(ids) = true;

            otherwise
                error('geom_hole_only:UnsupportedHoleType', ...
                    'Unsupported hole type "%s".', hk.type);
        end
    end

    edgeSets.hole = find(hole_ids);
end


function idx = nearest_node(p, xq)
%NEAREST_NODE Index of nearest node to query point xq.

    d2 = sum((p - xq).^2, 2);
    [~, idx] = min(d2);
end


function ids = nodes_near_polygon_edges(p, V, tol)
%NODES_NEAR_POLYGON_EDGES Return node ids near polygon edges.

    if isempty(V)
        ids = [];
        return;
    end

    if norm(V(end,:) - V(1,:), inf) ~= 0
        Vc = [V; V(1,:)];
    else
        Vc = V;
    end

    keep = false(size(p,1),1);

    for k = 1:size(Vc,1)-1
        A = Vc(k,:);
        B = Vc(k+1,:);
        keep = keep | point_near_segment(p, A, B, tol);
    end

    ids = find(keep);
end


function mask = point_near_segment(P, A, B, tol)
%POINT_NEAR_SEGMENT Logical mask of points near segment AB.

    AB = B - A;
    L2 = max(dot(AB,AB), 1e-30);

    t = ((P - A) * AB.') / L2;
    t = max(0, min(1, t));

    Q = A + t .* AB;
    d2 = sum((P - Q).^2, 2);

    mask = d2 <= tol^2;
end


function must(S, field)
%MUST Error if field does not exist or is empty.

    if ~isfield(S, field) || isempty(S.(field))
        error('geom_hole_only:MissingField', ...
            'Required field C.%s is missing or empty.', field);
    end
end


function must_field(S, field, errid)
%MUST_FIELD Error if struct field does not exist or is empty.

    if ~isfield(S, field) || isempty(S.(field))
        error(errid, 'Required field "%s" is missing or empty.', field);
    end
end


function v = getf(S, field, default)
%GETF Get struct field or default.

    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        v = S.(field);
    else
        v = default;
    end
end