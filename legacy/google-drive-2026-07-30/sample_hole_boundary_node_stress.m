function H = sample_hole_boundary_node_stress(C, G, S1)
%SAMPLE_HOLE_BOUNDARY_NODE_STRESS Tangential stress at hole-boundary mesh nodes.
%
%   H = sample_hole_boundary_node_stress(C, G, S1)
%
% This samples only actual T6 mesh nodes on the polygonal hole boundary.
% No eps_shift is used. Stresses are recovered by the same
% integration-point extrapolation used in StressExt, then averaged over
% the material elements adjacent to each boundary node.
%
% Output fields:
%   .node_id       T6 node ids on the hole boundary
%   .x             node coordinates [n x 2]
%   .phi           polar angle about the circular hole center, radians
%   .phi_deg       same angle, degrees
%   .sxx,.syy,.sxy Cartesian stress components at the node
%   .sig_tt        tangential stress at the node
%   .sig_nn        radial normal stress at the node
%   .sig_nt        radial-tangential shear stress at the node
%   .sig_tt_eff    effective tangential stress used in initiation criterion
%   .n_out         radial normals from hole center into material
%   .t_hat         tangents (CCW orientation)
%   .hole          circular hole specification
%   .stress_recovery label for the recovery method

    must(C,  'stage1');
    must(S1, 'mesh');
    must(S1, 'U');
    must(S1, 'mat');
    must(S1, 'quad');

    hole = get_single_circular_hole(C, G);

    coord = S1.mesh.coord;
    connect = S1.mesh.connect;

    if size(connect,2) ~= 6
        error('sample_hole_boundary_node_stress:NotT6', ...
            'S1.mesh.connect must be a T6 connectivity array.');
    end

    nodeIDs = hole_boundary_node_ids_from_loops(coord, G.holeLoops);
    if isempty(nodeIDs)
        error('sample_hole_boundary_node_stress:NoHoleNodes', ...
            'Could not identify T6 nodes on the polygonal hole boundary.');
    end

    elemStress = element_extrapolated_stresses_T6(S1.mesh, S1.U, S1.mat, S1.quad);
    nodeStress = average_extrapolated_stress_at_nodes(connect, elemStress, nodeIDs);

    c = hole.center(:).';
    x = coord(nodeIDs,:);
    phi = mod(atan2(x(:,2) - c(2), x(:,1) - c(1)), 2*pi);

    sig_tt = nan(numel(nodeIDs),1);
    sig_nn = nan(numel(nodeIDs),1);
    sig_nt = nan(numel(nodeIDs),1);

    for k = 1:numel(nodeIDs)
        n = [cos(phi(k)); sin(phi(k))];
        t = [-sin(phi(k)); cos(phi(k))];

        S = [nodeStress(k,1), nodeStress(k,3);
             nodeStress(k,3), nodeStress(k,2)];

        sig_nn(k) = n.' * S * n;
        sig_tt(k) = t.' * S * t;
        sig_nt(k) = n.' * S * t;
    end

    [phi, order] = sort(phi);

    H = struct();
    H.node_id = nodeIDs(order);
    H.x = x(order,:);
    H.phi = phi;
    H.phi_deg = rad2deg(phi);

    nodeStress = nodeStress(order,:);
    H.sxx = nodeStress(:,1);
    H.syy = nodeStress(:,2);
    H.sxy = nodeStress(:,3);

    H.sig_tt = sig_tt(order);
    H.sig_nn = sig_nn(order);
    H.sig_nt = sig_nt(order);
    H.sig_tt_eff = H.sig_tt;
    H.n_out = [cos(phi), sin(phi)];
    H.t_hat = [-sin(phi), cos(phi)];

    H.hole = hole;
    H.eps_shift = 0.0;
    H.stress_recovery = 'element_extrapolated_boundary_node_average';
end


function nodeStress = average_extrapolated_stress_at_nodes(connect, elemStress, nodeIDs)
%AVERAGE_EXTRAPOLATED_STRESS_AT_NODES Average adjacent element contributions.

    nodeIDs = nodeIDs(:);
    nodeStress = nan(numel(nodeIDs), 3);

    for i = 1:numel(nodeIDs)
        nid = nodeIDs(i);
        [elemIdx, localIdx] = find(connect == nid);

        if isempty(elemIdx)
            error('sample_hole_boundary_node_stress:NodeNotInElements', ...
                'Node %d is not referenced by S1.mesh.connect.', nid);
        end

        vals = nan(numel(elemIdx), 3);
        for k = 1:numel(elemIdx)
            vals(k,:) = squeeze(elemStress(elemIdx(k), localIdx(k), :)).';
        end

        nodeStress(i,:) = mean(vals, 1);
    end
end


function ids = hole_boundary_node_ids_from_loops(coord, holeLoops)
%HOLE_BOUNDARY_NODE_IDS_FROM_LOOPS Find nodes on polygonal hole edges.

    ids = [];

    if isempty(holeLoops)
        return;
    end

    scale = max(1, max(abs(coord), [], 'all'));
    tol = 1e-8 * scale;

    keep = false(size(coord,1), 1);

    for ih = 1:numel(holeLoops)
        H = holeLoops{ih};
        if isempty(H)
            continue;
        end

        if norm(H(end,:) - H(1,:), inf) ~= 0
            Hc = [H; H(1,:)];
        else
            Hc = H;
        end

        for k = 1:size(Hc,1)-1
            d = point_segment_distance(coord, Hc(k,:), Hc(k+1,:));
            keep = keep | (d <= tol);
        end
    end

    ids = find(keep);
end


function elemStress = element_extrapolated_stresses_T6(mesh, U, mat, quad)
%ELEMENT_EXTRAPOLATED_STRESSES_T6 Stress extrapolated to each element's nodes.

    coord = mesh.coord;
    connect = mesh.connect;
    Dmat = mat.D;

    nip2 = quad.nip2;
    xip2 = quad.xip2;
    Nextr = quad.Nextr;

    if size(Nextr,1) ~= 6 || size(Nextr,2) ~= nip2
        error('sample_hole_boundary_node_stress:BadNextr', ...
            'Expected quad.Nextr to be 6 x quad.nip2.');
    end
    if ~isequal(size(xip2), [2, nip2])
        error('sample_hole_boundary_node_stress:BadXip2', ...
            'Expected quad.xip2 to be 2 x quad.nip2.');
    end

    nelem = size(connect,1);
    elemStress = zeros(nelem, 6, 3);

    for elem = 1:nelem
        nodes = connect(elem,:);
        X = coord(nodes,:);

        eldof = zeros(12,1);
        for n = 1:6
            eldof(2*n-1:2*n) = [2*nodes(n)-1; 2*nodes(n)];
        end
        u_elem = U(eldof);

        stress_gp = zeros(nip2,3);
        for ig = 1:nip2
            [Bmat, DetJ] = BN_local(xip2(:,ig), X);
            if DetJ <= 0
                error('sample_hole_boundary_node_stress:BadElement', ...
                    'Inverted or degenerate T6 element e=%d.', elem);
            end

            strain = Bmat * u_elem;
            stress_gp(ig,:) = (Dmat * strain).';
        end

        for s = 1:3
            elemStress(elem,:,s) = (Nextr * stress_gp(:,s)).';
        end
    end
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


function [B, Det] = BN_local(xi0, X)
%BN_LOCAL T6 strain-displacement matrix.

    xi  = [xi0; 1 - sum(xi0)];
    Nap = [ 4*xi(1)-1, 0,          1-4*xi(3), 4*xi(2),          -4*xi(2),         4*xi(3)-4*xi(1);
            0,          4*xi(2)-1, 1-4*xi(3), 4*xi(1),           4*xi(3)-4*xi(2), -4*xi(1) ];

    dxdxi = Nap * X;
    Det = det(dxdxi);
    N1 = [ dxdxi(2,2), -dxdxi(1,2);
          -dxdxi(2,1),  dxdxi(1,1) ] / Det * Nap;

    eldf = 12;
    inx = (2:2:eldf)' - 1;
    iny = inx + 1;

    B = zeros(3, eldf);
    B(1, inx) = N1(1,:);
    B(2, iny) = N1(2,:);
    B(3, inx) = N1(2,:);
    B(3, iny) = N1(1,:);
end


function hole = get_single_circular_hole(C, G)
%GET_SINGLE_CIRCULAR_HOLE Return the single circular hole spec.

    if isfield(G, 'hole') && ~isempty(G.hole)
        hole = G.hole;
    elseif isfield(C, 'hole') && ~isempty(C.hole)
        hole = C.hole;
    elseif isfield(C, 'holes') && numel(C.holes) == 1
        hole = C.holes{1};
    else
        error('sample_hole_boundary_node_stress:HoleSpec', ...
            'This first draft supports exactly one circular hole.');
    end

    if ~isstruct(hole) || ~isfield(hole, 'type')
        error('sample_hole_boundary_node_stress:BadHoleSpec', ...
            'Hole specification must be a struct with field "type".');
    end

    if ~strcmpi(strtrim(hole.type), 'circle')
        error('sample_hole_boundary_node_stress:UnsupportedHoleType', ...
            'Only circular holes are supported in this first draft.');
    end

    must(hole, 'center');
    must(hole, 'r');
end


function must(S, field)
%MUST Error if field does not exist or is empty.

    if ~isfield(S, field) || isempty(S.(field))
        error('sample_hole_boundary_node_stress:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end
