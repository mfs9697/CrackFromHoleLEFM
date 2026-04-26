function ids = identify_pencil_ids_from_temp_mesh(mdl, msh, D, varargin)
%IDENTIFY_PENCIL_IDS_FROM_TEMP_MESH
% Recover TipVertexID and the two sharp-pencil EdgeIDs from a temporary mesh.
%
% This is a fallback when pure geometry-edge evaluation is unreliable.

    ip = inputParser;
    addParameter(ip, 'Tol', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(ip, 'Verbose', false, @(x)islogical(x) && isscalar(x));
    parse(ip, varargin{:});

    tolIn   = ip.Results.Tol;
    verbose = ip.Results.Verbose;

    must_field(D, 'channelGeom');
    must_field(D.channelGeom, 'append');

    Gapp = D.channelGeom.append;
    must_field(Gapp, 'xtip');
    must_field(Gapp, 'face_upper');
    must_field(Gapp, 'face_lower');

    xtip = Gapp.xtip(:).';
    segU = Gapp.face_upper;
    segL = Gapp.face_lower;

    mouthUpper = farther_endpoint(segU, xtip);
    mouthLower = farther_endpoint(segL, xtip);

    p = msh.Nodes.';
    geom = mdl.Geometry;

    V = get_geom_vertices_2d(geom);
    if isempty(V)
        error('identify_pencil_ids_from_temp_mesh:NoVertices', ...
            'Could not access geometry vertices from mdl.Geometry.');
    end

    [v_tip, d_tip] = nearest_vertex_to_point(V, xtip);
    [v_up,  d_up ] = nearest_vertex_to_point(V, mouthUpper);
    [v_lo,  d_lo ] = nearest_vertex_to_point(V, mouthLower);

    if isempty(tolIn)
        scale = max([1, norm(xtip), norm(mouthUpper-xtip), norm(mouthLower-xtip)]);
        tol = 1e-8 * scale;
    else
        tol = tolIn;
    end

    n_tip = unique(findNodes(msh, 'region', 'Vertex', v_tip));
    if isempty(n_tip)
        n_tip = nearest_node(p, xtip);
    else
        n_tip = n_tip(1);
    end

    nEdges = get_num_edges(geom);

    candE   = [];
    candFar = [];

    for eid = 1:nEdges
        nEdge = unique(findNodes(msh, 'region', 'Edge', eid));
        if isempty(nEdge)
            continue;
        end

        isIncident = any(nEdge == n_tip) || min(vecnorm(p(nEdge,:) - xtip, 2, 2)) <= 50*tol;
        if ~isIncident
            continue;
        end

        [~, imax] = max(vecnorm(p(nEdge,:) - xtip, 2, 2));
        Pfar = p(nEdge(imax), :);

        candE(end+1,1)   = eid;  %#ok<AGROW>
        candFar(end+1,:) = Pfar; %#ok<AGROW>
    end

    if numel(candE) < 2
        error('identify_pencil_ids_from_temp_mesh:TooFewCandidateEdges', ...
            'Could not find two tip-incident candidate edges in the temporary mesh.');
    end

    costU = vecnorm(candFar - mouthUpper, 2, 2);
    costL = vecnorm(candFar - mouthLower, 2, 2);

    bestCost = inf;
    bestPair = [NaN NaN];

    for i = 1:numel(candE)
        for j = 1:numel(candE)
            if i == j
                continue;
            end
            cij = costU(i) + costL(j);
            if cij < bestCost
                bestCost = cij;
                bestPair = [i j];
            end
        end
    end

    e_upper = candE(bestPair(1));
    e_lower = candE(bestPair(2));

    if verbose
        fprintf('identify_pencil_ids_from_temp_mesh: v_tip = %d, dist = %.3e\n', v_tip, d_tip);
        fprintf('identify_pencil_ids_from_temp_mesh: v_up  = %d, dist = %.3e\n', v_up,  d_up);
        fprintf('identify_pencil_ids_from_temp_mesh: v_lo  = %d, dist = %.3e\n', v_lo,  d_lo);
        fprintf('identify_pencil_ids_from_temp_mesh: e_upper = %d\n', e_upper);
        fprintf('identify_pencil_ids_from_temp_mesh: e_lower = %d\n', e_lower);
    end

    ids = struct();
    ids.v_tip  = v_tip;
    ids.v_up   = v_up;
    ids.v_lo   = v_lo;
    ids.e_upper = e_upper;
    ids.e_lower = e_lower;
    ids.e_tip   = [e_upper e_lower];
    ids.xtip    = xtip;
    ids.mouthUpper = mouthUpper;
    ids.mouthLower = mouthLower;
    ids.vertex_dist_tip = d_tip;
    ids.vertex_dist_up  = d_up;
    ids.vertex_dist_lo  = d_lo;
end


function V = get_geom_vertices_2d(geom)
    V = [];
    if isprop(geom, 'Vertices')
        VV = geom.Vertices;
        if isnumeric(VV)
            if size(VV,1) == 2
                V = VV.';
            elseif size(VV,2) == 2
                V = VV;
            end
        end
    end
end


function nEdges = get_num_edges(geom)
    if isprop(geom, 'NumEdges')
        nEdges = geom.NumEdges;
    else
        error('identify_pencil_ids_from_temp_mesh:NoNumEdges', ...
            'Could not access geom.NumEdges.');
    end
end


function [ivid, dmin] = nearest_vertex_to_point(V, x)
    d = vecnorm(V - x, 2, 2);
    [dmin, ivid] = min(d);
end


function Pfar = farther_endpoint(seg, xref)
    d1 = norm(seg(1,:) - xref);
    d2 = norm(seg(2,:) - xref);
    if d1 >= d2
        Pfar = seg(1,:);
    else
        Pfar = seg(2,:);
    end
end


function idx = nearest_node(p, xq)
    d2 = sum((p - xq).^2, 2);
    [~, idx] = min(d2);
end


function must_field(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('identify_pencil_ids_from_temp_mesh:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end