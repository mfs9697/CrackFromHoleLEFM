function Mc = collapse_pencil_faces_to_midline(M, D, varargin)
%COLLAPSE_PENCIL_FACES_TO_MIDLINE
% Collapse the two appended-pencil faces onto the crack midline by moving
% boundary-face node coordinates, while keeping node IDs/topology distinct.
%
% Usage:
%   Mc = collapse_pencil_faces_to_midline(M, D, ...
%       'EdgeIDs', [eUpper eLower], ...
%       'TipVertexID', vTip)
%
% Inputs:
%   M   mesh struct from mesh_hole_pencil_domain
%       required fields: .p, .t, .meshobj
%   D   geometry-description struct from build_domain_hole_pencil_polyline
%       required fields: .Pmid, .channelGeom.append
%
% Name-value options:
%   'EdgeIDs'      : [eUpper eLower] PDE geometry edge labels for the two
%                    pencil faces (required for current workflow)
%   'TipVertexID'  : PDE geometry vertex label of the sharp tip (optional)
%   'Tol'          : geometric tolerance (default = 1e-12)
%   'ProjectMode'  : 'segment' (default) or 'polyline'
%
% Output:
%   Mc  collapsed-mesh struct
%       .p              collapsed node coordinates
%       .t              unchanged element connectivity
%       .p0             original coordinates
%       .geom0          original PDE geometry object
%       .meshobj0       original PDE mesh object
%       .edgeSets       original edge sets plus crack face sets
%       .crack          metadata for the collapsed crack faces

    ip = inputParser;
    addParameter(ip, 'EdgeIDs', [], @(x)isnumeric(x) && numel(x)==2);
    addParameter(ip, 'TipVertexID', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'Tol', 1e-12, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'ProjectMode', 'segment', @(s)ischar(s) || isstring(s));
    parse(ip, varargin{:});

    edgeIDs  = ip.Results.EdgeIDs(:).';
    vTipID   = ip.Results.TipVertexID;
    tol      = ip.Results.Tol;
    projMode = lower(char(ip.Results.ProjectMode));

    % -------------------- checks --------------------
    must(M, 'p');
    must(M, 't');
    must(M, 'meshobj');

    must(D, 'Pmid');
    must(D, 'channelGeom');
    must(D.channelGeom, 'append');

    if isempty(edgeIDs)
        error('collapse_pencil_faces_to_midline:MissingEdgeIDs', ...
            'EdgeIDs = [eUpper eLower] must be provided.');
    end

    p0 = M.p;
    t  = M.t;
    msh = M.meshobj;

    % Midline definition:
    % prefer the appended-hole metadata if available
    Gapp = D.channelGeom.append;

    if isfield(Gapp, 'A') && ~isempty(Gapp.A)
        x0 = Gapp.A(:).';
    else
        x0 = D.Pmid(1,:);
    end

    if isfield(Gapp, 'xtip') && ~isempty(Gapp.xtip)
        xtip = Gapp.xtip(:).';
    else
        xtip = D.Pmid(end,:);
    end

    % Optional polyline support for future extension
    if strcmp(projMode, 'polyline')
        Pmid = D.Pmid;
    else
        Pmid = [x0; xtip];
    end

    % -------------------- recover crack-face nodes --------------------
    % PDE Toolbox edge labels from the meshed geometry
    eUpper = edgeIDs(1);
    eLower = edgeIDs(2);

    nUpper = unique(findNodes(msh, 'region', 'Edge', eUpper));
    nLower = unique(findNodes(msh, 'region', 'Edge', eLower));

    if isempty(nUpper) || isempty(nLower)
        error('collapse_pencil_faces_to_midline:EmptyFaceNodeSet', ...
            'Could not recover nodes on one or both crack-face edges.');
    end

    % Tip node: ideally from the PDE geometry vertex label
    if ~isempty(vTipID)
        nTip = unique(findNodes(msh, 'region', 'Vertex', vTipID));
        if isempty(nTip)
            nTip = nearest_node(p0, xtip);
        elseif numel(nTip) > 1
            nTip = nTip(1);
        end
    else
        nTip = nearest_node(p0, xtip);
    end

    % Sort face nodes by parameter along the midline
    sUpper = polyline_parameter(p0(nUpper,:), Pmid, tol);
    sLower = polyline_parameter(p0(nLower,:), Pmid, tol);

    [sUpper, iu] = sort(sUpper, 'ascend');
    [sLower, il] = sort(sLower, 'ascend');

    nUpper = nUpper(iu);
    nLower = nLower(il);

    % -------------------- collapse coordinates --------------------
    p = p0;

    targetUpper = polyline_points_from_parameter(Pmid, sUpper);
    targetLower = polyline_points_from_parameter(Pmid, sLower);

    p(nUpper,:) = targetUpper;
    p(nLower,:) = targetLower;
    p(nTip,:)   = xtip;   % enforce exact tip position

    % -------------------- pairing info (diagnostic only) --------------------
    % Pair by nearest parameter, without changing topology
    pairLowerForUpper = nearest_parameter_match(sUpper, sLower);

    crack = struct();
    crack.x0 = x0;
    crack.xtip = xtip;
    crack.Pmid = Pmid;
    crack.edgeIDs = edgeIDs;
    crack.tipVertexID = vTipID;
    crack.tipNode = nTip;

    crack.upperNodes = nUpper;
    crack.lowerNodes = nLower;
    crack.upperS = sUpper;
    crack.lowerS = sLower;

    crack.upperTarget = targetUpper;
    crack.lowerTarget = targetLower;
    crack.lowerMatchForUpper = pairLowerForUpper;

    crack.nUpper = numel(nUpper);
    crack.nLower = numel(nLower);
    crack.sameCount = (numel(nUpper) == numel(nLower));

    % -------------------- output --------------------
    Mc = struct();

    Mc.p  = p;
    Mc.t  = t;
    Mc.p0 = p0;

    % The original PDE geometry / mesh object correspond to the pre-collapse
    % slit geometry, so keep them as archival objects only.
    if isfield(M, 'geom')
        Mc.geom0 = M.geom;
    else
        Mc.geom0 = [];
    end
    Mc.meshobj0 = msh;

    Mc.dl = getfield_if_exists(M, 'dl', []);
    Mc.bt = getfield_if_exists(M, 'bt', []);
    Mc.gd = getfield_if_exists(M, 'gd', []);
    Mc.ns = getfield_if_exists(M, 'ns', []);
    Mc.sf = getfield_if_exists(M, 'sf', []);

    Mc.edgeSets = getfield_if_exists(M, 'edgeSets', struct());
    Mc.edgeSets.crackUpper = nUpper;
    Mc.edgeSets.crackLower = nLower;
    Mc.edgeSets.crackTip   = nTip;

    Mc.region = getfield_if_exists(M, 'region', struct());
    Mc.region.mode = 'collapsed_polyline_crack';

    Mc.crack = crack;
end


% =========================================================================
% helpers
% =========================================================================

function s = polyline_parameter(X, P, tol)
% Return normalized arc-length parameter s in [0,1] for each point in X
% by projection onto the nearest polyline segment.

    if size(P,1) < 2
        error('collapse_pencil_faces_to_midline:BadPolyline', ...
            'Polyline must contain at least two points.');
    end

    segLen = sqrt(sum(diff(P,1,1).^2, 2));
    Ltot = sum(segLen);
    if Ltot <= tol
        error('collapse_pencil_faces_to_midline:DegeneratePolyline', ...
            'Polyline length is too small.');
    end

    s = zeros(size(X,1),1);

    cumL = [0; cumsum(segLen)];

    for i = 1:size(X,1)
        xi = X(i,:);
        bestD = inf;
        bestS = 0;

        for k = 1:size(P,1)-1
            A = P(k,:);
            B = P(k+1,:);
            AB = B - A;
            L2 = dot(AB, AB);

            if L2 <= tol
                continue;
            end

            tk = dot(xi - A, AB) / L2;
            tk = max(0, min(1, tk));

            Q = A + tk*AB;
            dk = norm(xi - Q);

            if dk < bestD
                bestD = dk;
                bestS = (cumL(k) + tk*segLen(k)) / Ltot;
            end
        end

        s(i) = bestS;
    end
end


function X = polyline_points_from_parameter(P, s)
% Map normalized arc-length parameters s in [0,1] back to points on polyline P.

    segLen = sqrt(sum(diff(P,1,1).^2, 2));
    Ltot = sum(segLen);

    cumL = [0; cumsum(segLen)];
    X = zeros(numel(s), 2);

    for i = 1:numel(s)
        si = max(0, min(1, s(i)));
        Li = si * Ltot;

        % last point special case
        if abs(Li - Ltot) <= 1e-14*max(1,Ltot)
            X(i,:) = P(end,:);
            continue;
        end

        k = find(cumL <= Li, 1, 'last');
        if k >= numel(cumL)
            k = numel(cumL)-1;
        end

        if segLen(k) <= 0
            X(i,:) = P(k,:);
            continue;
        end

        tk = (Li - cumL(k)) / segLen(k);
        X(i,:) = P(k,:) + tk*(P(k+1,:) - P(k,:));
    end
end


function idx = nearest_parameter_match(sA, sB)
% For each value in sA, return the index of the closest value in sB.

    idx = zeros(size(sA));
    for i = 1:numel(sA)
        [~, idx(i)] = min(abs(sB - sA(i)));
    end
end


function idx = nearest_node(p, xq)
    d2 = sum((p - xq).^2, 2);
    [~, idx] = min(d2);
end


function must(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('collapse_pencil_faces_to_midline:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end


function v = getfield_if_exists(S, field, default)
    if isstruct(S) && isfield(S, field)
        v = S.(field);
    else
        v = default;
    end
end