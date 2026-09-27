function M = build_crack_path_polyline_LEFM_mesh(C)
%BUILD_CRACK_PATH_POLYLINE_LEFM_MESH
% Build a T3/T6 mesh for a traction-free polyline crack using the pencil
% channel construction from Crack-Path, then collapse both channel faces to
% the crack midline while preserving distinct node IDs/topology.
%
% This is a verification-only extraction/adaptation of the geometry logic in
% mfs9697/Crack-Path:geom_pencil.m.  Unlike the CZM routine, EVERY segment of
% C.Pmid is treated as a physical traction-free crack segment.
%
% Required C fields:
%   A, B, Pmid, chw, ncoh
%
% Optional C fields:
%   hgrad, hmax_ratio, join, miter_limit, corner_tol, tip,
%   plotGeom, plotMesh
%
% Output M:
%   .coord3, .connect3  collapsed parent T3 mesh
%   .coord,  .connect   T6 mesh
%   .coord3_pre         pre-collapse T3 nodes
%   .Bound, .Ggeom      pencil geometry diagnostics
%   .crack              crack-face diagnostic data
%   .meshobj            original PDE Toolbox mesh object

    must(C, 'A');
    must(C, 'B');
    must(C, 'Pmid');
    must(C, 'chw');
    must(C, 'ncoh');

    A = C.A;
    B = C.B;
    P0 = C.Pmid;
    w = C.chw;
    ncoh = max(1, round(C.ncoh));

    if size(P0,1) < 2 || size(P0,2) ~= 2
        error('build_crack_path_polyline_LEFM_mesh:BadPolyline', ...
            'C.Pmid must be [nPts x 2] with nPts >= 2.');
    end

    join        = getf(C, 'join', 'miter');
    miter_limit = getf(C, 'miter_limit', 6);
    corner_tol  = getf(C, 'corner_tol', 1e-10);
    tip_mode    = getf(C, 'tip', 'point');
    hgrad       = getf(C, 'hgrad', 1.15);
    hmax_ratio  = getf(C, 'hmax_ratio', 40);
    plotGeom    = logical(getf(C, 'plotGeom', false));
    plotMesh    = logical(getf(C, 'plotMesh', false));

    % Subdivision of the last leg is a meshing cue only; it does not change
    % the physical crack path.
    Psub = subdivide_last_leg(P0, ncoh);

    Llast = norm(P0(end,:) - P0(end-1,:));
    htip  = Llast / ncoh;

    [Bound, Ggeom] = build_domain_pencil_polyline( ...
        Psub, A, B, w, ...
        'join', join, ...
        'miter_limit', miter_limit, ...
        'corner_tol', corner_tol, ...
        'tip', tip_mode);

    if plotGeom
        figure(70); clf; hold on; axis equal; box on
        plot([Bound(:,1); Bound(1,1)], [Bound(:,2); Bound(1,2)], ...
            'k-', 'LineWidth', 1.2);
        plot(Psub(:,1), Psub(:,2), 'r.-', 'LineWidth', 1.3);
        plot(Ggeom.up_chain(:,1), Ggeom.up_chain(:,2), 'b--');
        plot(Ggeom.dn_chain(:,1), Ggeom.dn_chain(:,2), 'b--');
        title('SIF audit: pencil geometry before collapse');
        xlim([0 A]); ylim([-B B]);
    end

    % PDE geometry: non-perforated rectangular domain with the pencil slot.
    Gcol = [2; size(Bound,1); Bound(:,1); Bound(:,2)];
    gd = decsg(Gcol);

    mdl = createpde();
    geometryFromEdges(mdl, gd);

    Hmin  = max(0.98*htip, 1e-12);
    Hcap  = min(B/2, B/hmax_ratio);
    Hmax  = max(Hcap, 1.05*htip);
    Hgrad = max(1.01, hgrad);

    msh = generateMesh(mdl, ...
        'Hmin', Hmin, ...
        'Hmax', Hmax, ...
        'Hgrad', Hgrad, ...
        'GeometricOrder', 'linear');

    coord0   = msh.Nodes.';
    connect0 = msh.Elements.';

    % Recover nodes on the two pencil faces geometrically.  Both Ggeom
    % chains run mouth -> tip.
    tol_line = max(1e-3*Hmin, 1e-6*max(w,eps));

    idsU = face_ids_on_polyline(coord0, Ggeom.up_chain, tol_line);
    idsD = face_ids_on_polyline(coord0, Ggeom.dn_chain, tol_line);

    if isempty(idsU) || isempty(idsD)
        error('build_crack_path_polyline_LEFM_mesh:EmptyCrackFace', ...
            'Failed to identify one or both pencil faces.');
    end

    coord3 = coord0;

    [targetU, sU] = project_to_polyline(coord0(idsU,:), Psub);
    [targetD, sD] = project_to_polyline(coord0(idsD,:), Psub);

    coord3(idsU,:) = targetU;
    coord3(idsD,:) = targetD;

    % Enforce exact mouth/tip coordinates where the parameter indicates an
    % endpoint.  Opposing face nodes remain distinct IDs.
    endTol = 1e-10;
    coord3(idsU(sU <= endTol),:) = repmat(Psub(1,:), sum(sU <= endTol), 1);
    coord3(idsD(sD <= endTol),:) = repmat(Psub(1,:), sum(sD <= endTol), 1);
    coord3(idsU(sU >= 1-endTol),:) = repmat(Psub(end,:), sum(sU >= 1-endTol), 1);
    coord3(idsD(sD >= 1-endTol),:) = repmat(Psub(end,:), sum(sD >= 1-endTol), 1);

    % Remove zero-area triangles created by the collapse and orient the
    % surviving T3 mesh counter-clockwise.
    Atri = tri_areas_signed(connect0, coord3);
    areaScale = max(htip, 1e-12)^2;
    keep = abs(Atri) > 1e-12*areaScale;

    connect3 = connect0(keep,:);
    Atri = Atri(keep);

    cw = Atri < 0;
    if any(cw)
        tmp = connect3(cw,2);
        connect3(cw,2) = connect3(cw,3);
        connect3(cw,3) = tmp;
    end

    % Compact T3 node numbering so that no orphan DOFs enter the T6 solve.
    used = unique(connect3(:));
    map = zeros(size(coord3,1),1);
    map(used) = 1:numel(used);

    coord3_compact = coord3(used,:);
    connect3 = map(connect3);

    idsU_compact = map(idsU);
    idsD_compact = map(idsD);
    idsU_compact = unique(idsU_compact(idsU_compact > 0), 'stable');
    idsD_compact = unique(idsD_compact(idsD_compact > 0), 'stable');

    coord3 = coord3_compact;

    [coord6, connect6] = T3toT6_fast(coord3, connect3);

    if plotMesh
        figure(71); clf; hold on; axis equal; box on
        triplot(connect3, coord3(:,1), coord3(:,2));
        plot(P0(:,1), P0(:,2), 'r-', 'LineWidth', 2);
        title('SIF audit: collapsed traction-free crack mesh');
        xlim([0 A]); ylim([-B B]);
    end

    M = struct();
    M.coord3 = coord3;
    M.connect3 = connect3;
    M.coord = coord6;
    M.connect = connect6;

    M.coord3_pre = coord0;
    M.connect3_pre = connect0;

    M.Bound = Bound;
    M.Ggeom = Ggeom;
    M.model = mdl;
    M.meshobj = msh;

    M.Hmin = Hmin;
    M.Hmax = Hmax;
    M.Hgrad = Hgrad;
    M.htip = htip;

    M.crack = struct();
    M.crack.Pmid0 = P0;
    M.crack.Pmid = Psub;
    M.crack.upperNodesT3 = idsU_compact;
    M.crack.lowerNodesT3 = idsD_compact;
    M.crack.upperS = sU;
    M.crack.lowerS = sD;
    M.crack.tip = P0(end,:);

    M.audit = struct();
    M.audit.source = 'adapted from mfs9697/Crack-Path:geom_pencil.m';
    M.audit.all_segments_are_physical_crack = true;
end


% =========================================================================
function ids = face_ids_on_polyline(P, Chain, tol)
    ids = [];

    for k = 1:size(Chain,1)-1
        ik = face_ids_on_segment(P, Chain(k,:), Chain(k+1,:), tol);
        ids = [ids; ik(:)]; %#ok<AGROW>
    end

    ids = unique(ids, 'stable');

    if isempty(ids)
        return;
    end

    [~, s] = project_to_polyline(P(ids,:), Chain);
    [~, I] = sort(s, 'ascend');
    ids = ids(I);
end


function ids = face_ids_on_segment(P, A, B, tol)
    AB = B - A;
    L2 = dot(AB, AB);

    if L2 <= eps
        ids = [];
        return;
    end

    n = [-AB(2), AB(1)];
    nn = norm(n);

    d = abs((P - A) * (n.'/nn));
    t = ((P - A) * AB.') / L2;

    mask = (d <= tol) & (t >= -1e-10) & (t <= 1+1e-10);
    ids = find(mask);
end


function [Q, s] = project_to_polyline(X, P)
% Project each row of X to the nearest segment of P.
% s is normalized arc length along P in [0,1].

    seg = diff(P,1,1);
    segLen = sqrt(sum(seg.^2,2));
    Ltot = sum(segLen);

    if Ltot <= eps
        error('build_crack_path_polyline_LEFM_mesh:DegeneratePolyline', ...
            'Crack polyline has zero total length.');
    end

    cumL = [0; cumsum(segLen)];

    Q = zeros(size(X));
    s = zeros(size(X,1),1);

    for i = 1:size(X,1)
        bestD2 = inf;
        bestQ = P(1,:);
        bestS = 0;

        for k = 1:size(seg,1)
            AB = seg(k,:);
            L2 = dot(AB,AB);
            if L2 <= eps
                continue;
            end

            t = dot(X(i,:) - P(k,:), AB) / L2;
            t = max(0, min(1, t));

            q = P(k,:) + t*AB;
            d2 = sum((X(i,:) - q).^2);

            if d2 < bestD2
                bestD2 = d2;
                bestQ = q;
                bestS = (cumL(k) + t*segLen(k))/Ltot;
            end
        end

        Q(i,:) = bestQ;
        s(i) = bestS;
    end
end


function A = tri_areas_signed(T, X)
    v1 = X(T(:,1),:);
    v2 = X(T(:,2),:);
    v3 = X(T(:,3),:);

    A = 0.5*((v2(:,1)-v1(:,1)).*(v3(:,2)-v1(:,2)) - ...
             (v2(:,2)-v1(:,2)).*(v3(:,1)-v1(:,1)));
end


function must(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('build_crack_path_polyline_LEFM_mesh:MissingField', ...
            'Required field C.%s is missing or empty.', field);
    end
end


function v = getf(S, field, default)
    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        v = S.(field);
    else
        v = default;
    end
end
