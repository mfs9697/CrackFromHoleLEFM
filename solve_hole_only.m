function S1 = solve_hole_only(C, G, varargin)
%SOLVE_HOLE_ONLY  Solve the elastic plate-with-hole problem.
%
%   S1 = solve_hole_only(C, G)
%   S1 = solve_hole_only(C, G, 'lambda', lambda)
%
% Input:
%   C   configuration struct from cfg_hole_initiation
%   G   geometry struct from geom_hole_only
%
% Name-value pairs:
%   'lambda'    scalar load factor (default 1.0)
%
% Output:
%   S1  struct with fields:
%       .lambda      applied load factor
%       .U           displacement vector [ndof x 1]
%       .K           global stiffness matrix
%       .F           global load vector
%       .R           reaction vector = K*U - F
%       .stress      nodal stresses [nnod x 3] = [sxx syy sxy]
%       .coord_def   deformed coordinates from StressExt
%       .mat         material struct
%       .quad        quadrature struct
%       .mesh        T6 mesh struct with fields .coord, .connect
%       .bc          boundary-condition info
%       .load        load info
%
% Notes:
%   - This function expects the helper files to be on the MATLAB path:
%         T3toT6_fast.m
%         integr.m
%         stif_assem.m
%         StressExt.m
%         edge_loads_T6.m
%   - Current loading option implemented:
%         C.load.type = 'remote_tension_y'
%     which applies +sig0 on the top edge and -sig0 on the bottom edge
%     in the y-direction, scaled by lambda.
%
%   - Rigid-body constraints:
%         C.bc.anchor_mode = 'minimal'
%     means:
%         left-bottom corner : ux = 0, uy = 0
%         right-bottom corner: uy = 0

    p = inputParser;
    addParameter(p, 'lambda', 1.0, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});
    lambda = p.Results.lambda;

    % ------------------------------------------------------------
    % checks
    % ------------------------------------------------------------
    must(G, 'p');
    must(G, 't');
    must(C, 'E');
    must(C, 'nu');

    % ------------------------------------------------------------
    % build T6 mesh
    % ------------------------------------------------------------
    coord = G.p;
    tria  = G.t;

    if size(tria,2) == 3
        [coord6, tria6] = T3toT6_fast(coord, tria);
    elseif size(tria,2) == 6
        coord6 = coord;
        tria6  = tria;
    else
        error('solve_hole_only:BadMesh', ...
            'G.t must have 3 or 6 columns. Got %d.', size(tria,2));
    end

    mesh = struct();
    mesh.coord   = coord6;
    mesh.connect = tria6;

    ndof = 2 * size(mesh.coord,1);

    % ------------------------------------------------------------
    % material matrix
    % ------------------------------------------------------------
    E  = C.E;
    nu = C.nu;
    ps = getf(C, 'ps', 1);   % 1 = plane strain, 0 = plane stress

    mat = struct();
    mat.E  = E;
    mat.nu = nu;
    mat.ps = ps;

    if ps == 1
        % plane strain
        coef = E / ((1 + nu) * (1 - 2*nu));
        D = coef * [ ...
            1 - nu,   nu,          0;
            nu,       1 - nu,      0;
            0,        0,   (1 - 2*nu)/2 ];
    else
        % plane stress
        coef = E / (1 - nu^2);
        D = coef * [ ...
            1,    nu,   0;
            nu,   1,    0;
            0,    0,  (1 - nu)/2 ];
    end
    mat.D = D;

    % ------------------------------------------------------------
    % quadrature
    % ------------------------------------------------------------
    [nip2, xip2, w2, Nextr] = integr();

    quad = struct();
    quad.nip2  = nip2;
    quad.xip2  = xip2;
    quad.w2    = w2;
    quad.Nextr = Nextr;

    % ------------------------------------------------------------
    % boundary conditions (minimal anchoring)
    % ------------------------------------------------------------
    bc = struct();
    bc.anchor_mode = getf(getf(C, 'bc', struct()), 'anchor_mode', 'minimal');

    switch lower(strtrim(bc.anchor_mode))
        case 'minimal'
            iLB = get_corner_node(G, 'left_bottom');
            iRB = get_corner_node(G, 'right_bottom');

            fixvar = [ ...
                2*iLB - 1;   % ux(left-bottom) = 0
                2*iLB;       % uy(left-bottom) = 0
                2*iRB        % uy(right-bottom) = 0
                ];

            bc.corner_nodes.left_bottom  = iLB;
            bc.corner_nodes.right_bottom = iRB;

        otherwise
            error('solve_hole_only:UnknownAnchorMode', ...
                'Unsupported C.bc.anchor_mode = "%s".', bc.anchor_mode);
    end

    bc.fixvar = unique(fixvar(:));

    % ------------------------------------------------------------
    % stiffness matrix
    % ------------------------------------------------------------
    K = stif_assem(mesh, mat, quad, bc.fixvar);

    % ------------------------------------------------------------
    % load vector
    % ------------------------------------------------------------
    loadC = getf(C, 'load', struct());
    loadType = getf(loadC, 'type', 'remote_tension_y');
    sig0     = getf(loadC, 'sig0', 1.0);

    F = zeros(ndof,1);

    switch lower(strtrim(loadType))
        case 'remote_tension_y'
            % The helper returns signed nodal line-load weights for the
            % top (+) and bottom (-) edges for unit traction magnitude.
            %
            % We then multiply by (lambda * sig0).
            A = getf(C, 'A', max(mesh.coord(:,1)));
            B = getf(C, 'B', max(abs(mesh.coord(:,2))));
            eps1 = 1e-8 * max(A, 2*B);

            elod = edge_loads_T6(mesh.coord, B, eps1);

            qy = lambda * sig0;
            for k = 1:size(elod,1)
                node = elod(k,1);
                w    = elod(k,2);
                F(2*node) = F(2*node) + qy * w;
            end

            loadInfo = struct();
            loadInfo.type = loadType;
            loadInfo.sig0 = sig0;
            loadInfo.qy   = qy;
            loadInfo.elod = elod;

        otherwise
            error('solve_hole_only:UnknownLoadType', ...
                'Unsupported C.load.type = "%s".', loadType);
    end

    % enforce homogeneous Dirichlet values in RHS
    F(bc.fixvar) = 0;

    % ------------------------------------------------------------
    % solve
    % ------------------------------------------------------------
    solverC = getf(C, 'solver', struct());
    linear_solver = getf(solverC, 'linear_solver', 'backslash');
    verbose       = getf(solverC, 'verbose', 1);

    switch lower(strtrim(linear_solver))
        case 'backslash'
            U = K \ F;
        otherwise
            error('solve_hole_only:UnknownLinearSolver', ...
                'Unsupported solver "%s".', linear_solver);
    end

    % reaction vector (with clamped-row K)
    R = K*U - F;

    % ------------------------------------------------------------
    % stress recovery
    % ------------------------------------------------------------
    [coord_def, stress] = StressExt(mesh, U, mat, quad, 1.0);

    if verbose
        umax = max(abs(U));
        smax = max(abs(stress), [], 'all');
        fprintf('solve_hole_only: ndof=%d, nnod=%d, nelem=%d | ||U||_inf=%.6e | max|stress|=%.6e\n', ...
            ndof, size(mesh.coord,1), size(mesh.connect,1), umax, smax);
    end

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    S1 = struct();

    S1.lambda   = lambda;
    S1.U        = U;
    S1.K        = K;
    S1.F        = F;
    S1.R        = R;

    S1.stress   = stress;
    S1.coord_def = coord_def;

    S1.mat      = mat;
    S1.quad     = quad;
    S1.mesh     = mesh;
    S1.bc       = bc;
    S1.load     = loadInfo;

end


% =========================================================================
% helpers
% =========================================================================

function idx = get_corner_node(G, name)
%GET_CORNER_NODE Return a corner node index from G.edgeSets.corners
% with fallback to nearest geometric point.

    if isfield(G, 'edgeSets') && isfield(G.edgeSets, 'corners') ...
            && isfield(G.edgeSets.corners, name) ...
            && ~isempty(G.edgeSets.corners.(name))
        idx = G.edgeSets.corners.(name);
        return;
    end

    % fallback
    must(G, 'p');
    must(G, 'meta');
    must(G.meta, 'A');
    must(G.meta, 'B');

    A = G.meta.A;
    B = G.meta.B;

    switch lower(strrep(name, '_', ''))
        case 'leftbottom'
            xq = [0, -B];
        case 'rightbottom'
            xq = [A, -B];
        case 'lefttop'
            xq = [0,  B];
        case 'righttop'
            xq = [A,  B];
        otherwise
            error('solve_hole_only:UnknownCornerName', ...
                'Unknown corner name "%s".', name);
    end

    idx = nearest_node(G.p, xq);
end


function idx = nearest_node(p, xq)
%NEAREST_NODE Index of nearest node to point xq.

    d2 = sum((p - xq).^2, 2);
    [~, idx] = min(d2);
end


function must(S, field)
%MUST Error if field does not exist or is empty.

    if ~isfield(S, field) || isempty(S.(field))
        error('solve_hole_only:MissingField', ...
            'Required field "%s" is missing or empty.', field);
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