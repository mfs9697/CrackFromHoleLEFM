function Sol = solve_crack_path_polyline_field(C, sigma0)
%SOLVE_CRACK_PATH_POLYLINE_FIELD
% Solve ONE elastic FEM field for the Crack-Path polyline-crack control.
%
% The purpose is to separate the FEM solution from SIF extraction.  The
% returned Sol is then passed unchanged to both:
%   SIF_LEFM_circle2_debug
%   SIF_LEFM_interaction_EDI
%
% This preserves the essential scientific requirement that old/new SIF
% methods see exactly the same mesh and displacement field.
%
% The rigid-body constraints follow mfs9697/Crack-Path:geom_pencil.m:
%   uy(A,0)=0, ux(0,B)=0, ux(0,-B)=0.
%
% Loading follows the Crack-Path LEFM workflow:
%   uniform remote tension on y=+B and y=-B through edge_loads_T6.

    if nargin < 2 || isempty(sigma0)
        sigma0 = getf(C, 'sigma0', 1.0);
    end

    M = build_crack_path_polyline_LEFM_mesh(C);

    mesh = struct();
    mesh.coord3 = M.coord3;
    mesh.connect3 = M.connect3;
    mesh.coord = M.coord;
    mesh.connect = M.connect;

    E = get_material_E(C);
    nu = C.nu;
    ps = getf(C, 'ps', 1);

    if ps == 1
        coef = E / ((1+nu)*(1-2*nu));
        D = coef * [ ...
            1-nu, nu, 0; ...
            nu, 1-nu, 0; ...
            0, 0, (1-2*nu)/2 ];
    else
        coef = E/(1-nu^2);
        D = coef * [ ...
            1, nu, 0; ...
            nu, 1, 0; ...
            0, 0, (1-nu)/2 ];
    end

    mat = struct();
    mat.E = E;
    mat.nu = nu;
    mat.ps = ps;
    mat.G12 = E/(2*(1+nu));
    mat.D = D;
    mat.Dmat = D;

    [nip2, xip2, w2, Nextr] = integr();
    quad = struct('nip2',nip2,'xip2',xip2,'w2',w2,'Nextr',Nextr);

    coord = mesh.coord;
    A = C.A;
    B = C.B;

    fix_pts = [A,0; 0,B; 0,-B];
    fix = zeros(3,1);

    for i = 1:3
        [~,fix(i)] = min((coord(:,1)-fix_pts(i,1)).^2 + ...
                         (coord(:,2)-fix_pts(i,2)).^2);
    end

    if numel(unique(fix)) < 3
        error('solve_crack_path_polyline_field:BCNodeCollision', ...
            'Rigid-body anchor points mapped to fewer than three nodes.');
    end

    fixvar = [ ...
        2*fix(1); ...       % uy(A,0)=0
        2*fix(2)-1; ...     % ux(0,B)=0
        2*fix(3)-1 ];       % ux(0,-B)=0

    K = stif_assem(mesh, mat, quad, fixvar);

    eps1 = getf(C, 'eps1', 1e-9);
    elod = edge_loads_T6(coord, B, eps1);

    ndof = 2*size(coord,1);
    F = zeros(ndof,1);

    if ~isempty(elod)
        ids = elod(:,1);
        w = elod(:,2);
        F(2*ids) = F(2*ids) + sigma0*w;
    end

    F(fixvar) = 0;

    U = K \ F;

    zeroRows = find(sum(abs(K),2) == 0);
    if ~isempty(zeroRows)
        warning('solve_crack_path_polyline_field:ZeroRows', ...
            '%d zero row(s) remain in the stiffness matrix.', numel(zeroRows));
    end

    Sol = struct();
    Sol.mesh = mesh;
    Sol.U = U;
    Sol.K = K;
    Sol.F = F;

    Sol.V = C.Pmid;
    Sol.mat = mat;
    Sol.quad = quad;

    Sol.sigma0 = sigma0;
    Sol.fixvar = fixvar;
    Sol.fixNodes = fix;
    Sol.elod = elod;

    Sol.geometry = M;
    Sol.C = C;

    Sol.audit = struct();
    Sol.audit.fem_source = 'mfs9697/Crack-Path geometry/load/BC logic';
    Sol.audit.sif_extraction_performed = false;
end


function E = get_material_E(C)
    if isfield(C,'E2') && ~isempty(C.E2)
        E = C.E2;
    elseif isfield(C,'E') && ~isempty(C.E)
        E = C.E;
    else
        error('solve_crack_path_polyline_field:MissingE', ...
            'C must contain E2 or E.');
    end
end


function v = getf(S, field, default)
    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        v = S.(field);
    else
        v = default;
    end
end
