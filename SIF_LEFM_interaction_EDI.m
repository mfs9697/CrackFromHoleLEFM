function [KI, KII, Aux] = SIF_LEFM_interaction_EDI(mesh, U, V, mat, domain, varargin)
%SIF_LEFM_INTERACTION_EDI
% Prototype interaction-integral SIF extractor for 2D isotropic LEFM.
%
% Purpose:
%   Compute signed KI and KII from a general, non-symmetric FEM mesh using
%   an equivalent-domain interaction integral with analytical auxiliary
%   mode-I and mode-II crack-tip fields.
%
% Usage:
%   [KI, KII, Aux] = SIF_LEFM_interaction_EDI(mesh, U, V, mat, domain)
%
% Inputs:
%   mesh    struct with fields
%             .coord    [nT6 x 2]
%             .connect  [nElem x 6]
%   U       displacement vector [ux1; uy1; ux2; uy2; ...]
%   V       crack polyline, last point is crack tip
%   mat     struct with .E, .nu, and .Dmat or .D
%   domain  struct with fields:
%             .r_inner
%             .r_outer
%
% Name-value options:
%   'AuxK'          auxiliary SIF amplitude, default 1.0
%   'UsePlaneStrain' true/false. If omitted, uses mat.ps if present;
%                    otherwise assumes plane strain.
%   'Verbose'       true/false
%
% Outputs:
%   KI, KII signed stress intensity factors
%   Aux     diagnostic struct
%
% Notes:
%   This is a first validation implementation. The sign convention must be
%   checked using a centered-hole parity benchmark:
%       KII(0) ~ 0 and KII(-theta) ~ -KII(+theta).

    %% ------------------------------------------------------------
    % 0. Options
    %% ------------------------------------------------------------
    ip = inputParser;
    addParameter(ip, 'AuxK', 1.0, @(x)isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'UsePlaneStrain', [], @(x)islogical(x) || isnumeric(x) || isempty(x));
    addParameter(ip, 'Verbose', false, @(x)islogical(x) || isnumeric(x));
    parse(ip, varargin{:});

    Kaux = ip.Results.AuxK;
    verbose = logical(ip.Results.Verbose);

    %% ------------------------------------------------------------
    % 1. Checks and unpacking
    %% ------------------------------------------------------------
    must_have(mesh, 'coord');
    must_have(mesh, 'connect');

    coord = mesh.coord;
    connect = mesh.connect;

    if size(connect,2) ~= 6
        error('SIF_LEFM_interaction_EDI:BadConnect', ...
            'mesh.connect must be T6 connectivity with 6 columns.');
    end

    if numel(U) ~= 2*size(coord,1)
        error('SIF_LEFM_interaction_EDI:BadU', ...
            'numel(U) must equal 2*size(mesh.coord,1).');
    end

    must_have(mat, 'E');
    must_have(mat, 'nu');

    E = mat.E;
    nu = mat.nu;

    if isfield(mat, 'Dmat') && ~isempty(mat.Dmat)
        Dmat = mat.Dmat;
    elseif isfield(mat, 'D') && ~isempty(mat.D)
        Dmat = mat.D;
    else
        error('SIF_LEFM_interaction_EDI:MissingDmat', ...
            'mat must contain .Dmat or .D.');
    end

    if isempty(ip.Results.UsePlaneStrain)
        if isfield(mat, 'ps') && ~isempty(mat.ps)
            planeStrain = logical(mat.ps == 1);
        else
            planeStrain = true;
        end
    else
        planeStrain = logical(ip.Results.UsePlaneStrain);
    end

    if planeStrain
        Eeff = E/(1 - nu^2);
        kappa = 3 - 4*nu;
    else
        Eeff = E;
        kappa = (3 - nu)/(1 + nu);
    end

    mu = E/(2*(1 + nu));

    must_have(domain, 'r_inner');
    must_have(domain, 'r_outer');

    r_inner = domain.r_inner;
    r_outer = domain.r_outer;

    if ~(r_inner >= 0 && r_outer > r_inner)
        error('SIF_LEFM_interaction_EDI:BadDomain', ...
            'domain must satisfy 0 <= r_inner < r_outer.');
    end

    %% ------------------------------------------------------------
    % 2. Crack-tip frame
    %% ------------------------------------------------------------
    if size(V,1) < 2 || size(V,2) ~= 2
        error('SIF_LEFM_interaction_EDI:BadCrackPolyline', ...
            'V must be [nPts x 2], nPts >= 2.');
    end

    x_tip = V(end,:).';
    p_prev = V(end-1,:).';

    e1 = x_tip - p_prev;
    Le = norm(e1);

    if Le <= eps
        error('SIF_LEFM_interaction_EDI:DegenerateLastSegment', ...
            'Last crack segment has nearly zero length.');
    end

    e1 = e1 / Le;
    e2 = [-e1(2); e1(1)];

    R_gl = [e1, e2];
    R_loc = R_gl.';

    %% ------------------------------------------------------------
    % 3. Quadrature
    %% ------------------------------------------------------------
    [nip2, xip2, w2] = local_integr_T6();

    %% ------------------------------------------------------------
    % 4. Loop over elements and Gauss points
    %% ------------------------------------------------------------
    I_modeI = 0.0;
    I_modeII = 0.0;

    nGP_total = 0;
    nGP_used = 0;
    nElem_used = 0;

    rows = [];

    for e = 1:size(connect,1)

        nodes = connect(e,:);
        X = coord(nodes,:);

        % Quick centroid test for element participation.
        Xc = mean(X(1:3,:), 1).';
        xc_loc = R_loc * (Xc - x_tip);
        rc = norm(xc_loc);

        if rc > r_outer + local_element_radius(X)
            continue;
        end

        elemUsed = false;

        uel = [U(2*nodes - 1), U(2*nodes)];

        for igp = 1:nip2

            nGP_total = nGP_total + 1;

            xi = xip2(:,igp);

            [N, DetJ, dNdx] = local_T6_shape_grad(xi, X);

            if DetJ <= 0
                error('SIF_LEFM_interaction_EDI:BadElement', ...
                    'Inverted or degenerate T6 element e=%d.', e);
            end

            xg = (N(:).' * X).';
            xl = R_loc * (xg - x_tip);

            x1 = xl(1);
            x2 = xl(2);
            r = hypot(x1, x2);

            if r <= r_inner || r >= r_outer
                continue;
            end

            % Exclude points extremely close to the crack tip.
            if r < 1e-12
                continue;
            end

            nGP_used = nGP_used + 1;
            elemUsed = true;

            % Actual FEM gradient in global coordinates.
            dux_dx = dNdx(1,:) * uel(:,1);
            dux_dy = dNdx(2,:) * uel(:,1);
            duy_dx = dNdx(1,:) * uel(:,2);
            duy_dy = dNdx(2,:) * uel(:,2);

            GradU_gl = [dux_dx, dux_dy;
                        duy_dx, duy_dy];

            % Rotate actual gradient to local crack-tip frame.
            GradU1 = R_loc * GradU_gl * R_gl;

            eps1 = [ ...
                GradU1(1,1);
                GradU1(2,2);
                GradU1(1,2) + GradU1(2,1)];

            sig1 = Dmat * eps1;

            du1_dx1 = GradU1(:,1);

            % Weight function q and gradient q_,j in local coordinates.
            % q = 1 at r_inner, q = 0 at r_outer.
            qgrad = local_qgrad_radial(xl, r_inner, r_outer);

            % Auxiliary mode I, normalized by Kaux.
            auxI = local_aux_LEFM_fields(x1, x2, Kaux, 0.0, E, nu, mu, kappa, Dmat);

            % Auxiliary mode II, normalized by Kaux.
            auxII = local_aux_LEFM_fields(x1, x2, 0.0, Kaux, E, nu, mu, kappa, Dmat);

            % Interaction integral densities.
            densI  = local_interaction_density(sig1, eps1, du1_dx1, auxI,  qgrad);
            densII = local_interaction_density(sig1, eps1, du1_dx1, auxII, qgrad);

            dA = w2(igp) * (DetJ/2);

            I_modeI  = I_modeI  + densI  * dA;
            I_modeII = I_modeII + densII * dA;

            rows = [rows; e, igp, x1, x2, r, densI, densII, dA]; %#ok<AGROW>
        end

        if elemUsed
            nElem_used = nElem_used + 1;
        end
    end

    %% ------------------------------------------------------------
    % 5. Convert interaction integrals to SIFs
    %% ------------------------------------------------------------
    % With Kaux = 1:
    %   I_modeI  = KI_actual  * Kaux / Eeff
    %   I_modeII = KII_actual * Kaux / Eeff
    %
    % The sign may need one global convention correction after validation.
    KI  = Eeff * I_modeI  / Kaux;
    KII = Eeff * I_modeII / Kaux;

    %% ------------------------------------------------------------
    % 6. Diagnostics
    %% ------------------------------------------------------------
    Aux = struct();

    Aux.method = 'interaction_EDI_prototype';
    Aux.KI = KI;
    Aux.KII = KII;

    Aux.I_modeI = I_modeI;
    Aux.I_modeII = I_modeII;

    Aux.Kaux = Kaux;
    Aux.Eeff = Eeff;
    Aux.mu = mu;
    Aux.kappa = kappa;
    Aux.planeStrain = planeStrain;

    Aux.r_inner = r_inner;
    Aux.r_outer = r_outer;

    Aux.x_tip = x_tip.';
    Aux.e1 = e1.';
    Aux.e2 = e2.';
    Aux.R_gl = R_gl;
    Aux.R_loc = R_loc;

    Aux.nGP_total = nGP_total;
    Aux.nGP_used = nGP_used;
    Aux.nElem_used = nElem_used;

    if ~isempty(rows)
        Aux.gpTable = array2table(rows, ...
            'VariableNames', {'elem','igp','x1','x2','r','densI','densII','dA'});
    else
        Aux.gpTable = table();
    end

    if verbose
        fprintf('\nSIF_LEFM_interaction_EDI summary:\n');
        fprintf('  r_inner, r_outer = %.6e, %.6e\n', r_inner, r_outer);
        fprintf('  GP used / total   = %d / %d\n', nGP_used, nGP_total);
        fprintf('  elements used     = %d\n', nElem_used);
        fprintf('  I_modeI           = %.8e\n', I_modeI);
        fprintf('  I_modeII          = %.8e\n', I_modeII);
        fprintf('  KI                = %.8e\n', KI);
        fprintf('  KII               = %.8e\n', KII);
        fprintf('  KII/KI            = %.8e\n', KII/KI);
    end
end


% =========================================================================
% Interaction density
% =========================================================================

function dens = local_interaction_density(sig1, eps1, du1_dx1, aux, qgrad)
% Domain form used here:
%
% I = int_A [ -sigma2:eps1 * delta_1j
%             + sigma1_ij u2_i,1
%             + sigma2_ij u1_i,1 ] q_,j dA
%
% All quantities are in local crack coordinates.

    sig2 = aux.sig;
    eps2 = aux.eps;
    du2_dx1 = aux.du_dx1;

    % sigma2:eps1 with engineering shear convention:
    % sig:eps = sig11*eps11 + sig22*eps22 + sig12*gamma12
    interaction_energy = sig2(1)*eps1(1) + sig2(2)*eps1(2) + sig2(3)*eps1(3);

    % Build stress matrices.
    S1 = [sig1(1), sig1(3);
          sig1(3), sig1(2)];

    S2 = [sig2(1), sig2(3);
          sig2(3), sig2(2)];

    % Vector A_j = -Wint*delta_1j + sigma1_ij*u2_i,1 + sigma2_ij*u1_i,1
    A = zeros(2,1);

    A(1) = -interaction_energy ...
           + dot(S1(:,1), du2_dx1) ...
           + dot(S2(:,1), du1_dx1);

    A(2) = dot(S1(:,2), du2_dx1) ...
           + dot(S2(:,2), du1_dx1);

    dens = A.' * qgrad(:);
end


% =========================================================================
% Radial q-gradient
% =========================================================================

function qgrad = local_qgrad_radial(xl, r_inner, r_outer)

    x1 = xl(1);
    x2 = xl(2);
    r = hypot(x1, x2);

    if r <= r_inner || r >= r_outer || r <= eps
        qgrad = [0; 0];
        return;
    end

    % q = (r_outer - r)/(r_outer - r_inner)
    dqdr = -1/(r_outer - r_inner);

    qgrad = dqdr * [x1; x2] / r;
end


% =========================================================================
% Auxiliary LEFM fields by finite-difference derivatives of displacements
% =========================================================================

function aux = local_aux_LEFM_fields(x1, x2, KI, KII, E, nu, mu, kappa, Dmat)
% Return auxiliary stress, strain, and du/dx1 in local crack coordinates.
%
% This prototype evaluates stresses from standard near-tip formulas and
% computes strain/du_dx1 from numerical derivatives of the analytical
% displacement field. This is robust enough for first validation but can be
% replaced later by closed-form derivatives.

    r = hypot(x1, x2);
    th = atan2(x2, x1);

    if r <= 1e-14
        aux.sig = [0;0;0];
        aux.eps = [0;0;0];
        aux.du_dx1 = [0;0];
        return;
    end

    sig = local_aux_stress_polar_to_cart(r, th, KI, KII);

    % Numerical derivative of auxiliary displacement wrt local x1.
    h = max(1e-7, 1e-5*r);

    u0 = local_aux_displacement(x1, x2, KI, KII, mu, kappa); %#ok<NASGU>
    up = local_aux_displacement(x1 + h, x2, KI, KII, mu, kappa);
    um = local_aux_displacement(x1 - h, x2, KI, KII, mu, kappa);

    du_dx1 = (up - um)/(2*h);

    % To obtain eps, compute full displacement gradient numerically.
    vp = local_aux_displacement(x1, x2 + h, KI, KII, mu, kappa);
    vm = local_aux_displacement(x1, x2 - h, KI, KII, mu, kappa);

    du_dx2 = (vp - vm)/(2*h);

    Grad = [du_dx1, du_dx2];

    eps = [ ...
        Grad(1,1);
        Grad(2,2);
        Grad(1,2) + Grad(2,1)];

    % For consistency, one could use eps = inv(Dmat)*sig.
    % But using displacement derivatives makes du_dx1 and eps compatible.
    % If needed, compare with:
    % eps_from_sig = Dmat \ sig;

    aux = struct();
    aux.sig = sig;
    aux.eps = eps;
    aux.du_dx1 = du_dx1;
end


% =========================================================================
% Auxiliary stresses
% =========================================================================

function sig = local_aux_stress_polar_to_cart(r, th, KI, KII)
% Standard isotropic LEFM stresses near a crack tip in local coordinates.

    c = cos(th/2);
    s = sin(th/2);

    ct = cos(th);
    st = sin(th);

    fac = 1/sqrt(2*pi*r);

    % Mode I polar stresses
    srr_I = KI*fac*c*(1 - s*sin(3*th/2));
    stt_I = KI*fac*c*(1 + s*sin(3*th/2));
    srt_I = KI*fac*c*s*cos(3*th/2);

    % Mode II polar stresses
    srr_II = -KII*fac*s*(2 + c*cos(3*th/2));
    stt_II =  KII*fac*s*c*cos(3*th/2);
    srt_II =  KII*fac*c*(1 - s*sin(3*th/2));

    srr = srr_I + srr_II;
    stt = stt_I + stt_II;
    srt = srt_I + srt_II;

    % Polar to Cartesian in local x1,x2.
    s11 = srr*ct^2 + stt*st^2 - 2*srt*st*ct;
    s22 = srr*st^2 + stt*ct^2 + 2*srt*st*ct;
    s12 = (srr - stt)*st*ct + srt*(ct^2 - st^2);

    sig = [s11; s22; s12];
end


% =========================================================================
% Auxiliary displacements
% =========================================================================

function u = local_aux_displacement(x1, x2, KI, KII, mu, kappa)
% Standard isotropic LEFM displacement field in local crack coordinates.
%
% u1,u2 are expressed in local Cartesian coordinates.
%
% The formulas are conventional near-tip expressions. Sign convention for
% mode II should be validated by the centered-hole parity test.

    r = hypot(x1, x2);
    th = atan2(x2, x1);

    if r <= 1e-14
        u = [0;0];
        return;
    end

    fac = sqrt(r/(2*pi))/(2*mu);

    c = cos(th/2);
    s = sin(th/2);

    % Mode I
    u1_I = KI * fac * c * (kappa - 1 + 2*s^2);
    u2_I = KI * fac * s * (kappa + 1 - 2*c^2);

    % Mode II
    u1_II = KII * fac * s * (kappa + 1 + 2*c^2);
    u2_II = -KII * fac * c * (kappa - 1 - 2*s^2);

    u = [u1_I + u1_II;
         u2_I + u2_II];
end


% =========================================================================
% T6 shape functions and gradients
% =========================================================================

function [N, DetJ, dNdx] = local_T6_shape_grad(xi, X)
% xi = [L1; L2], L3 = 1-L1-L2.
% T6 node order assumed:
%   1=L1 vertex, 2=L2 vertex, 3=L3 vertex,
%   4=edge 1-2, 5=edge 2-3, 6=edge 3-1.
%
% This matches the common T3toT6 convention used by the repo.

    L1 = xi(1);
    L2 = xi(2);
    L3 = 1 - L1 - L2;

    N = [ ...
        L1*(2*L1 - 1);
        L2*(2*L2 - 1);
        L3*(2*L3 - 1);
        4*L1*L2;
        4*L2*L3;
        4*L3*L1];

    dN_dL1 = [ ...
        4*L1 - 1;
        0;
        -(4*L3 - 1);
        4*L2;
        -4*L2;
        4*(L3 - L1)];

    dN_dL2 = [ ...
        0;
        4*L2 - 1;
        -(4*L3 - 1);
        4*L1;
        4*(L3 - L2);
        -4*L1];

    dNdxi = [dN_dL1.'; dN_dL2.'];

    J = dNdxi * X;

    DetJ = det(J);

    dNdx = J \ dNdxi;
end


% =========================================================================
% Quadrature for T6 triangle
% =========================================================================

function [nip, xip, w] = local_integr_T6()
% Simple 7-point Dunavant rule on reference triangle.
% Weights sum to 1 for the parent triangle area convention used with DetJ/2.

    nip = 7;

    xip = zeros(2,nip);
    w = zeros(1,nip);

    % centroid
    xip(:,1) = [1/3; 1/3];
    w(1) = 0.225;

    a1 = 0.059715871789770;
    b1 = 0.470142064105115;
    w1 = 0.132394152788506;

    xip(:,2) = [a1; b1];
    xip(:,3) = [b1; a1];
    xip(:,4) = [b1; b1];
    w(2:4) = w1;

    a2 = 0.797426985353087;
    b2 = 0.101286507323456;
    w2 = 0.125939180544827;

    xip(:,5) = [a2; b2];
    xip(:,6) = [b2; a2];
    xip(:,7) = [b2; b2];
    w(5:7) = w2;
end


% =========================================================================
% Utility
% =========================================================================

function rad = local_element_radius(X)

    xc = mean(X(1:3,:), 1);
    rad = max(sqrt(sum((X - xc).^2, 2)));
end


function must_have(S, field)

    if ~isstruct(S) || ~isfield(S, field) || isempty(S.(field))
        error('SIF_LEFM_interaction_EDI:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end