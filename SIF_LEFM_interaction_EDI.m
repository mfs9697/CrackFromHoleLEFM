function [KI, KII, Aux] = SIF_LEFM_interaction_EDI(mesh, U, V, mat, domain, varargin)
%SIF_LEFM_INTERACTION_EDI  Equivalent-domain interaction integral for 2D LEFM.
%
% This extractor does not pair mirrored FEM sampling points and therefore
% does not require a mirror-symmetric mesh around the crack tip.
%
% Inputs
%   mesh    struct with .coord [n x 2], .connect [ne x 6]
%   U       displacement vector [ux1;uy1;ux2;uy2;...]
%   V       crack polyline, last point = crack tip
%   mat     .E, .nu, and .Dmat or .D; optional .ps (1 plane strain)
%   domain  .r_inner, .r_outer
%
% Name-value options
%   'AuxK'              auxiliary SIF amplitude (default 1)
%   'UsePlaneStrain'    [] -> infer mat.ps, otherwise true by default
%   'AuxDerivativeMode' 'analytic' (default) or 'finite_difference'
%   'AuxFDRelStep'      relative FD step if requested (default 1e-5)
%   'AuxFDAbsStep'      absolute FD floor (default 1e-10)
%   'Verbose'           print diagnostics
%
% Normalization
%   For homogeneous isotropic elasticity,
%
%     I = (2/Eeff) (KI*Kaux_I + KII*Kaux_II),
%
%   where Eeff = E/(1-nu^2) in plane strain and E in plane stress.
%   Therefore a unit pure-mode auxiliary field gives
%
%     K = (Eeff/2) I / Kaux.
%
% Sign convention
%   Local x1 is the direction of the final crack segment and
%   x2 = [-e1_y,e1_x].  The sign of KII must still be validated by an
%   independent mirror/parity benchmark before production use.

    ip = inputParser;
    addParameter(ip, 'AuxK', 1.0, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'UsePlaneStrain', [], @(x)islogical(x) || isnumeric(x) || isempty(x));
    addParameter(ip, 'AuxDerivativeMode', 'analytic', @(x)ischar(x) || isstring(x));
    addParameter(ip, 'AuxFDRelStep', 1e-5, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'AuxFDAbsStep', 1e-10, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'Verbose', false, @(x)islogical(x) || isnumeric(x));
    parse(ip, varargin{:});

    Kaux = ip.Results.AuxK;
    verbose = logical(ip.Results.Verbose);

    must_have(mesh,'coord');
    must_have(mesh,'connect');
    must_have(mat,'E');
    must_have(mat,'nu');
    must_have(domain,'r_inner');
    must_have(domain,'r_outer');

    coord = mesh.coord;
    connect = mesh.connect;

    if size(connect,2) ~= 6
        error('SIF_LEFM_interaction_EDI:BadConnect', ...
            'mesh.connect must contain T6 elements.');
    end
    if numel(U) ~= 2*size(coord,1)
        error('SIF_LEFM_interaction_EDI:BadU', ...
            'numel(U) must equal 2*size(mesh.coord,1).');
    end

    E = mat.E;
    nu = mat.nu;

    if isfield(mat,'Dmat') && ~isempty(mat.Dmat)
        Dmat = mat.Dmat;
    elseif isfield(mat,'D') && ~isempty(mat.D)
        Dmat = mat.D;
    else
        error('SIF_LEFM_interaction_EDI:MissingD', ...
            'mat must contain Dmat or D.');
    end

    if isempty(ip.Results.UsePlaneStrain)
        if isfield(mat,'ps') && ~isempty(mat.ps)
            planeStrain = logical(mat.ps == 1);
        else
            planeStrain = true;
        end
    else
        planeStrain = logical(ip.Results.UsePlaneStrain);
    end

    if planeStrain
        Eeff = E/(1-nu^2);
    else
        Eeff = E;
    end

    r_inner = domain.r_inner;
    r_outer = domain.r_outer;
    if ~(isfinite(r_inner) && isfinite(r_outer) && ...
            r_inner >= 0 && r_outer > r_inner)
        error('SIF_LEFM_interaction_EDI:BadDomain', ...
            'Require 0 <= r_inner < r_outer.');
    end

    if size(V,1) < 2 || size(V,2) ~= 2
        error('SIF_LEFM_interaction_EDI:BadCrackPolyline', ...
            'V must be [nPts x 2], nPts >= 2.');
    end

    x_tip = V(end,:).';
    e1 = x_tip - V(end-1,:).';
    if norm(e1) <= eps
        error('SIF_LEFM_interaction_EDI:DegenerateLastSegment', ...
            'The final crack segment has zero length.');
    end
    e1 = e1/norm(e1);
    e2 = [-e1(2); e1(1)];
    R_gl = [e1,e2];
    R_loc = R_gl.';

    [nip,xip,w] = local_integr_T6();

    I_modeI = 0;
    I_modeII = 0;
    nGP_total = 0;
    nGP_used = 0;
    usedElem = false(size(connect,1),1);

    rows = zeros(0,10);
    auxMismatchI = zeros(0,1);
    auxMismatchII = zeros(0,1);

    matAux = mat;
    matAux.Dmat = Dmat;
    matAux.ps = double(planeStrain);

    for e = 1:size(connect,1)
        nodes = connect(e,:);
        X = coord(nodes,:);

        % Cheap rejection based on the vertex centroid and element radius.
        xc = mean(X(1:3,:),1).';
        rc = norm(R_loc*(xc-x_tip));
        if rc > r_outer + local_element_radius(X)
            continue;
        end

        uel = [U(2*nodes-1), U(2*nodes)];

        for ig = 1:nip
            nGP_total = nGP_total + 1;
            xi = xip(:,ig);

            [N,DetJ,dNdx] = local_T6_shape_grad(xi,X);
            if DetJ <= 0
                error('SIF_LEFM_interaction_EDI:BadElement', ...
                    'Inverted or degenerate T6 element e=%d.', e);
            end

            xg = (N(:).'*X).';
            xl = R_loc*(xg-x_tip);
            r = hypot(xl(1),xl(2));

            if r <= r_inner || r >= r_outer || r <= 1e-14
                continue;
            end

            nGP_used = nGP_used + 1;
            usedElem(e) = true;

            % Actual FEM field.
            Grad_gl = [ ...
                dNdx(1,:)*uel(:,1), dNdx(2,:)*uel(:,1); ...
                dNdx(1,:)*uel(:,2), dNdx(2,:)*uel(:,2)];
            Grad1 = R_loc*Grad_gl*R_gl;

            eps1 = [Grad1(1,1); Grad1(2,2); Grad1(1,2)+Grad1(2,1)];
            sig1 = Dmat*eps1;
            du1_dx1 = Grad1(:,1);

            qgrad = local_qgrad_radial(xl,r_inner,r_outer);

            auxI = SIF_LEFM_auxiliary_fields( ...
                xl(1),xl(2),Kaux,0,matAux, ...
                'UsePlaneStrain',planeStrain, ...
                'DerivativeMode',ip.Results.AuxDerivativeMode, ...
                'FDRelStep',ip.Results.AuxFDRelStep, ...
                'FDAbsStep',ip.Results.AuxFDAbsStep);

            auxII = SIF_LEFM_auxiliary_fields( ...
                xl(1),xl(2),0,Kaux,matAux, ...
                'UsePlaneStrain',planeStrain, ...
                'DerivativeMode',ip.Results.AuxDerivativeMode, ...
                'FDRelStep',ip.Results.AuxFDRelStep, ...
                'FDAbsStep',ip.Results.AuxFDAbsStep);

            densI = local_interaction_density(sig1,eps1,du1_dx1,auxI,qgrad);
            densII = local_interaction_density(sig1,eps1,du1_dx1,auxII,qgrad);

            dA = w(ig)*(DetJ/2);
            I_modeI = I_modeI + densI*dA;
            I_modeII = I_modeII + densII*dA;

            auxMismatchI(end+1,1) = auxI.constitutive_mismatch; %#ok<AGROW>
            auxMismatchII(end+1,1) = auxII.constitutive_mismatch; %#ok<AGROW>
            rows(end+1,:) = [e,ig,xl(1),xl(2),r,densI,densII,dA, ...
                auxI.constitutive_mismatch,auxII.constitutive_mismatch]; %#ok<AGROW>
        end
    end

    % IMPORTANT: interaction-integral normalization contains the factor 2.
    KI = 0.5*Eeff*I_modeI/Kaux;
    KII = 0.5*Eeff*I_modeII/Kaux;

    Aux = struct();
    Aux.method = 'interaction_EDI';
    Aux.normalization = 'K=(Eeff/2)*I/Kaux';
    Aux.KI = KI;
    Aux.KII = KII;
    Aux.I_modeI = I_modeI;
    Aux.I_modeII = I_modeII;
    Aux.Kaux = Kaux;
    Aux.Eeff = Eeff;
    Aux.planeStrain = planeStrain;
    Aux.auxDerivativeMode = char(ip.Results.AuxDerivativeMode);
    Aux.r_inner = r_inner;
    Aux.r_outer = r_outer;
    Aux.x_tip = x_tip.';
    Aux.e1 = e1.';
    Aux.e2 = e2.';
    Aux.R_gl = R_gl;
    Aux.R_loc = R_loc;
    Aux.nGP_total = nGP_total;
    Aux.nGP_used = nGP_used;
    Aux.nElem_used = nnz(usedElem);

    if isempty(auxMismatchI)
        Aux.auxConsistency = struct('maxI',NaN,'maxII',NaN, ...
            'meanI',NaN,'meanII',NaN);
    else
        Aux.auxConsistency = struct( ...
            'maxI',max(auxMismatchI), ...
            'maxII',max(auxMismatchII), ...
            'meanI',mean(auxMismatchI), ...
            'meanII',mean(auxMismatchII));
    end

    if isempty(rows)
        Aux.gpTable = table();
    else
        Aux.gpTable = array2table(rows, 'VariableNames', { ...
            'elem','igp','x1','x2','r','densI','densII','dA', ...
            'auxMismatchI','auxMismatchII'});
    end

    if verbose
        fprintf('\nSIF_LEFM_interaction_EDI summary:\n');
        fprintf('  domain r/a-unscaled = [%.8e, %.8e]\n',r_inner,r_outer);
        fprintf('  derivative mode     = %s\n',Aux.auxDerivativeMode);
        fprintf('  I_modeI / I_modeII  = %.10e / %.10e\n',I_modeI,I_modeII);
        fprintf('  KI / KII            = %.10e / %.10e\n',KI,KII);
        fprintf('  KII/KI              = %.10e\n',KII/KI);
        fprintf('  GP used / total     = %d / %d\n',nGP_used,nGP_total);
        fprintf('  elements used       = %d\n',Aux.nElem_used);
        fprintf('  aux mismatch max I/II = %.3e / %.3e\n', ...
            Aux.auxConsistency.maxI,Aux.auxConsistency.maxII);
    end
end


function dens = local_interaction_density(sig1,eps1,du1_dx1,aux,qgrad)
% Interaction flux with the sign chosen for q=1 at the inner boundary and
% q=0 at the outer boundary:
%
% A_j = -Wint*delta_1j + sigma1_ij*u2_i,1 + sigma2_ij*u1_i,1
% I   = integral_A A_j q_,j dA.

    sig2 = aux.sig;
    du2_dx1 = aux.du_dx1;

    Wint = sig2(1)*eps1(1) + sig2(2)*eps1(2) + sig2(3)*eps1(3);

    S1 = [sig1(1),sig1(3); sig1(3),sig1(2)];
    S2 = [sig2(1),sig2(3); sig2(3),sig2(2)];

    A = [ ...
        -Wint + dot(S1(:,1),du2_dx1) + dot(S2(:,1),du1_dx1); ...
                 dot(S1(:,2),du2_dx1) + dot(S2(:,2),du1_dx1)];

    dens = A.'*qgrad(:);
end


function qgrad = local_qgrad_radial(xl,r_inner,r_outer)
    r = hypot(xl(1),xl(2));
    if r <= r_inner || r >= r_outer || r <= eps
        qgrad = [0;0];
        return;
    end
    qgrad = -(xl(:)/r)/(r_outer-r_inner);
end


function [N,DetJ,dNdx] = local_T6_shape_grad(xi,X)
    L1 = xi(1);
    L2 = xi(2);
    L3 = 1-L1-L2;

    N = [ ...
        L1*(2*L1-1); ...
        L2*(2*L2-1); ...
        L3*(2*L3-1); ...
        4*L1*L2; ...
        4*L2*L3; ...
        4*L3*L1];

    dN_dL1 = [4*L1-1;0;-(4*L3-1);4*L2;-4*L2;4*(L3-L1)];
    dN_dL2 = [0;4*L2-1;-(4*L3-1);4*L1;4*(L3-L2);-4*L1];

    dNdxi = [dN_dL1.';dN_dL2.'];
    J = dNdxi*X;
    DetJ = det(J);
    dNdx = J\dNdxi;
end


function [nip,xip,w] = local_integr_T6()
% Seven-point Dunavant rule. Weights sum to one; hence dA=w*DetJ/2.
    nip = 7;
    xip = zeros(2,nip);
    w = zeros(1,nip);

    xip(:,1) = [1/3;1/3];
    w(1) = 0.225;

    a1 = 0.059715871789770;
    b1 = 0.470142064105115;
    w1 = 0.132394152788506;
    xip(:,2) = [a1;b1];
    xip(:,3) = [b1;a1];
    xip(:,4) = [b1;b1];
    w(2:4) = w1;

    a2 = 0.797426985353087;
    b2 = 0.101286507323456;
    w2 = 0.125939180544827;
    xip(:,5) = [a2;b2];
    xip(:,6) = [b2;a2];
    xip(:,7) = [b2;b2];
    w(5:7) = w2;
end


function rad = local_element_radius(X)
    xc = mean(X(1:3,:),1);
    rad = max(sqrt(sum((X-xc).^2,2)));
end


function must_have(S,f)
    if ~isstruct(S) || ~isfield(S,f) || isempty(S.(f))
        error('SIF_LEFM_interaction_EDI:MissingField', ...
            'Required field "%s" is missing or empty.',f);
    end
end
