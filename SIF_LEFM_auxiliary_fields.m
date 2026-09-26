function aux = SIF_LEFM_auxiliary_fields(x1, x2, KI, KII, mat, varargin)
%SIF_LEFM_AUXILIARY_FIELDS  Analytical 2D isotropic crack-tip auxiliary fields.
%
% The local frame is x1 along the crack and x2 normal to the crack.  The
% sign convention is the one used by SIF_LEFM_interaction_EDI.
%
% Inputs
%   x1,x2   local crack-tip coordinates
%   KI,KII  auxiliary stress-intensity amplitudes
%   mat     struct with E, nu and Dmat (or D)
%
% Name-value options
%   'UsePlaneStrain'  [] (infer from mat.ps, otherwise true)
%   'DerivativeMode'  'analytic' (default) or 'finite_difference'
%   'FDRelStep'       1e-5
%   'FDAbsStep'       1e-10
%
% Outputs in aux
%   sig              [sigma11;sigma22;sigma12]
%   eps              [eps11;eps22;gamma12] from displacement gradient
%   eps_from_sig     D\sig
%   GradU            2x2 displacement gradient
%   du_dx1           first column of GradU
%   constitutive_mismatch
%
% The analytical gradient is obtained from u_i = C sqrt(r) f_i(theta):
%   u_i,1 = C/sqrt(r) [0.5 cos(theta) f_i - sin(theta) f_i']
%   u_i,2 = C/sqrt(r) [0.5 sin(theta) f_i + cos(theta) f_i'].

    ip = inputParser;
    addParameter(ip, 'UsePlaneStrain', [], @(x)islogical(x) || isnumeric(x) || isempty(x));
    addParameter(ip, 'DerivativeMode', 'analytic', @(x)ischar(x) || isstring(x));
    addParameter(ip, 'FDRelStep', 1e-5, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'FDAbsStep', 1e-10, @(x)isnumeric(x) && isscalar(x) && x>0);
    parse(ip, varargin{:});

    must_have(mat, 'E');
    must_have(mat, 'nu');

    E = mat.E;
    nu = mat.nu;

    if isfield(mat,'Dmat') && ~isempty(mat.Dmat)
        Dmat = mat.Dmat;
    elseif isfield(mat,'D') && ~isempty(mat.D)
        Dmat = mat.D;
    else
        error('SIF_LEFM_auxiliary_fields:MissingD', ...
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

    mu = E/(2*(1+nu));
    if planeStrain
        kappa = 3 - 4*nu;
    else
        kappa = (3-nu)/(1+nu);
    end

    r = hypot(x1,x2);
    if r <= 1e-14
        aux = zero_aux(mu,kappa,planeStrain);
        return;
    end

    th = atan2(x2,x1);

    sig = local_stress(r,th,KI,KII);

    mode = lower(char(ip.Results.DerivativeMode));
    switch mode
        case 'analytic'
            GradU = local_grad_analytic(r,th,KI,KII,mu,kappa);
            h = NaN;
        case {'finite_difference','fd'}
            h = max(ip.Results.FDAbsStep, ip.Results.FDRelStep*r);
            up = local_displacement(x1+h,x2,KI,KII,mu,kappa);
            um = local_displacement(x1-h,x2,KI,KII,mu,kappa);
            vp = local_displacement(x1,x2+h,KI,KII,mu,kappa);
            vm = local_displacement(x1,x2-h,KI,KII,mu,kappa);
            GradU = [(up-um)/(2*h), (vp-vm)/(2*h)];
        otherwise
            error('SIF_LEFM_auxiliary_fields:BadDerivativeMode', ...
                'DerivativeMode must be analytic or finite_difference.');
    end

    eps_grad = [GradU(1,1); GradU(2,2); GradU(1,2)+GradU(2,1)];
    eps_sig  = Dmat \ sig;

    denom = max(norm(eps_sig), eps);
    mismatch = norm(eps_grad-eps_sig)/denom;

    aux = struct();
    aux.sig = sig;
    aux.eps = eps_grad;
    aux.eps_from_sig = eps_sig;
    aux.GradU = GradU;
    aux.du_dx1 = GradU(:,1);
    aux.constitutive_mismatch = mismatch;
    aux.derivative_mode = mode;
    aux.fd_step = h;
    aux.r = r;
    aux.theta = th;
    aux.mu = mu;
    aux.kappa = kappa;
    aux.planeStrain = planeStrain;
end


function sig = local_stress(r,th,KI,KII)
    c  = cos(th/2);
    s  = sin(th/2);
    c3 = cos(3*th/2);
    s3 = sin(3*th/2);
    fac = 1/sqrt(2*pi*r);

    s11I = KI*fac*c*(1-s*s3);
    s22I = KI*fac*c*(1+s*s3);
    s12I = KI*fac*c*s*c3;

    s11II = -KII*fac*s*(2+c*c3);
    s22II =  KII*fac*s*c*c3;
    s12II =  KII*fac*c*(1-s*s3);

    sig = [s11I+s11II; s22I+s22II; s12I+s12II];
end


function GradU = local_grad_analytic(r,th,KI,KII,mu,kappa)
    c = cos(th/2);
    s = sin(th/2);
    ct = cos(th);
    st = sin(th);

    % Mode I angular functions.
    f1I = c*(kappa-ct);
    f2I = s*(kappa-ct);
    f1Ip = -0.5*s*(kappa-ct) + c*st;
    f2Ip =  0.5*c*(kappa-ct) + s*st;

    % Mode II angular functions.
    f1II = s*(kappa+2+ct);
    f2II = -c*(kappa-2+ct);
    f1IIp = 0.5*c*(kappa+2+ct) - s*st;
    f2IIp = 0.5*s*(kappa-2+ct) + c*st;

    f1  = KI*f1I  + KII*f1II;
    f2  = KI*f2I  + KII*f2II;
    f1p = KI*f1Ip + KII*f1IIp;
    f2p = KI*f2Ip + KII*f2IIp;

    pref = 1/(2*mu*sqrt(2*pi*r));

    du1dx1 = pref*(0.5*ct*f1 - st*f1p);
    du2dx1 = pref*(0.5*ct*f2 - st*f2p);

    du1dx2 = pref*(0.5*st*f1 + ct*f1p);
    du2dx2 = pref*(0.5*st*f2 + ct*f2p);

    GradU = [du1dx1,du1dx2; du2dx1,du2dx2];
end


function u = local_displacement(x1,x2,KI,KII,mu,kappa)
    r = hypot(x1,x2);
    th = atan2(x2,x1);

    if r <= 1e-14
        u = [0;0];
        return;
    end

    fac = sqrt(r/(2*pi))/(2*mu);
    c = cos(th/2);
    s = sin(th/2);

    u1I  = KI*fac*c*(kappa-1+2*s^2);
    u2I  = KI*fac*s*(kappa+1-2*c^2);
    u1II = KII*fac*s*(kappa+1+2*c^2);
    u2II = -KII*fac*c*(kappa-1-2*s^2);

    u = [u1I+u1II; u2I+u2II];
end


function A = zero_aux(mu,kappa,planeStrain)
    A = struct('sig',[0;0;0], 'eps',[0;0;0], ...
        'eps_from_sig',[0;0;0], 'GradU',zeros(2), ...
        'du_dx1',[0;0], 'constitutive_mismatch',0, ...
        'derivative_mode','analytic', 'fd_step',NaN, ...
        'r',0, 'theta',0, 'mu',mu, 'kappa',kappa, ...
        'planeStrain',planeStrain);
end


function must_have(S,f)
    if ~isstruct(S) || ~isfield(S,f) || isempty(S.(f))
        error('SIF_LEFM_auxiliary_fields:MissingField', ...
            'Required field "%s" is missing.', f);
    end
end
