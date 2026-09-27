function Out = validate_EDI_Williams_fields(varargin)
%VALIDATE_EDI_WILLIAMS_FIELDS
% Independent synthetic-field check of SIF_LEFM_interaction_EDI.
%
% A polar annulus is meshed directly (no PDE Toolbox, no crack solve).
% The two crack faces theta=-pi and theta=+pi have distinct node IDs.
% Exact leading-order Williams displacements are prescribed at all T6
% nodes, and the existing interaction-EDI extractor is then asked to
% recover the imposed KI/KII values.
%
% This test is intended to diagnose:
%   1) interaction-integral normalization;
%   2) mode-II sign convention;
%   3) cross-mode leakage;
%   4) convergence with mesh refinement.
%
% IMPORTANT:
%   The function does NOT modify or compensate the production EDI result.
%   Both the raw recovery ratio and the diagnostic 0.5*raw ratio are
%   reported so a possible factor-of-two normalization error is visible.
%
% Usage:
%   O = validate_EDI_Williams_fields();
%
% Name-value options:
%   'E'          Young modulus, default 4e3
%   'nu'         Poisson ratio, default 0.30
%   'ps'         1 plane strain, 0 plane stress, default 1
%   'NrList'     radial element counts, default [8 16]
%   'NthList'    angular element counts, default [64 128]
%   'rMeshInner' mesh inner radius, default 0.01
%   'rMeshOuter' mesh outer radius, default 0.20
%   'rInner'     EDI q=1 radius, default 0.04
%   'rOuter'     EDI q=0 radius, default 0.16
%   'Verbose'    print table, default true

    ip = inputParser;
    addParameter(ip, 'E', 4e3, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'nu', 0.30, @(x)isnumeric(x) && isscalar(x) && x>0 && x<0.5);
    addParameter(ip, 'ps', 1, @(x)isnumeric(x) && isscalar(x) && any(x==[0 1]));
    addParameter(ip, 'NrList', [8 16], @(x)isnumeric(x) && isvector(x) && all(x>=2));
    addParameter(ip, 'NthList', [64 128], @(x)isnumeric(x) && isvector(x) && all(x>=8));
    addParameter(ip, 'rMeshInner', 0.01, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'rMeshOuter', 0.20, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'rInner', 0.04, @(x)isnumeric(x) && isscalar(x) && x>=0);
    addParameter(ip, 'rOuter', 0.16, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'Verbose', true, @(x)islogical(x) || isnumeric(x));
    parse(ip, varargin{:});
    S = ip.Results;

    NrList = round(S.NrList(:));
    NthList = round(S.NthList(:));

    if numel(NrList) ~= numel(NthList)
        error('validate_EDI_Williams_fields:MeshListSize', ...
            'NrList and NthList must have equal length.');
    end

    if ~(S.rMeshInner < S.rInner && S.rInner < S.rOuter && S.rOuter < S.rMeshOuter)
        error('validate_EDI_Williams_fields:BadRadii', ...
            'Require rMeshInner < rInner < rOuter < rMeshOuter.');
    end

    E = S.E;
    nu = S.nu;
    ps = S.ps;

    if ps == 1
        coef = E/((1+nu)*(1-2*nu));
        D = coef * [ ...
            1-nu, nu, 0; ...
            nu, 1-nu, 0; ...
            0, 0, (1-2*nu)/2 ];
        kappa = 3 - 4*nu;
    else
        coef = E/(1-nu^2);
        D = coef * [ ...
            1, nu, 0; ...
            nu, 1, 0; ...
            0, 0, (1-nu)/2 ];
        kappa = (3-nu)/(1+nu);
    end

    mu = E/(2*(1+nu));

    mat = struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);

    % Tip at origin; previous crack point on negative x-axis makes the
    % local crack direction e1 = +x.
    V = [-1,0; 0,0];
    domain = struct('r_inner',S.rInner,'r_outer',S.rOuter);

    % [KI, KII] imposed exact fields.
    Kcases = [ ...
        1.0, 0.0; ...
        0.0, 1.0; ...
        1.0, 0.35 ];

    caseName = ["pure_I"; "pure_II"; "mixed_I_II"];

    rows = [];
    details = cell(numel(NrList), size(Kcases,1));

    for im = 1:numel(NrList)
        Nr = NrList(im);
        Nth = NthList(im);

        mesh = build_polar_crack_annulus( ...
            S.rMeshInner, S.rMeshOuter, Nr, Nth);

        for ic = 1:size(Kcases,1)
            KIin = Kcases(ic,1);
            KIIin = Kcases(ic,2);

            U = exact_williams_displacement_vector( ...
                mesh.coord, KIin, KIIin, mu, kappa);

            [KIraw, KIIraw, Aux] = SIF_LEFM_interaction_EDI( ...
                mesh, U, V, mat, domain, ...
                'UsePlaneStrain', ps==1, ...
                'Verbose', false);

            if abs(KIin) > 0
                ratioKI = KIraw/KIin;
                ratioKIhalf = 0.5*KIraw/KIin;
            else
                ratioKI = NaN;
                ratioKIhalf = NaN;
            end

            if abs(KIIin) > 0
                ratioKII = KIIraw/KIIin;
                ratioKIIhalf = 0.5*KIIraw/KIIin;
            else
                ratioKII = NaN;
                ratioKIIhalf = NaN;
            end

            rows = [rows; ...
                im, ic, Nr, Nth, KIin, KIIin, KIraw, KIIraw, ...
                ratioKI, ratioKII, ratioKIhalf, ratioKIIhalf, ...
                Aux.nElem_used, Aux.nGP_used]; %#ok<AGROW>

            details{im,ic} = Aux;
        end
    end

    T = array2table(rows, 'VariableNames', { ...
        'meshLevel','caseID','Nr','Nth','KI_input','KII_input', ...
        'KI_raw','KII_raw','KI_raw_over_input','KII_raw_over_input', ...
        'half_KI_raw_over_input','half_KII_raw_over_input', ...
        'nElem_used','nGP_used'});

    % Add readable case names without relying on categorical ordering.
    T.caseName = caseName(T.caseID);

    % Reorder for console readability.
    T = movevars(T, 'caseName', 'After', 'caseID');

    Out = struct();
    Out.table = T;
    Out.details = details;
    Out.settings = S;
    Out.material = mat;
    Out.domain = domain;
    Out.Kcases = Kcases;
    Out.caseName = caseName;

    if logical(S.Verbose)
        fprintf('\n============================================================\n');
        fprintf('EDI SYNTHETIC WILLIAMS-FIELD VALIDATION\n');
        fprintf('============================================================\n');
        fprintf('plane state     : %s\n', ternary(ps==1,'plane strain','plane stress'));
        fprintf('EDI annulus     : [%.6g, %.6g]\n', S.rInner, S.rOuter);
        fprintf('mesh annulus    : [%.6g, %.6g]\n\n', S.rMeshInner, S.rMeshOuter);

        disp(T(:, { ...
            'meshLevel','caseName','Nr','Nth', ...
            'KI_input','KII_input','KI_raw','KII_raw', ...
            'KI_raw_over_input','KII_raw_over_input', ...
            'half_KI_raw_over_input','half_KII_raw_over_input'}));

        fprintf(['\nInterpretation rule:\n', ...
            '  correct raw normalization -> raw/input tends to 1;\n', ...
            '  missing factor 1/2       -> raw/input tends to 2 while ', ...
            '0.5*raw/input tends to 1.\n', ...
            'For pure modes, the non-imposed mode should converge to zero.\n']);
    end
end


% =========================================================================
function mesh = build_polar_crack_annulus(r0, r1, Nr, Nth)
% Structured T3 annulus with duplicated theta=-pi/+pi crack-face nodes.

    rv = linspace(r0, r1, Nr+1);
    tv = linspace(-pi, pi, Nth+1);

    nR = numel(rv);
    nT = numel(tv);

    coord3 = zeros(nR*nT,2);

    id = @(ir,it) (it-1)*nR + ir;

    for it = 1:nT
        th = tv(it);
        for ir = 1:nR
            r = rv(ir);
            coord3(id(ir,it),:) = r*[cos(th), sin(th)];
        end
    end

    connect3 = zeros(2*Nr*Nth,3);
    e = 0;

    for it = 1:Nth
        for ir = 1:Nr
            A = id(ir,   it);
            B = id(ir+1, it);
            C = id(ir+1, it+1);
            D = id(ir,   it+1);

            e = e+1;
            connect3(e,:) = [A B C];

            e = e+1;
            connect3(e,:) = [A C D];
        end
    end

    % Numerical safety: enforce CCW.
    area = tri_area_signed(connect3, coord3);
    cw = area < 0;
    if any(cw)
        tmp = connect3(cw,2);
        connect3(cw,2) = connect3(cw,3);
        connect3(cw,3) = tmp;
    end

    [coord6, connect6] = T3toT6_fast(coord3, connect3);

    mesh = struct();
    mesh.coord3 = coord3;
    mesh.connect3 = connect3;
    mesh.coord = coord6;
    mesh.connect = connect6;
end


function U = exact_williams_displacement_vector(coord, KI, KII, mu, kappa)
% Leading-order isotropic Williams displacement field in local Cartesian
% crack coordinates. The crack lies on x<0 and faces are theta=+/-pi.

    n = size(coord,1);
    U = zeros(2*n,1);

    for i = 1:n
        x = coord(i,1);
        y = coord(i,2);

        r = hypot(x,y);
        th = atan2(y,x);

        fac = sqrt(r/(2*pi))/(2*mu);
        c = cos(th/2);
        s = sin(th/2);

        u1I = KI * fac * c * (kappa - 1 + 2*s^2);
        u2I = KI * fac * s * (kappa + 1 - 2*c^2);

        u1II = KII * fac * s * (kappa + 1 + 2*c^2);
        u2II = -KII * fac * c * (kappa - 1 - 2*s^2);

        U(2*i-1) = u1I + u1II;
        U(2*i)   = u2I + u2II;
    end
end


function A = tri_area_signed(T, X)
    v1 = X(T(:,1),:);
    v2 = X(T(:,2),:);
    v3 = X(T(:,3),:);

    A = 0.5*((v2(:,1)-v1(:,1)).*(v3(:,2)-v1(:,2)) - ...
             (v2(:,2)-v1(:,2)).*(v3(:,1)-v1(:,1)));
end


function out = ternary(tf,a,b)
    if tf
        out = a;
    else
        out = b;
    end
end
