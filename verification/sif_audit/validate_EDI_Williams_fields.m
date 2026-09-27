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
%   The production EDI normalization has now been corrected by the factor
%   1/2 established by the first synthetic run. This routine is therefore
%   a regression/convergence check: recovered/input should tend to 1.
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
%   'AssertFine'  assert fine-mesh recovery within tolerance, default true
%   'FineTol'     fine-mesh relative tolerance, default 0.02
%   'MeshTopology' 'mirror_reflected' (default) or 'legacy_same_diagonal'

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
    addParameter(ip, 'AssertFine', true, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'FineTol', 0.02, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'MeshTopology', 'mirror_reflected', ...
        @(x)ischar(x) || (isstring(x) && isscalar(x)));
    parse(ip, varargin{:});
    S = ip.Results;

    NrList = round(S.NrList(:));
    NthList = round(S.NthList(:));
    meshTopology = char(S.MeshTopology);

    validTopology = strcmpi(meshTopology,'mirror_reflected') || ...
                    strcmpi(meshTopology,'legacy_same_diagonal');
    if ~validTopology
        error('validate_EDI_Williams_fields:BadMeshTopology', ...
            'MeshTopology must be mirror_reflected or legacy_same_diagonal.');
    end

    if strcmpi(meshTopology,'mirror_reflected') && any(mod(NthList,2)~=0)
        error('validate_EDI_Williams_fields:OddNth', ...
            'mirror_reflected topology requires every Nth value to be even.');
    end

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
    meshAudit = cell(numel(NrList),1);

    for im = 1:numel(NrList)
        Nr = NrList(im);
        Nth = NthList(im);

        mesh = build_polar_crack_annulus( ...
            S.rMeshInner, S.rMeshOuter, Nr, Nth, meshTopology);

        meshAudit{im} = mesh.audit;

        for ic = 1:size(Kcases,1)
            KIin = Kcases(ic,1);
            KIIin = Kcases(ic,2);

            U = exact_williams_displacement_vector( ...
                mesh.coord, KIin, KIIin, mu, kappa);

            [KIrec, KIIrec, Aux] = SIF_LEFM_interaction_EDI( ...
                mesh, U, V, mat, domain, ...
                'UsePlaneStrain', ps==1, ...
                'Verbose', false);

            if abs(KIin) > 0
                ratioKI = KIrec/KIin;
                relErrKI = abs(ratioKI - 1);
            else
                ratioKI = NaN;
                relErrKI = abs(KIrec);
            end

            if abs(KIIin) > 0
                ratioKII = KIIrec/KIIin;
                relErrKII = abs(ratioKII - 1);
            else
                ratioKII = NaN;
                relErrKII = abs(KIIrec);
            end

            rows = [rows; ...
                im, ic, Nr, Nth, KIin, KIIin, KIrec, KIIrec, ...
                ratioKI, ratioKII, relErrKI, relErrKII, ...
                Aux.nElem_used, Aux.nGP_used]; %#ok<AGROW>

            details{im,ic} = Aux;
        end
    end

    T = array2table(rows, 'VariableNames', { ...
        'meshLevel','caseID','Nr','Nth','KI_input','KII_input', ...
        'KI_recovered','KII_recovered','KI_recovered_over_input','KII_recovered_over_input', ...
        'KI_error_metric','KII_error_metric', ...
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
    Out.meshAudit = meshAudit;
    Out.meshTopology = meshTopology;

    if logical(S.Verbose)
        fprintf('\n============================================================\n');
        fprintf('EDI SYNTHETIC WILLIAMS-FIELD VALIDATION\n');
        fprintf('============================================================\n');
        fprintf('plane state     : %s\n', ternary(ps==1,'plane strain','plane stress'));
        fprintf('EDI annulus     : [%.6g, %.6g]\n', S.rInner, S.rOuter);
        fprintf('mesh annulus    : [%.6g, %.6g]\n', S.rMeshInner, S.rMeshOuter);
        fprintf('mesh topology   : %s\n', meshTopology);
        if strcmpi(meshTopology,'mirror_reflected')
            fprintf('mirror coord err: %.3e\n\n', meshAudit{end}.maxMirrorCoordError);
        else
            fprintf('\n');
        end

        disp(T(:, { ...
            'meshLevel','caseName','Nr','Nth', ...
            'KI_input','KII_input','KI_recovered','KII_recovered', ...
            'KI_recovered_over_input','KII_recovered_over_input', ...
            'KI_error_metric','KII_error_metric'}));

        fprintf(['\nInterpretation rule:\n', ...
            '  imposed-mode recovered/input should tend to 1;\n', ...
            '  non-imposed pure-mode component should tend to zero;\n', ...
            '  pure mode II should recover with positive sign.\n']);
    end
    % Fine-mesh regression gate.
    fine = T(T.meshLevel == max(T.meshLevel), :);

    tol = S.FineTol;

    pureI = fine(fine.caseName == "pure_I", :);
    pureII = fine(fine.caseName == "pure_II", :);
    mixed = fine(fine.caseName == "mixed_I_II", :);

    checks = struct();
    checks.pureI_ratio = pureI.KI_recovered_over_input;
    checks.pureI_cross = abs(pureI.KII_recovered);
    checks.pureII_ratio = pureII.KII_recovered_over_input;
    checks.pureII_cross = abs(pureII.KI_recovered);
    checks.mixed_KI_ratio = mixed.KI_recovered_over_input;
    checks.mixed_KII_ratio = mixed.KII_recovered_over_input;

    scale = max([1; abs(fine.KI_input); abs(fine.KII_input)]);

    checks.pass = ...
        abs(checks.pureI_ratio - 1) <= tol && ...
        abs(checks.pureII_ratio - 1) <= tol && ...
        checks.pureI_cross <= tol*scale && ...
        checks.pureII_cross <= tol*scale && ...
        abs(checks.mixed_KI_ratio - 1) <= tol && ...
        abs(checks.mixed_KII_ratio - 1) <= tol;

    Out.checks = checks;
    Out.fineTolerance = tol;

    if logical(S.AssertFine) && ~checks.pass
        error('validate_EDI_Williams_fields:RegressionFailed', ...
            'Fine-mesh Williams-field recovery failed the %.3g tolerance.', tol);
    end
end


% =========================================================================
function mesh = build_polar_crack_annulus(r0, r1, Nr, Nth, topology)
% Structured polar crack annulus.
%
% topology = 'mirror_reflected'
%   Build only the upper half (0 <= theta <= pi), then reflect the complete
%   T3 mesh across x2=0. Reflection maps both node coordinates and element
%   connectivity. Triangle orientation is reversed after reflection so all
%   elements remain CCW. The theta=0 radial line is shared, while the
%   theta=+pi and theta=-pi crack faces use distinct node IDs.
%
% topology = 'legacy_same_diagonal'
%   Historical synthetic builder: generate theta=-pi..pi directly and use
%   the same logical A-C diagonal in every polar quadrilateral.

    if strcmpi(topology,'mirror_reflected')
        mesh = build_mirror_reflected_annulus(r0,r1,Nr,Nth);
    else
        mesh = build_legacy_same_diagonal_annulus(r0,r1,Nr,Nth);
    end
end


function mesh = build_mirror_reflected_annulus(r0,r1,Nr,Nth)

    if mod(Nth,2)~=0
        error('validate_EDI_Williams_fields:OddNthInternal', ...
            'Mirror-reflected annulus requires even Nth.');
    end

    rv = linspace(r0,r1,Nr+1);
    nR = numel(rv);
    Nh = Nth/2;

    % Upper half only: theta = 0 ... pi.
    tvU = linspace(0,pi,Nh+1);
    idU = zeros(nR,Nh+1);
    coord3 = zeros(nR*(Nh+1),2);

    id = 0;
    for jt = 1:Nh+1
        th = tvU(jt);
        for ir = 1:nR
            id = id+1;
            idU(ir,jt)=id;
            coord3(id,:) = rv(ir)*[cos(th),sin(th)];
        end
    end

    % Upper T3 cells use one consistent A-C diagonal.
    Tup = zeros(2*Nr*Nh,3);
    e = 0;
    for jt = 1:Nh
        for ir = 1:Nr
            A = idU(ir,  jt);
            B = idU(ir+1,jt);
            C = idU(ir+1,jt+1);
            D = idU(ir,  jt+1);

            e=e+1; Tup(e,:)=[A B C];
            e=e+1; Tup(e,:)=[A C D];
        end
    end

    % Enforce CCW on the upper half before reflection.
    areaU = tri_area_signed(Tup,coord3);
    cw = areaU<0;
    if any(cw)
        tmp=Tup(cw,2);
        Tup(cw,2)=Tup(cw,3);
        Tup(cw,3)=tmp;
    end

    % Mirror map. theta=0 nodes are shared. All theta>0 nodes, including
    % theta=pi crack-face nodes, receive distinct reflected node IDs.
    mirrorMap = zeros(size(coord3,1),1);
    mirrorMap(idU(:,1)) = idU(:,1);

    originalUpperCount = size(coord3,1);
    next = originalUpperCount;

    for jt = 2:Nh+1
        for ir = 1:nR
            iu = idU(ir,jt);
            next = next+1;
            coord3(next,:) = [coord3(iu,1), -coord3(iu,2)]; %#ok<AGROW>
            mirrorMap(iu)=next;
        end
    end

    % Reflect upper connectivity. A geometric reflection changes triangle
    % orientation, so swap local vertices 2 and 3 to restore CCW.
    Tlo = mirrorMap(Tup);
    Tlo = Tlo(:,[1 3 2]);

    connect3 = [Tup;Tlo];

    % Safety checks.
    area = tri_area_signed(connect3,coord3);
    if any(area<=0)
        error('validate_EDI_Williams_fields:MirrorMeshOrientation', ...
            'Mirror-reflected mesh contains non-positive T3 area.');
    end

    % Coordinate reflection must be exact to roundoff for every upper node.
    upperIDs = (1:originalUpperCount).';
    mapped = mirrorMap(upperIDs);
    reflected = [coord3(upperIDs,1), -coord3(upperIDs,2)];
    mirrorErr = sqrt(sum((coord3(mapped,:)-reflected).^2,2));

    [coord6,connect6]=T3toT6_fast(coord3,connect3);

    mesh=struct();
    mesh.coord3=coord3;
    mesh.connect3=connect3;
    mesh.coord=coord6;
    mesh.connect=connect6;

    mesh.audit=struct();
    mesh.audit.topology='mirror_reflected';
    mesh.audit.Nr=Nr;
    mesh.audit.Nth=Nth;
    mesh.audit.nUpperT3=size(Tup,1);
    mesh.audit.nLowerT3=size(Tlo,1);
    mesh.audit.nT3=size(connect3,1);
    mesh.audit.nT3Vertices=size(coord3,1);
    mesh.audit.nT6Nodes=size(coord6,1);
    mesh.audit.maxMirrorCoordError=max(mirrorErr);
    mesh.audit.upperVertexMirrorMap=mirrorMap;
    mesh.audit.upperCrackFaceIDs=idU(:,end);
    mesh.audit.lowerCrackFaceIDs=mirrorMap(idU(:,end));
    mesh.audit.crackFacesDistinct=all(idU(:,end) ~= mirrorMap(idU(:,end)));

    if ~mesh.audit.crackFacesDistinct
        error('validate_EDI_Williams_fields:CrackFacesJoined', ...
            'Upper/lower negative-x crack-face node IDs must be distinct.');
    end
end


function mesh = build_legacy_same_diagonal_annulus(r0,r1,Nr,Nth)

    rv=linspace(r0,r1,Nr+1);
    tv=linspace(-pi,pi,Nth+1);

    nR=numel(rv);
    nT=numel(tv);

    coord3=zeros(nR*nT,2);
    id=@(ir,it) (it-1)*nR+ir;

    for it=1:nT
        th=tv(it);
        for ir=1:nR
            r=rv(ir);
            coord3(id(ir,it),:)=r*[cos(th),sin(th)];
        end
    end

    connect3=zeros(2*Nr*Nth,3);
    e=0;

    for it=1:Nth
        for ir=1:Nr
            A=id(ir,it);
            B=id(ir+1,it);
            C=id(ir+1,it+1);
            D=id(ir,it+1);

            e=e+1; connect3(e,:)=[A B C];
            e=e+1; connect3(e,:)=[A C D];
        end
    end

    area=tri_area_signed(connect3,coord3);
    cw=area<0;
    if any(cw)
        tmp=connect3(cw,2);
        connect3(cw,2)=connect3(cw,3);
        connect3(cw,3)=tmp;
    end

    [coord6,connect6]=T3toT6_fast(coord3,connect3);

    mesh=struct();
    mesh.coord3=coord3;
    mesh.connect3=connect3;
    mesh.coord=coord6;
    mesh.connect=connect6;
    mesh.audit=struct( ...
        'topology','legacy_same_diagonal', ...
        'Nr',Nr,'Nth',Nth, ...
        'nT3',size(connect3,1), ...
        'nT3Vertices',size(coord3,1), ...
        'nT6Nodes',size(coord6,1), ...
        'maxMirrorCoordError',NaN, ...
        'crackFacesDistinct',true);
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
