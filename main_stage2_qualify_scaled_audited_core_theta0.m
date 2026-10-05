function Q = main_stage2_qualify_scaled_audited_core_theta0(varargin)
%MAIN_STAGE2_QUALIFY_SCALED_AUDITED_CORE_THETA0
% Stage II-B1/B2 reference qualification, NO PHYSICAL FEM SOLVE.
%
% Construct the closed-audit reflection-paired crack-tip/core topology,
% rescaled from the historical a0=8 mm audit to the current frozen a0,
% place it at the frozen Stage-I initiation point with theta_1=0 deg,
% upgrade T3->T6, and qualify the topology/extractor using prescribed
% leading Williams displacement fields.
%
% Audited nondimensional transfer:
%   hTip/a0       = 0.00675308135
%   rInner/a0     = 0.10
%   rOuter/a0     = 0.65
%   rCore/a0      = 0.75
%   h(r)/a0       = s*[0.00675308135 + 0.028*(r/a0)]
%
% Cases:
%   pure I       : (KI,KII) = (1,0)
%   pure II      : (KI,KII) = (0,1)
%   tiny mixed   : (KI,KII) = (1,1e-4)
%
% Extraction:
%   SIF_LEFM_interaction_EDI
%   WeightFunction = 'fe_nodal'
%   QuadratureRule = 16
%
% The exact prescribed displacement uses topological crack-face labels:
% upper face theta=+pi, lower face theta=-pi. This is essential because
% collapsed upper/lower coordinates coincide.
%
% Usage:
%   Q = main_stage2_qualify_scaled_audited_core_theta0( ...
%       'FrozenState', R0);
%
% This driver never calls solve_cracked_LEFM.

    ip = inputParser;
    addParameter(ip,'FrozenState',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'StateFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'Scale',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<=1);
    addParameter(ip,'SaveCompact',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'OutputFile','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'VerboseEDI',false,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opt = ip.Results;

    repoRoot = fileparts(mfilename('fullpath'));
    [R0,sourceLabel,stateFile] = local_load_frozen_state( ...
        repoRoot,opt.FrozenState,char(opt.StateFile));

    assert(isfield(R0,'summary')&&istable(R0.summary)&&height(R0.summary)==1, ...
        'stage2core:FrozenSummary','Frozen R0.summary is required.');
    assert(isfield(R0,'C')&&isstruct(R0.C), ...
        'stage2core:FrozenConfig','Frozen R0.C is required.');

    T0 = R0.summary;
    req = {'a0_reserved_m','x_star_m','y_star_m', ...
        'nmat_x','nmat_y','that_x','that_y','stage1_pass'};
    local_require_table_variables(T0,req);
    row = T0(1,:);

    assert(logical(row.stage1_pass), ...
        'stage2core:Stage1NotPassed','Frozen Stage-I state did not pass.');

    C = R0.C;
    a0 = row.a0_reserved_m;
    xMouth = [row.x_star_m,row.y_star_m];
    e1 = [row.nmat_x,row.nmat_y];
    e1 = e1/norm(e1);
    e2 = [-e1(2),e1(1)];

    % Narrow experiment: theta_1 is exactly zero in the frozen local frame.
    xTip = xMouth + a0*e1;

    hTipOverA0 = 0.00675308135;
    rInnerOverA0 = 0.10;
    rOuterOverA0 = 0.65;
    rCoreOverA0 = 0.75;

    rInner = rInnerOverA0*a0;
    rOuter = rOuterOverA0*a0;
    rCore = rCoreOverA0*a0;

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('STAGE II: a0-SCALED AUDITED CORE, theta_1 = 0 deg\n');
    fprintf('============================================================\n');
    fprintf('  Frozen source : %s\n',sourceLabel);
    fprintf('  NO physical FEM solve is performed.\n');
    fprintf('  a0            = %.9f mm\n',1e3*a0);
    fprintf('  mouth         = [%.12g, %.12g] m\n',xMouth(1),xMouth(2));
    fprintf('  tip           = [%.12g, %.12g] m\n',xTip(1),xTip(2));
    fprintf('  e1            = [%.12g, %.12g]\n',e1(1),e1(2));
    fprintf('  scale s       = %.12g\n',opt.Scale);
    fprintf('  hTip/a0       = %.11g\n',hTipOverA0);
    fprintf('  hTip          = %.9f mm\n',1e3*opt.Scale*hTipOverA0*a0);
    fprintf('  EDI annulus   = [%.6f, %.6f] mm = [0.10,0.65] a0\n', ...
        1e3*rInner,1e3*rOuter);
    fprintf('  paired core   = %.6f mm = 0.75 a0\n\n',1e3*rCore);

    Core = build_stage2_scaled_audited_core( ...
        xTip,e1,a0, ...
        'Scale',opt.Scale, ...
        'HTipOverA0',hTipOverA0, ...
        'RCoreOverA0',rCoreOverA0);

    % ------------------------------------------------------------
    % Geometry / topology gates
    % ------------------------------------------------------------
    X3 = Core.local.coord3;
    T3 = Core.local.connect3;
    X6 = Core.local.coord;
    T6 = Core.local.connect;

    area2 = local_twice_signed_area(X3,T3);
    minArea2 = min(area2);

    [minDetJ,nBadJ] = local_min_t6_detj(X6,T6);

    tipNode = Core.crack.tipNode;
    tipElems = find(any(T3==tipNode,2));
    tipNeighbors = unique(T3(tipElems,:));
    tipNeighbors(tipNeighbors==tipNode) = [];

    up3 = Core.crack.upperT3(:);
    lo3 = Core.crack.lowerT3(:);
    up6 = Core.crack.upperT6(:);
    lo6 = Core.crack.lowerT6(:);

    shared3 = intersect(up3,lo3);
    shared6 = intersect(up6,lo6);

    rUp3 = sort(vecnorm(X3(up3,:),2,2));
    rLo3 = sort(vecnorm(X3(lo3,:),2,2));
    rUp6 = sort(vecnorm(X6(up6,:),2,2));
    rLo6 = sort(vecnorm(X6(lo6,:),2,2));

    mismatch3 = Inf;
    if numel(rUp3)==numel(rLo3)
        mismatch3 = max(abs(rUp3-rLo3));
    end

    mismatch6 = Inf;
    if numel(rUp6)==numel(rLo6)
        mismatch6 = max(abs(rUp6-rLo6));
    end

    nUpper3 = Core.mirror.nUpperT3Nodes;
    mirror3 = Core.mirror.T3(:);
    ref3 = [X3(1:nUpper3,1),-X3(1:nUpper3,2)];
    mirrorErr3 = max(vecnorm(X3(mirror3,:)-ref3,2,2));

    upperElems = (1:Core.mirror.nUpperT3Elements)';
    upper6 = unique(T6(upperElems,:));
    mirror6 = Core.mirror.T6;
    ref6 = [X6(upper6,1),-X6(upper6,2)];
    mirrorErr6 = max(vecnorm(X6(mirror6(upper6),:)-ref6,2,2));

    % Global transform check.
    tipGlobalErr = norm(Core.global.coord3(tipNode,:)-xTip);
    e2Frozen = [row.that_x,row.that_y];
    frameErr = norm(e2-e2Frozen/norm(e2Frozen));

    % Core must lie strictly inside material for theta_1=0. The nearest
    % circular-hole boundary is a0 behind the tip in this radial case.
    holeCenter = C.hole.center;
    holeR = C.hole.r;
    tipHoleClearance = norm(xTip-holeCenter)-holeR;
    coreClearOfHole = rCore < tipHoleClearance - 1e-12*max(1,a0);

    % EDI q-support must be strictly inside the paired core.
    ediInsideCore = rOuter < rCore;

    gates = struct();
    gates.frozenStage1Pass = logical(row.stage1_pass);
    gates.theta0FrameConsistent = frameErr <= 1e-10;
    gates.tipLocationCorrect = tipGlobalErr <= 1e-12*max(1,a0);
    gates.auditedRingCount59 = Core.design.nRings==59;
    gates.auditedT3NodeCount = size(X3,1)==6486;
    gates.auditedT3ElementCount = size(T3,1)==12678;
    gates.auditedT6NodeCount = size(X6,1)==25649;
    gates.auditedFaceNodeCounts = ...
        numel(up3)==60 && numel(lo3)==60 && ...
        numel(up6)==119 && numel(lo6)==119;
    gates.tipSixTriangles = numel(tipElems)==6;
    gates.tipSevenTopologicalEdges = numel(tipNeighbors)==7;
    gates.positiveT3Areas = minArea2>0;
    gates.positiveT6Jacobians = nBadJ==0 && minDetJ>0;
    gates.T3ReflectionPaired = mirrorErr3<=1e-13*max(1,a0);
    gates.T6ReflectionPaired = mirrorErr6<=1e-13*max(1,a0);
    gates.T3FacesShareOnlyTip = numel(shared3)==1 && shared3==tipNode;
    gates.T6FacesShareOnlyTip = numel(shared6)==1 && shared6==tipNode;
    gates.T3FaceAbscissaeMatched = mismatch3<=1e-13*max(1,a0);
    gates.T6FaceAbscissaeMatched = mismatch6<=1e-13*max(1,a0);
    gates.coreClearOfHole = coreClearOfHole;
    gates.EDIInsidePairedCore = ediInsideCore;

    topoPass = all(structfun(@logical,gates));

    fprintf('CORE TOPOLOGY\n');
    fprintf('  T3 nodes/elements = %d / %d\n',size(X3,1),size(T3,1));
    fprintf('  T6 nodes/elements = %d / %d\n',size(X6,1),size(T6,1));
    fprintf('  rings             = %d (audited base: 59)\n',Core.design.nRings);
    fprintf('  expected counts   = T3 6486 nodes / 12678 elements; T6 25649 nodes\n');
    fprintf('  tip T3 fan        = %d triangles / %d topological edges\n', ...
        numel(tipElems),numel(tipNeighbors));
    fprintf('  face T3 nodes     = %d / %d\n',numel(up3),numel(lo3));
    fprintf('  face T6 nodes     = %d / %d\n',numel(up6),numel(lo6));
    fprintf('  T3 face mismatch  = %.3e m\n',mismatch3);
    fprintf('  T6 face mismatch  = %.3e m\n',mismatch6);
    fprintf('  T3 mirror error   = %.3e m\n',mirrorErr3);
    fprintf('  T6 mirror error   = %.3e m\n',mirrorErr6);
    fprintf('  min T3 2*area     = %.6e m^2\n',minArea2);
    fprintf('  min T6 detJ       = %.6e m^2\n',minDetJ);
    fprintf('  tip-hole clearance= %.6f mm\n\n',1e3*tipHoleClearance);

    fprintf('TOPOLOGY GATES\n');
    local_print_gates(gates);

    if ~topoPass
        error('stage2core:TopologyFailed', ...
            'Scaled audited core failed topology qualification; no EDI was run.');
    end

    % ------------------------------------------------------------
    % Prescribed Williams fields on the exact qualified T6 topology
    % ------------------------------------------------------------
    mat = local_material(C);
    mu = mat.E/(2*(1+mat.nu));
    if mat.ps==1
        kappa = 3-4*mat.nu;
    else
        kappa = (3-mat.nu)/(1+mat.nu);
    end

    meshEDI = struct('coord',Core.global.coord,'connect',Core.global.connect);
    V = [xMouth;xTip];
    domain = struct('r_inner',rInner,'r_outer',rOuter);

    Kcase = [1 0;0 1;1 1e-4];
    caseName = ["pure_I";"pure_II";"tiny_mixed"];

    vals = nan(3,8);
    Ucase = cell(3,1);
    Aux = cell(3,1);

    for k = 1:3
        U = local_exact_williams_global( ...
            Core,Kcase(k,1),Kcase(k,2),mu,kappa);
        Ucase{k} = U;

        [KI,KII,A] = SIF_LEFM_interaction_EDI( ...
            meshEDI,U,V,mat,domain, ...
            'WeightFunction','fe_nodal', ...
            'QuadratureRule',16, ...
            'UsePlaneStrain',mat.ps==1, ...
            'StoreGPDiagnostics',false, ...
            'Verbose',opt.VerboseEDI);
        Aux{k} = A;

        vals(k,:) = [Kcase(k,1),Kcase(k,2),KI,KII, ...
            KI-Kcase(k,1),KII-Kcase(k,2),A.nElem_used,A.nGP_used];
    end

    Synthetic = array2table(vals,'VariableNames',{ ...
        'KI_input','KII_input','KI_recovered','KII_recovered', ...
        'KI_error','KII_error','nElem_used','nGP_used'});
    Synthetic.caseName = caseName;
    Synthetic = movevars(Synthetic,'caseName','Before','KI_input');

    pureI = Synthetic(1,:);
    pureII = Synthetic(2,:);
    mixed = Synthetic(3,:);

    M = [pureI.KI_recovered,pureII.KI_recovered; ...
         pureI.KII_recovered,pureII.KII_recovered];

    matrixError = norm(M-eye(2),'fro');
    mixedKIIrel = abs(mixed.KII_recovered-1e-4)/1e-4;
    mixedSuperposition = norm([mixed.KI_recovered;mixed.KII_recovered] - ...
        ([pureI.KI_recovered;pureI.KII_recovered] + ...
         1e-4*[pureII.KI_recovered;pureII.KII_recovered]));

    syntheticGates = struct();
    syntheticGates.pureIRecovery = abs(pureI.KI_recovered-1)<=2e-4;
    syntheticGates.pureICrossLeakage = abs(pureI.KII_recovered)<=1e-10;
    syntheticGates.pureIIRecovery = abs(pureII.KII_recovered-1)<=2e-4;
    syntheticGates.pureIICrossLeakage = abs(pureII.KI_recovered)<=1e-10;
    syntheticGates.recoveryMatrix = matrixError<=2e-4;
    syntheticGates.tinyMixedKII = mixedKIIrel<=2e-4;
    syntheticGates.superposition = mixedSuperposition<=1e-10;

    syntheticPass = all(structfun(@logical,syntheticGates));

    fprintf('\nPRESCRIBED WILLIAMS EDI QUALIFICATION\n');
    disp(Synthetic);
    fprintf('  recovery matrix ||M-I||_F = %.6e\n',matrixError);
    fprintf('  tiny-mixed KII relative error = %.6e\n',mixedKIIrel);
    fprintf('  superposition residual = %.6e\n\n',mixedSuperposition);
    fprintf('SYNTHETIC GATES\n');
    local_print_gates(syntheticGates);

    pass = topoPass && syntheticPass;

    if pass
        fprintf('\nSTAGE-II SCALED AUDITED CORE QUALIFICATION PASS.\n');
        fprintf('  T3/T6 topology and prescribed I/II/mixed EDI are qualified.\n');
        fprintf('  NO physical displacement field has been solved.\n');
    else
        fprintf('\nSTAGE-II SCALED AUDITED CORE QUALIFICATION FAIL.\n');
    end

    Summary = table( ...
        a0,opt.Scale,Core.hTip,rInner,rOuter,rCore, ...
        size(X3,1),size(T3,1),size(X6,1), ...
        Core.design.nRings,numel(tipElems),numel(tipNeighbors), ...
        mismatch3,mismatch6,mirrorErr3,mirrorErr6,minArea2,minDetJ, ...
        matrixError,mixedKIIrel,mixedSuperposition,pass, ...
        'VariableNames',{ ...
        'a0_m','scale','hTip_m','rInner_m','rOuter_m','rCore_m', ...
        'T3_nodes','T3_elements','T6_nodes','rings', ...
        'tip_triangles','tip_topological_edges', ...
        'T3_face_mismatch_m','T6_face_mismatch_m', ...
        'T3_mirror_error_m','T6_mirror_error_m', ...
        'min_T3_twice_area_m2','min_T6_detJ_m2', ...
        'recovery_matrix_error','tiny_mixed_KII_rel_error', ...
        'superposition_residual','pass'});

    Q = struct();
    Q.summary = Summary;
    Q.gates = gates;
    Q.syntheticGates = syntheticGates;
    Q.synthetic = Synthetic;
    Q.recoveryMatrix = M;
    Q.matrixError = matrixError;
    Q.mixedKIIrel = mixedKIIrel;
    Q.superpositionResidual = mixedSuperposition;
    Q.Core = Core;
    Q.domain = domain;
    Q.material = mat;
    Q.V = V;
    Q.source = sourceLabel;
    Q.pass = pass;

    if opt.SaveCompact
        outFile = char(opt.OutputFile);
        if isempty(outFile)
            outDir = fullfile(repoRoot,'verification','crack_path');
            if exist(outDir,'dir')~=7,mkdir(outDir);end
            outFile = fullfile(outDir,'stage2_scaled_audited_core_theta0.mat');
        else
            outDir = fileparts(outFile);
            if ~isempty(outDir)&&exist(outDir,'dir')~=7,mkdir(outDir);end
        end

        Qsave = Q;
        Qsave = rmfield(Qsave,'Core');
        Qsave.coreLocal = Core.local;
        Qsave.coreGlobal = Core.global;
        Qsave.crack = Core.crack;
        Qsave.mirror = Core.mirror;
        Qsave.design = Core.design;
        save(outFile,'Qsave','-v7');
        Q.outputFile = outFile;
        fprintf('  Compact MAT: %s\n',outFile);
    else
        Q.outputFile = '';
    end
end


% =========================================================================
function [R0,sourceLabel,stateFile] = local_load_frozen_state(repoRoot,Rin,stateFileIn)
    stateFile = '';

    if ~isempty(Rin)
        R0 = Rin;
        sourceLabel = '<in-memory FrozenState>';
        return;
    end

    if isempty(strtrim(stateFileIn))
        stateFile = fullfile(repoRoot,'verification','crack_path', ...
            'stage1_starting_state.mat');
    else
        stateFile = stateFileIn;
        if exist(stateFile,'file')~=2 && ~local_is_absolute_path(stateFile)
            q = fullfile(repoRoot,stateFile);
            if exist(q,'file')==2,stateFile=q;end
        end
    end

    if exist(stateFile,'file')~=2
        error('stage2core:MissingFrozenState', ...
            ['Frozen Stage-I MAT not found at:\n  %s\n', ...
             'Pass ''FrozenState'',R0 if the accepted struct is in memory.'], ...
            stateFile);
    end

    S = load(stateFile);
    assert(isfield(S,'R0')&&isstruct(S.R0), ...
        'stage2core:BadFrozenMAT','MAT file must contain R0.');
    R0 = S.R0;
    sourceLabel = stateFile;
end


function local_require_table_variables(T,names)
    missing = names(~ismember(names,T.Properties.VariableNames));
    if ~isempty(missing)
        error('stage2core:MissingFrozenFields', ...
            'Frozen summary missing variables: %s',strjoin(missing,', '));
    end
end


function mat = local_material(C)
    E = C.E;
    nu = C.nu;
    ps = C.ps;

    if ps==1
        coef = E/((1+nu)*(1-2*nu));
        D = coef*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
    else
        coef = E/(1-nu^2);
        D = coef*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
    end

    mat = struct('E',E,'nu',nu,'ps',ps,'D',D,'Dmat',D);
end


function U = local_exact_williams_global(Core,KI,KII,mu,kappa)
% Prescribe exact leading Williams displacement at T6 nodes. Upper/lower
% collapsed crack faces use topology to force theta=+pi/-pi respectively.

    X = Core.local.coord;
    n = size(X,1);
    th = atan2(X(:,2),X(:,1));
    r = hypot(X(:,1),X(:,2));

    up = setdiff(Core.crack.upperT6(:),Core.crack.tipNode);
    lo = setdiff(Core.crack.lowerT6(:),Core.crack.tipNode);
    th(up) = pi;
    th(lo) = -pi;

    c = cos(th/2);
    s = sin(th/2);
    fac = sqrt(r/(2*pi))/(2*mu);

    u1I = KI .* fac .* c .* (kappa-1+2*s.^2);
    u2I = KI .* fac .* s .* (kappa+1-2*c.^2);

    u1II = KII .* fac .* s .* (kappa+1+2*c.^2);
    u2II = -KII .* fac .* c .* (kappa-1-2*s.^2);

    uLocal = [u1I+u1II,u2I+u2II];
    uGlobal = uLocal*Core.Rgl.';

    U = zeros(2*n,1);
    U(1:2:end) = uGlobal(:,1);
    U(2:2:end) = uGlobal(:,2);
end


function A2 = local_twice_signed_area(X,T)
    a = X(T(:,2),:)-X(T(:,1),:);
    b = X(T(:,3),:)-X(T(:,1),:);
    A2 = a(:,1).*b(:,2)-a(:,2).*b(:,1);
end


function [minDet,nBad] = local_min_t6_detj(X,T)
    [xip,~] = local_rule16();
    minDet = Inf;
    nBad = 0;

    for e = 1:size(T,1)
        Xe = X(T(e,:),:);
        for g = 1:size(xip,2)
            detJ = local_t6_detj(xip(:,g),Xe);
            minDet = min(minDet,detJ);
            if ~(isfinite(detJ)&&detJ>0)
                nBad = nBad+1;
            end
        end
    end
end


function detJ = local_t6_detj(xi,X)
    L1 = xi(1);
    L2 = xi(2);
    L3 = 1-L1-L2;

    dN1 = [ ...
        4*L1-1;
        0;
        -(4*L3-1);
        4*L2;
        -4*L2;
        4*(L3-L1)];

    dN2 = [ ...
        0;
        4*L2-1;
        -(4*L3-1);
        4*L1;
        4*(L3-L2);
        -4*L1];

    J = [dN1.';dN2.']*X;
    detJ = det(J);
end


function [xip,w] = local_rule16()
    xip = zeros(2,16);
    w = zeros(1,16);

    xip(:,1) = [1/3;1/3];
    w(1) = 0.144315607677787;

    a = 0.170569307751760;
    b = 0.658861384496480;
    xip(:,2:4) = [a a b;a b a];
    w(2:4) = 0.103217370534718;

    a = 0.050547228317031;
    b = 0.898905543365938;
    xip(:,5:7) = [a a b;a b a];
    w(5:7) = 0.032458497623198;

    a = 0.459292588292723;
    b = 0.081414823414554;
    xip(:,8:10) = [a a b;a b a];
    w(8:10) = 0.095091634267285;

    a = 0.263112829634638;
    b = 0.728492392955404;
    d = 0.008394777409958;
    xip(:,11:16) = [a a b b d d;b d a d a b];
    w(11:16) = 0.027230314174435;
end


function local_print_gates(S)
    fn = fieldnames(S);
    for k = 1:numel(fn)
        fprintf('  %-32s : %d\n',fn{k},logical(S.(fn{k})));
    end
end


function tf = local_is_absolute_path(p)
    p = char(p);
    tf = startsWith(p,filesep) || ...
        ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once')) || ...
        startsWith(p,'\\');
end
