function Out = sweep_same_field_radii(C, sigma0, varargin)
%SWEEP_SAME_FIELD_RADII
% Radius/domain-sensitivity gate for the SIF audit.
%
% One FEM field is solved once. The OLD mirror/J contour and the corrected
% interaction EDI are then evaluated repeatedly on that exact same field.
%
% Two sweeps are performed:
%   A) matched outer-radius sweep:
%        rI_old = r_outer_EDI = f * lastLeg
%        r_inner_EDI = innerFactor * r_outer_EDI
%   B) EDI inner-radius sweep at fixed outer radius.
%
% This separates contour/domain sensitivity from any later change in mesh.
%
% Usage:
%   O = sweep_same_field_radii(C);
%
% Name-value options:
%   'radiusFractions'   default [0.2 0.3 0.4 0.5 0.6]
%   'innerFactor'       default 0.1
%   'innerFactorsSweep' default [0.05 0.1 0.2 0.3]
%   'fixedOuterFraction'default 0.5
%   'nthet'             default 100
%   'Verbose'           default true

    if nargin < 1 || isempty(C)
        C = cfg_crack_path_two_leg_control();
    end
    if nargin < 2 || isempty(sigma0)
        if isfield(C,'sigma0') && ~isempty(C.sigma0)
            sigma0 = C.sigma0;
        else
            sigma0 = 1.0;
        end
    end

    ip = inputParser;
    addParameter(ip,'radiusFractions',[0.2 0.3 0.4 0.5 0.6], ...
        @(x)isnumeric(x) && isvector(x) && all(x>0));
    addParameter(ip,'innerFactor',0.1, ...
        @(x)isnumeric(x) && isscalar(x) && x>0 && x<1);
    addParameter(ip,'innerFactorsSweep',[0.05 0.1 0.2 0.3], ...
        @(x)isnumeric(x) && isvector(x) && all(x>0) && all(x<1));
    addParameter(ip,'fixedOuterFraction',0.5, ...
        @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip,'nthet',100, ...
        @(x)isnumeric(x) && isscalar(x) && x>=10);
    addParameter(ip,'Verbose',true, ...
        @(x)islogical(x) || isnumeric(x));
    parse(ip,varargin{:});
    S = ip.Results;

    % ------------------------------------------------------------
    % ONE FEM solve
    % ------------------------------------------------------------
    Sol = solve_crack_path_polyline_field(C,sigma0);

    matSIF = struct( ...
        'E',Sol.mat.E, ...
        'nu',Sol.mat.nu, ...
        'Dmat',Sol.mat.Dmat, ...
        'D',Sol.mat.D, ...
        'ps',Sol.mat.ps);

    V = Sol.V;
    lastLeg = norm(V(end,:) - V(end-1,:));

    % ------------------------------------------------------------
    % A. Matched outer-radius sweep
    % ------------------------------------------------------------
    rf = S.radiusFractions(:);
    nR = numel(rf);
    A = nan(nR,18);

    DbgOldA = cell(nR,1);
    DbgEDIA = cell(nR,1);

    for i = 1:nR
        r = rf(i)*lastLeg;
        rin = S.innerFactor*r;

        [KIold,KIIold,DO] = SIF_LEFM_circle2_debug( ...
            Sol.mesh,Sol.U,V,matSIF,r, ...
            'nthet',S.nthet,'plot',false,'verbose',false);

        domain = struct('r_inner',rin,'r_outer',r);
        [KIedi,KIIedi,DE] = SIF_LEFM_interaction_EDI( ...
            Sol.mesh,Sol.U,V,matSIF,domain, ...
            'UsePlaneStrain',Sol.mat.ps==1,'Verbose',false);

        dKI = KIold-KIedi;
        dKII = KIIold-KIIedi;
        kn = hypot(KIedi,KIIedi);

        detDen = 0.5*(abs(DO.detJP)+abs(DO.detJQ));
        goodDet = detDen > 0;
        detMismatch = NaN;
        if any(goodDet)
            detMismatch = median(abs(DO.detJP(goodDet)-DO.detJQ(goodDet)) ...
                ./ detDen(goodDet));
        end

        baryMismatch = median(abs(DO.baryMinP-DO.baryMinQ));

        A(i,:) = [ ...
            rf(i), r, rin, ...
            KIold, KIIold, KIedi, KIIedi, ...
            dKI, dKII, safe_div(abs(dKI),abs(KIedi)), ...
            safe_div(abs(dKII),abs(KIIedi)), safe_div(hypot(dKI,dKII),kn), ...
            DO.JI_over_absint, DO.JII_over_absint, ...
            DO.diagnostics.J1_sign_changes, DO.diagnostics.J2_sign_changes, ...
            detMismatch, baryMismatch];

        DbgOldA{i} = DO;
        DbgEDIA{i} = DE;
    end

    Touter = array2table(A,'VariableNames',{ ...
        'r_over_lastLeg','r_outer','r_inner', ...
        'KI_old','KII_old','KI_EDI','KII_EDI', ...
        'dKI_old_minus_EDI','dKII_old_minus_EDI', ...
        'abs_dKI_over_abs_KI_EDI','abs_dKII_over_abs_KII_EDI', ...
        'vector_difference_rel', ...
        'old_JI_over_absint','old_JII_over_absint', ...
        'old_J1_sign_changes','old_J2_sign_changes', ...
        'median_detJ_PQ_mismatch','median_baryMin_PQ_mismatch'});

    % ------------------------------------------------------------
    % B. EDI inner-radius sweep at fixed outer radius
    % ------------------------------------------------------------
    fixedOuter = S.fixedOuterFraction*lastLeg;
    infac = S.innerFactorsSweep(:);
    nI = numel(infac);
    B = nan(nI,7);
    DbgEDIB = cell(nI,1);

    for i = 1:nI
        rin = infac(i)*fixedOuter;
        domain = struct('r_inner',rin,'r_outer',fixedOuter);

        [KIedi,KIIedi,DE] = SIF_LEFM_interaction_EDI( ...
            Sol.mesh,Sol.U,V,matSIF,domain, ...
            'UsePlaneStrain',Sol.mat.ps==1,'Verbose',false);

        B(i,:) = [ ...
            infac(i), rin, fixedOuter, ...
            KIedi, KIIedi, DE.nElem_used, DE.nGP_used];

        DbgEDIB{i} = DE;
    end

    Tinner = array2table(B,'VariableNames',{ ...
        'inner_over_outer','r_inner','r_outer', ...
        'KI_EDI','KII_EDI','nElem_used','nGP_used'});

    Out = struct();
    Out.solution = Sol;
    Out.settings = S;
    Out.lastLeg = lastLeg;
    Out.outerSweep = Touter;
    Out.innerSweep = Tinner;
    Out.debugOldOuter = DbgOldA;
    Out.debugEDIOuter = DbgEDIA;
    Out.debugEDIInner = DbgEDIB;

    if logical(S.Verbose)
        fprintf('\n============================================================\n');
        fprintf('SIF AUDIT STEP 3A: SAME-FIELD OUTER-RADIUS SWEEP\n');
        fprintf('============================================================\n');
        fprintf('last-leg length = %.8e\n',lastLeg);
        fprintf('EDI inner/outer = %.4g in sweep A\n\n',S.innerFactor);
        disp(Touter(:,{ ...
            'r_over_lastLeg','KI_old','KI_EDI','KII_old','KII_EDI', ...
            'abs_dKI_over_abs_KI_EDI','abs_dKII_over_abs_KII_EDI', ...
            'vector_difference_rel','old_JII_over_absint', ...
            'median_detJ_PQ_mismatch','median_baryMin_PQ_mismatch'}));

        fprintf('\n============================================================\n');
        fprintf('SIF AUDIT STEP 3B: EDI INNER-RADIUS SWEEP\n');
        fprintf('============================================================\n');
        fprintf('fixed r_outer/lastLeg = %.4g\n\n',S.fixedOuterFraction);
        disp(Tinner);
    end
end


function y = safe_div(a,b)
    if isfinite(b) && abs(b)>0
        y = a/b;
    else
        y = NaN;
    end
end
