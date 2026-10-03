function Results = run_crack_path_old_vs_edi(C, sigma0, varargin)
%RUN_CRACK_PATH_OLD_VS_EDI
% Compute old mirror/J and interaction-EDI SIFs from the SAME FEM solution.
%
% Usage:
%   C = cfg_crack_path_two_leg_control(2.0);
%   R = run_crack_path_old_vs_edi(C);
%
% Name-value options:
%   'rI'              old circular contour radius
%   'r_inner'         EDI inner radius
%   'r_outer'         EDI outer radius
%   'innerFactorEDI'  used only when r_inner is omitted (default 0.1)
%   'nthet'           old contour samples (default 100)
%   'Verbose'         print summary (default true)
%
% IMPORTANT:
%   SIF_LEFM_interaction_EDI is still a prototype under independent audit.
%   Therefore old-vs-EDI differences are recorded as differences, NOT as
%   errors with EDI treated as ground truth.

    if nargin < 1 || isempty(C)
        C = cfg_crack_path_two_leg_control();
    end
    if nargin < 2 || isempty(sigma0)
        sigma0 = getf(C, 'sigma0', 1.0);
    end

    V = C.Pmid;
    lastLeg = norm(V(end,:) - V(end-1,:));

    ip = inputParser;
    addParameter(ip, 'rI', 0.5*lastLeg, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(ip, 'r_inner', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));
    addParameter(ip, 'r_outer', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
    addParameter(ip, 'innerFactorEDI', 0.1, @(x)isnumeric(x) && isscalar(x) && x>=0 && x<1);
    addParameter(ip, 'nthet', 100, @(x)isnumeric(x) && isscalar(x) && x>=10);
    addParameter(ip, 'Verbose', true, @(x)islogical(x) || isnumeric(x));
    parse(ip, varargin{:});
    S = ip.Results;

    rI = S.rI;

    if isempty(S.r_outer)
        r_outer = rI;
    else
        r_outer = S.r_outer;
    end

    if isempty(S.r_inner)
        r_inner = S.innerFactorEDI*r_outer;
    else
        r_inner = S.r_inner;
    end

    if r_inner >= r_outer
        error('run_crack_path_old_vs_edi:BadEDIDomain', ...
            'Require 0 <= r_inner < r_outer.');
    end

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);

    % ------------------------------------------------------------
    % ONE FEM solve only
    % ------------------------------------------------------------
    Sol = solve_crack_path_polyline_field(C, sigma0);

    matSIF = struct( ...
        'E', Sol.mat.E, ...
        'nu', Sol.mat.nu, ...
        'Dmat', Sol.mat.Dmat, ...
        'D', Sol.mat.D, ...
        'ps', Sol.mat.ps);

    % ------------------------------------------------------------
    % Historical mirror/J extractor
    % ------------------------------------------------------------
    [KI_old, KII_old, DbgOld] = SIF_LEFM_circle2_debug( ...
        Sol.mesh, Sol.U, Sol.V, matSIF, rI, ...
        'nthet', S.nthet, ...
        'plot', false, ...
        'verbose', false);

    % ------------------------------------------------------------
    % Interaction equivalent-domain extractor
    % ------------------------------------------------------------
    domain = struct('r_inner', r_inner, 'r_outer', r_outer);

    [KI_edi, KII_edi, DbgEDI] = SIF_LEFM_interaction_EDI( ...
        Sol.mesh, Sol.U, Sol.V, matSIF, domain, ...
        'UsePlaneStrain', Sol.mat.ps == 1, ...
        'Verbose', false);

    dKI = KI_old - KI_edi;
    dKII = KII_old - KII_edi;

    KediNorm = hypot(KI_edi, KII_edi);
    if KediNorm > 0
        vectorDifferenceRel = hypot(dKI, dKII)/KediNorm;
    else
        vectorDifferenceRel = NaN;
    end

    T = table( ...
        ["old_mirror_J"; "interaction_EDI"], ...
        [KI_old; KI_edi], ...
        [KII_old; KII_edi], ...
        [rI; r_outer], ...
        'VariableNames', {'method','KI','KII','outer_radius'});

    Results = struct();
    Results.caseConfig = C;
    Results.solution = Sol;

    Results.old = struct('KI',KI_old,'KII',KII_old,'Dbg',DbgOld,'rI',rI);
    Results.edi = struct('KI',KI_edi,'KII',KII_edi,'Dbg',DbgEDI, ...
        'r_inner',r_inner,'r_outer',r_outer);

    Results.difference = struct( ...
        'KI_old_minus_EDI', dKI, ...
        'KII_old_minus_EDI', dKII, ...
        'vector_relative_to_EDI_norm', vectorDifferenceRel);

    Results.table = T;

    Results.status = struct();
    Results.status.same_FEM_field = true;
    Results.status.old_method_role = 'historical/control';
    Results.status.edi_role = 'prototype under independent verification';
    Results.status.reference_truth_assigned = false;

    if logical(S.Verbose)
        fprintf('\n============================================================\n');
        fprintf('CRACK-PATH SIF AUDIT: SAME FEM FIELD\n');
        fprintf('============================================================\n');
        fprintf('last-leg length = %.8e\n', lastLeg);
        fprintf('old contour rI  = %.8e\n', rI);
        fprintf('EDI annulus     = [%.8e, %.8e]\n', r_inner, r_outer);
        fprintf('\n');
        disp(T);
        fprintf('old - EDI: dKI = %+ .8e, dKII = %+ .8e\n', dKI, dKII);
        fprintf('vector difference / ||K_EDI|| = %.8e\n', vectorDifferenceRel);
        fprintf(['NOTE: EDI is not yet a validated reference. ', ...
                 'This is a same-field comparison only.\n']);
    end
end


function v = getf(S, field, default)
    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        v = S.(field);
    else
        v = default;
    end
end
