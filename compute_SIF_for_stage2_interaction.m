function [KI, KII, Aux] = compute_SIF_for_stage2_interaction(C, G2, Mc, S2, varargin)
%COMPUTE_SIF_FOR_STAGE2_INTERACTION
% Adapter for SIF_LEFM_interaction_EDI in the Stage-II crack-from-hole workflow.

    ip = inputParser;
    addParameter(ip, 'r_inner', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
    addParameter(ip, 'r_outer', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(ip, 'inner_factor', 0.25, @(x)isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
    addParameter(ip, 'outer_factor', 0.70, @(x)isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'Verbose', false, @(x)islogical(x) || isnumeric(x));
    parse(ip, varargin{:});

    meshSIF = S2.mesh;

    matSIF = S2.mat;
    if ~isfield(matSIF, 'Dmat')
        if isfield(matSIF, 'D')
            matSIF.Dmat = matSIF.D;
        else
            error('compute_SIF_for_stage2_interaction:MissingDmat', ...
                'Material struct must contain Dmat or D.');
        end
    end

    if isfield(Mc, 'crack') && isfield(Mc.crack, 'Pmid') && ~isempty(Mc.crack.Pmid)
        V = Mc.crack.Pmid;
    else
        V = G2.crack.polyline;
    end

    if isempty(ip.Results.r_outer)
        r_outer = ip.Results.outer_factor * G2.crack.a0;
    else
        r_outer = ip.Results.r_outer;
    end

    if isempty(ip.Results.r_inner)
        r_inner = ip.Results.inner_factor * r_outer;
    else
        r_inner = ip.Results.r_inner;
    end

    domain = struct();
    domain.r_inner = r_inner;
    domain.r_outer = r_outer;

    [KI, KII, Aux] = SIF_LEFM_interaction_EDI( ...
        meshSIF, S2.U, V, matSIF, domain, ...
        'AuxK', 1.0, ...
        'Verbose', logical(ip.Results.Verbose));

    Aux.meshSIF = meshSIF;
    Aux.matSIF = matSIF;
    Aux.V = V;
    Aux.domain = domain;

    Aux.r_inner = r_inner;
    Aux.r_outer = r_outer;
    Aux.r_inner_over_a0 = r_inner / G2.crack.a0;
    Aux.r_outer_over_a0 = r_outer / G2.crack.a0;
end