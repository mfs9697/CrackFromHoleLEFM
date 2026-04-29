function [KI, KII, Aux] = compute_SIF_for_stage2(C, G2, Mc, S2, varargin)
%COMPUTE_SIF_FOR_STAGE2
% Adapter for SIF_LEFM_circle2 in the hole-crack workflow.

    ip = inputParser;
    addParameter(ip, 'rI', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    parse(ip, varargin{:});

    rI = ip.Results.rI;

    meshSIF = S2.mesh;

    if ~isfield(meshSIF, 'coord3') || ~isfield(meshSIF, 'connect3')
        error('compute_SIF_for_stage2:MissingParentT3', ...
            'S2.mesh must contain coord3 and connect3 for SIF_LEFM_circle2.');
    end

    matSIF = S2.mat;

    if ~isfield(matSIF, 'Dmat')
        if isfield(matSIF, 'D')
            matSIF.Dmat = matSIF.D;
        else
            error('compute_SIF_for_stage2:MissingDmat', ...
                'Material struct must contain Dmat or D.');
        end
    end

    % Prefer collapsed crack midline because it is the actual zero-thickness
    % crack used by the final mesh.
    if isfield(Mc, 'crack') && isfield(Mc.crack, 'Pmid') && ~isempty(Mc.crack.Pmid)
        V = Mc.crack.Pmid;
    else
        V = G2.crack.polyline;
    end

    if isempty(rI)
        rI = G2.tip.radiusJ;
    end

    %[KI, KII] = SIF_LEFM_circle2(meshSIF, S2.U, V, matSIF, rI);
    [KI, KII, Dbg] = SIF_LEFM_circle2_debug(meshSIF, S2.U, V, matSIF, rI, ...
    'nthet', 200, ...
    'plot', true, ...
    'verbose', true);

Aux.Dbg = Dbg;

    Aux = struct();
    Aux.meshSIF = meshSIF;
    Aux.matSIF  = matSIF;
    Aux.V       = V;
    Aux.rI      = rI;
end