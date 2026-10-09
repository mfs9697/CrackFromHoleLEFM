function tf=crack_candidate_mesh_controls_match(candidate,requested)
% Candidate-cache identity, with historical ExteriorScale=CoreScale fallback.
tf=false;
try
    if isfield(candidate,'coreMeshControls')
        candidateScale=candidate.coreMeshControls.scale;
    else
        candidateScale=1;
    end
    if isfield(candidate,'exteriorMeshControls')
        exterior=candidate.exteriorMeshControls;
    else
        exterior=crack_exterior_identity(struct(),candidateScale);
    end
    tf=isnumeric(candidateScale)&&isscalar(candidateScale)&&isfinite(candidateScale)&& ...
        abs(candidateScale-requested.coreScale)<=1e-14 && ...
        crack_exterior_controls_match(exterior,requested,candidateScale,requested.coreScale);
catch
    % Corrupt metadata cannot authorize cached reuse.
end
end
