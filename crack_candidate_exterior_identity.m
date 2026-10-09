function identity=crack_candidate_exterior_identity(candidate)
scale=1;
if isfield(candidate,'coreMeshControls'),scale=candidate.coreMeshControls.scale;end
controls=struct();
if isfield(candidate,'exteriorMeshControls'),controls=candidate.exteriorMeshControls;end
identity=crack_exterior_identity(controls,scale);
end
