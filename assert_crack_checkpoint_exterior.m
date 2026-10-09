function assert_crack_checkpoint_exterior(saved,candidate,id)
% Exact geometry/physics checks remain in the callers. This adds nominal family identity.
wanted=crack_candidate_exterior_identity(candidate);
coreScale=1;
if isfield(candidate,'coreMeshControls'),coreScale=candidate.coreMeshControls.scale;end
assert(isfield(saved,'meta'),id,'Checkpoint lacks mesh-family metadata.');
if isfield(saved.meta,'coreScale')
    assert(abs(saved.meta.coreScale-coreScale)<=1e-14,id, ...
        'Checkpoint structured core scale differs.');
end
if isfield(saved.meta,'exteriorMeshControls')
    actual=saved.meta.exteriorMeshControls;
    assert(crack_exterior_controls_match(actual,wanted,coreScale,coreScale),id, ...
        'Checkpoint exterior mesh controls differ.');
else
    % Legacy checkpoint geometry is still validated exactly by the caller.
    % Before independent scaling existed, its exterior always used CoreScale.
    assert(abs(wanted.exteriorScale-coreScale)<=1e-14,id, ...
        'Legacy checkpoint used ExteriorScale=CoreScale; isolated exterior cannot reuse it.');
end
end
