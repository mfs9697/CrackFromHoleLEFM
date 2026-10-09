function fp=stage2_native_sampling_fingerprint(candidate,mesh,mat,id)
% Independent T3 face-chain and canonical-core checks of the native T6 grid.
% Full-face counts depend on the exterior; core counts depend only on CoreScale.
if nargin<4,id='stage2phys:NativeSampling';end
a0=candidate.a0;cr=candidate.crack;
direction=diff(cr.Pmid,1,1);direction=direction/norm(direction);
local=(candidate.p-cr.Pmid(end,:))*[direction(:),[-direction(2);direction(1)]];
edges=sort([candidate.t(:,[1,2]);candidate.t(:,[2,3]);candidate.t(:,[3,1])],2);
[edges,~,group]=unique(edges,'rows');incidence=accumarray(group,1);
upper=expected_face(cr.upperNodes);lower=expected_face(cr.lowerNodes);
assert(numel(upper)==numel(lower)&&max(abs(upper-lower))<=1e-12,id, ...
    'Qualified T3 crack-face chains do not pair.');
[r,~,diag]=native_COD_audit(mesh,zeros(2*size(mesh.coord,1),1),mat,cr,8);
assert(diag.nUpper==numel(upper)&&diag.nLower==numel(lower)&& ...
    isfinite(diag.gridMismatch)&&diag.gridMismatch<=1e-12&& ...
    numel(r)==numel(upper)&&max(abs(r-upper))<=1e-12,id, ...
    'Native T6 full-face grid differs from the qualified T3 face-chain fingerprint.');
core=build_stage2_scaled_audited_core([0,0],[1,0],a0, ...
    'Scale',candidate.coreMeshControls.scale);
canonical=sort(-core.local.coord(core.crack.upperT6,1));
canonical=canonical(canonical>max(1e-12,1e-8*a0));
coreR=r(r<=candidate.coreMeshControls.rCore_m+1e-12);
assert(numel(coreR)==numel(canonical)&&max(abs(coreR-canonical))<=1e-12,id, ...
    'Native T6 core-face grid differs from the canonical structured core.');
windows=[.04,.20;.04,.30;.08,.30;.12,.30];counts=zeros(4,1);
for k=1:4,counts(k)=nnz(r/a0>=windows(k,1)&r/a0<=windows(k,2));end
assert(isequal(counts,candidate.coreMeshControls.expectedNativeSamples(:)),id, ...
    'Exact qualified native COD window counts changed.');
fp=struct('schemaVersion',1,'nUpper',numel(upper),'nLower',numel(lower), ...
    'coreFaceNodes',numel(canonical),'fullFaceR_m',upper, ...
    'coreFaceR_m',canonical,'windows',windows,'nativePoints',counts);
if isfield(candidate,'nativeSamplingFingerprint')
    saved=candidate.nativeSamplingFingerprint;
    for name={'schemaVersion','nUpper','nLower','coreFaceNodes','windows','nativePoints'}
        key=name{1};assert(isfield(saved,key)&&isequal(saved.(key),fp.(key)),id, ...
            'Saved native sampling fingerprint differs: %s.',key);
    end
    for name={'fullFaceR_m','coreFaceR_m'}
        key=name{1};assert(isfield(saved,key)&&isequal(size(saved.(key)),size(fp.(key)))&& ...
            max(abs(saved.(key)-fp.(key)))<=1e-12,id, ...
            'Saved native face abscissae differ: %s.',key);
    end
end

    function radii=expected_face(nodes)
        nodes=unique(nodes(:));[positions,order]=sort(-local(nodes,1));nodes=nodes(order);
        assert(nodes(1)==cr.tipNode&&abs(positions(1))<=1e-12&& ...
            abs(positions(end)-a0)<=1e-12&&all(diff(positions)>0)&& ...
            all(abs(local(nodes,2))<=1e-12),id, ...
            'Qualified T3 face is not a complete straight tip-to-mouth chain.');
        chain=sort([nodes(1:end-1),nodes(2:end)],2);
        [found,where]=ismember(chain,edges,'rows');
        assert(all(found)&&all(incidence(where)==1),id, ...
            'Qualified T3 face-chain edges are not distinct domain boundaries.');
        % Expected midsides derive from T3 endpoints, independently of T6 IDs/order.
        radii=sort([positions(2:end);.5*(positions(1:end-1)+positions(2:end))]);
        assert(all(diff(radii)>0),id,'Qualified face contains duplicate abscissae.');
    end
end
