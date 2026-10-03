function [mesh,info,parent] = build_literal_ring_crack_cut_mesh()
%BUILD_LITERAL_RING_CRACK_CUT_MESH Cut the approved parent, then upgrade to T6.
% Fixed geometry gate: Ntheta=64 on all 47 rings, 46 radial intervals,
% r0=0.005, r1=0.20, and alternating 0/half-sector phase. The historical
% parent builder and experimental mesh/SIF drivers are deliberately untouched.
% No FEM solve, Williams field, or SIF extraction is performed here.
[parent,parentInfo] = build_literal_ring_lattice( ...
    'r0',0.005,'r1',0.20,'Ntheta',64);
assert(parentInfo.Nr==46 && parentInfo.Ntheta==64 && ...
    abs(parentInfo.qActual-40^(1/46))<1e-14, ...
    'build_literal_ring_crack_cut_mesh:ParentChanged', ...
    'The approved literal parent lattice has changed.');
[mesh,cut] = cut_negative_x_ray_t3(parent);

% This hard gate must precede creation of any quadratic nodes.
auditT3 = validate_literal_crack_cut(parent,mesh,cut);
[mesh.coord,mesh.connect] = T3toT6_fast(mesh.coord3,mesh.connect3);
mesh.crackUpperT6IDs = quadratic_face(mesh,cut.crackUpperIDs);
mesh.crackLowerT6IDs = quadratic_face(mesh,cut.crackLowerIDs);
auditT6 = validate_literal_crack_cut(parent,mesh,cut);
info = struct('parent',parentInfo,'cut',cut, ...
    'auditT3',auditT3,'auditT6',auditT6);
end

function ids = quadratic_face(mesh,corners)
edges = [mesh.connect(:,[1 2]);mesh.connect(:,[2 3]);mesh.connect(:,[3 1])];
mids = [mesh.connect(:,4);mesh.connect(:,5);mesh.connect(:,6)];
onFace = all(ismember(edges,corners),2);
ids = unique([corners(:);mids(onFace)]);
[~,order] = sort(mesh.coord(ids,1));
ids = ids(order);
end
