function report=test_literal_crack_cut()
%TEST_LITERAL_CRACK_CUT Geometry-only regression and destructive-input checks.
% No stiffness assembly, Williams displacement field, or SIF solve is called.
% Run from the repository root with:
%   addpath('verification/sif_audit'); report=test_literal_crack_cut();

[mesh,info,parent]=build_literal_ring_crack_cut_mesh();
cut=info.cut;
assert(info.auditT3.passed && strcmp(info.auditT3.stage,'T3'));
assert(info.auditT6.passed && strcmp(info.auditT6.stage,'T6'));

% Independently lock the approved parent geometry. A parent drift must not
% become acceptable simply because the cut preserves that changed parent.
assert(isequal(size(parent.coord3),[3008 2]));
assert(isequal(size(parent.connect3),[5888 3]));
assert(numel(parent.rings)==47 && all(cellfun(@numel,parent.rings)==64));
for j=1:47
    ids=parent.rings{j};
    assert(isequal(ids,(j-1)*64+(1:64)));
    r=0.005*40^((j-1)/46);
    theta=(0:63).'*pi/32+mod(j-1,2)*pi/64;
    expected=r*[cos(theta),sin(theta)];
    actual=parent.coord3(ids,:);
    assert(max(abs(actual(:)-expected(:)))<2e-14*r);
    chord=sqrt(sum((actual-actual([2:end 1],:)).^2,2));
    assert((max(chord)-min(chord))/mean(chord)<1e-12);
end
assert(abs(info.parent.qActual-40^(1/46))<1e-15);

assert(isequal(size(mesh.coord3),[3078 2]));
assert(isequal(size(mesh.connect3),[5934 3]));
assert(isequal(size(mesh.coord),[12089 2]));
assert(isequal(size(mesh.connect),[5934 6]));
assert(numel(cut.originalCrackNodeIDs)==24);
assert(numel(cut.insertedEdgeNodeIDs)==23);
assert(numel(cut.crackUpperIDs)==47 && numel(mesh.crackUpperT6IDs)==93);
assert(numel(cut.cutNeighborhoodParentIDs)==138);
bands=(0:45).';
assert(isequal(cut.splitParentIDs,128*bands+64-mod(bands,2)));
assert(info.auditT3.nBoundaryEdges==220);
assert(abs(info.auditT3.qualityCutMin-0.751266192754879)<1e-12);
assert(info.auditT3.minAngleCutDeg>30 && info.auditT3.minAngleCutDeg<30.1);
assert(info.auditT3.areaRelativeErrorMax<1e-12);
% Intersections belong to straight phase-half chords, not projected circles.
insertedX=sort(-mesh.coord3(cut.insertedEdgeNodeIDs,1));
expectedX=0.005*40.^((1:2:45).'/46)*cos(pi/64);
assert(max(abs(insertedX-expectedX))<1e-14);
assert(mesh.coord3(cut.crackUpperIDs(1),1)==-0.20);
assert(mesh.coord3(cut.crackUpperIDs(end),1)==-0.005);

t3=rmfield(mesh,{'coord','connect','crackUpperT6IDs','crackLowerT6IDs'});
failures={};

bad=t3; bad.coord3(1,1)=bad.coord3(1,1)+1e-6;
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'moved original off-seam node','OutsideNodes');

outside=setdiff((1:size(parent.connect3,1)).',cut.cutNeighborhoodParentIDs);
bad=t3; bad.connect3(outside(1),[2 3])=bad.connect3(outside(1),[3 2]);
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'modified element outside cut neighborhood','OutsideConnectivity');

bad=t3; bad.coord3(cut.crackLowerIDs(10),2)=1e-7;
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'off-axis crack-face node','FaceCoordinates');

badCut=cut; badCut.crackLowerIDs(1)=badCut.crackUpperIDs(1);
failures{end+1}=reject(@()validate_literal_crack_cut(parent,t3,badCut), ...
    'same upper and lower face ID','DistinctFaces');

bad=t3;
for k=1:numel(cut.crackLowerIDs)
    bad.connect3(bad.connect3==cut.crackLowerIDs(k))=cut.crackUpperIDs(k);
end
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'welded complete crack seam','UnusedNodes');

bad=t3; e=cut.splitParentIDs(1); bad.connect3(e,:)=parent.connect3(e,:);
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'element bridging the crack','Crossing');

bad=t3; e=cut.splitParentIDs(1); bad.connect3(e,[2 3])=bad.connect3(e,[3 2]);
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'negative signed triangle area','PositiveArea');

badCut=cut; badCut.parentElementID(end)=0;
failures{end+1}=reject(@()validate_literal_crack_cut(parent,t3,badCut), ...
    'invalid parent provenance','ParentMap');

badCut=cut; badCut.cutNeighborhoodParentIDs(end)=[];
failures{end+1}=reject(@()validate_literal_crack_cut(parent,t3,badCut), ...
    'underreported cut neighborhood','CutNeighborhood');

bad=t3; pair=find(cut.parentElementID==cut.splitParentIDs(1));
bad.connect3(pair(2),:)=bad.connect3(pair(1),:);
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'overlapping children and missing lower child', ...
    {'Manifold','UnusedNodes','AreaConservation','SeamBoundary','FaceTopology','BoundaryPreservation'});

failures{end+1}=reject(@()sliver_parent(parent,cut), ...
    'near-cut sliver despite valid positive areas','CutQuality');

bad=mesh; bad.coord(bad.connect(1,4),2)=bad.coord(bad.connect(1,4),2)+1e-6;
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'incorrect T6 midpoint coordinates','T6Midpoints');

bad=mesh; bad.connect(1,4)=bad.connect(1,5);
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'wrong T6 midpoint ID','T6EdgeIDs');

bad=mesh;
for k=2:2:numel(mesh.crackLowerT6IDs)
    bad.connect(bad.connect==mesh.crackLowerT6IDs(k))=mesh.crackUpperT6IDs(k);
end
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'welded T6 seam midpoints','T6EdgeIDs');

bad=mesh; bad.crackUpperT6IDs(1)=[]; bad.crackLowerT6IDs(1)=[];
failures{end+1}=reject(@()validate_literal_crack_cut(parent,bad,cut), ...
    'missing crack endpoint in complete T6 face list','T6FaceCompleteness');

report=struct('passed',true,'nRejectedMutations',numel(failures), ...
    'rejectedMutations',{failures},'auditT3',info.auditT3,'auditT6',info.auditT6);
fprintf('Literal crack geometry regression PASSED: approved parent, T3/T6 gates, %d rejected mutations.\n', ...
    report.nRejectedMutations);
end

function result=reject(action,label,expected)
if ischar(expected), expected={expected}; end
try
    action();
catch err
    identifiers=cellfun(@(s)['literalCrackCut:',s],expected,'UniformOutput',false);
    assert(any(strcmp(err.identifier,identifiers)), ...
        'literalCrackCutTest:UnexpectedFailure','%s failed unexpectedly: %s (%s)', ...
        label,err.message,err.identifier);
    result=struct('case',label,'identifier',err.identifier);
    return;
end
error('literalCrackCutTest:UndetectedMutation','Validator accepted %s.',label);
end

function sliver_parent(parent,cut)
% Deliberately perturb the INPUT ONLY for this negative test. Preserve the
% annulus topology and positive parent areas but flatten the innermost cut
% triangle so only the independent quality gate can reject its children.
e=cut.splitParentIDs(1); tri=parent.connect3(e,:);
on=find(ismember(tri,cut.originalCrackNodeIDs));
off=setdiff(1:3,on);
parent.coord3(tri(on),1)=mean(parent.coord3(tri(off),1))+1e-8;
[mesh,badCut]=cut_negative_x_ray_t3(parent);
validate_literal_crack_cut(parent,mesh,badCut);
end
