function audit=validate_literal_crack_cut(parent,mesh,cut)
%VALIDATE_LITERAL_CRACK_CUT Independent geometry gate for the annular ray cut.
% Checks T3 before any conversion; when coord/connect are present also checks
% every T6 edge and the complete quadratic crack faces. Throws on failure.
% Parent element rows must be retained first, with extra children appended.
% This gate deliberately checks the approved annulus, not a general crack tip.

P=parent.coord3; PT=parent.connect3;
X=mesh.coord3; T=mesh.connect3;
check_arrays(P,PT,3,'Parent'); check_arrays(X,T,3,'T3');
nP=size(P,1); nPT=size(PT,1); nX=size(X,1); nT=size(T,1);
tol=128*eps(max(1,max(abs(P(:)))));
need(nX>=nP && nT>=nPT,'Size','The cut must retain all parent nodes and element rows.');
need(all(areas(P,PT)>0),'ParentArea','Parent triangles must be positively oriented.');
parentSeam=find(P(:,1)<-tol & abs(P(:,2))<=tol);
Pc=P; Pc(parentSeam,2)=0;
[touched,split]=ray_hits(Pc,PT,tol);
outside=~touched;

% Derive preservation and intersected parents directly from parent geometry.
offSeam=true(nP,1); offSeam(parentSeam)=false;
need(isequal(X(offSeam,:),P(offSeam,:)),'OutsideNodes','An original node away from the seam changed.');
need(isequal(X(parentSeam,1),P(parentSeam,1)) && all(X(parentSeam,2)==0), ...
    'OriginalSeam','Original seam vertices must retain x and canonicalize y to zero.');
need(isequal(T(outside,:),PT(outside,:)),'OutsideConnectivity', ...
    'Original element connectivity away from the ray changed.');
fields={'parentElementID','splitParentIDs','cutNeighborhoodParentIDs', ...
    'originalCrackNodeIDs','insertedEdgeNodeIDs','insertedSourceEdges', ...
    'crackUpperIDs','crackLowerIDs'};
for k=1:numel(fields)
    need(isfield(cut,fields{k}),'Metadata',['Missing cut metadata: ',fields{k}]);
end
pid=cut.parentElementID(:);
need(numel(pid)==nT && all(pid>=1 & pid<=nPT & pid==fix(pid)), ...
    'ParentMap','Invalid parent-to-child mapping.');
need(isequal(pid(1:nPT),(1:nPT).'),'ParentMap','Original element rows must be retained first.');
need(isequal(sort(cut.splitParentIDs(:)),find(split)) && ...
     isequal(sort(cut.cutNeighborhoodParentIDs(:)),find(touched)), ...
    'CutNeighborhood','Reported cut neighborhood disagrees with independent ray intersections.');
childCount=accumarray(pid,1,[nPT 1]);
need(all(childCount(~split)==1) && all(childCount(split)>=2), ...
    'SplitParents','Only crossed parents may be split, and every crossed parent must be split.');
need(isequal(sort(cut.originalCrackNodeIDs(:)),parentSeam), ...
    'OriginalSeam','Original crack-node metadata disagrees with parent coordinates.');

% The inserted vertices must be exactly the parent edges crossed by the ray.
[parentEdges,parentEdgeCount]=edge_table(PT);
pe1=Pc(parentEdges(:,1),:); pe2=Pc(parentEdges(:,2),:);
strictCross=pe1(:,2).*pe2(:,2)<0;
alpha=zeros(size(parentEdges,1),1);
alpha(strictCross)=-pe1(strictCross,2)./(pe2(strictCross,2)-pe1(strictCross,2));
hitX=pe1(:,1)+alpha.*(pe2(:,1)-pe1(:,1));
expectedEdges=parentEdges(strictCross & hitX<-tol,:);
ins=cut.insertedEdgeNodeIDs(:); source=cut.insertedSourceEdges;
need(size(source,2)==2 && size(source,1)==numel(ins) && ...
    valid_ids(ins,nX) && all(ins>nP) && numel(unique(ins))==numel(ins), ...
    'InsertedNodes','Invalid inserted edge-node metadata.');
source=sort(source,2);
need(isequal(sortrows(source),expectedEdges), ...
    'InsertedEdges','Inserted source edges disagree with independently crossed parent edges.');
for k=1:numel(ins)
    a=Pc(source(k,1),:); b=Pc(source(k,2),:);
    t=-a(2)/(b(2)-a(2)); expected=a+t*(b-a); expected(2)=0;
    need(norm(X(ins(k),:)-expected,inf)<=tol && X(ins(k),2)==0, ...
        'IntersectionPosition','An inserted vertex is not at its parent-edge ray intersection.');
end

upper=cut.crackUpperIDs(:); lower=cut.crackLowerIDs(:);
check_faces(X,upper,lower,tol,'T3');
need(isequal(sort(upper),sort([parentSeam;ins])), ...
    'FaceCompleteness','The upper face must contain every original and inserted seam vertex.');
need(isequal(sort([ins;lower]),(nP+1:nX).'), ...
    'ExtraNodes','Only parent-edge intersections and lower-face copies may be appended.');

A=areas(X,T);
need(all(A>0),'PositiveArea','A cut triangle has zero or negative signed area.');
[~,crossed]=ray_hits(X,T,tol);
need(~any(crossed),'Crossing','A cut triangle crosses the negative-x crack ray.');
need(numel(unique(T(:)))==nX,'UnusedNodes','The final T3 mesh contains unused vertices.');
parentA=areas(Pc,PT);
childA=accumarray(pid,A,[nPT 1]);
relativeAreaError=abs(childA-parentA)./parentA;
need(all(relativeAreaError<1e-10),'AreaConservation','Children do not conserve their parent areas.');
% Area conservation alone permits an overlap balanced by a gap. Containment
% plus conforming edge incidence and the preserved boundary exclude that.
for corner=1:3
    V=X(T(:,corner),:);
    for edge=1:3
        a=Pc(PT(pid,edge),:); b=Pc(PT(pid,mod(edge,3)+1),:);
        side=(b(:,1)-a(:,1)).*(V(:,2)-a(:,2))-(b(:,2)-a(:,2)).*(V(:,1)-a(:,1));
        need(all(side>=-2e-10*parentA(pid)),'Containment','A child triangle leaves its parent.');
    end
end

[edges,edgeCount,edgeGroup,edgeOwner]=edge_table(T);
need(all(edgeCount<=2),'Manifold','An edge belongs to more than two triangles.');
boundary=find(edgeCount==1);
seamEdge=all(X(edges(:,1),2)==0 & X(edges(:,2),2)==0,2) & ...
    X(edges(:,1),1)<-tol & X(edges(:,2),1)<-tol;
need(all(edgeCount(seamEdge)==1),'WeldedSeam','A crack edge is shared by two elements.');
seamBoundary=boundary(seamEdge(boundary));
need(numel(seamBoundary)==2*(numel(upper)-1),'SeamBoundary','Wrong number of crack boundary edges.');
ownerY=mean(reshape(X(T(edgeOwner(seamBoundary),:),2),[],3),2);
need(all(abs(ownerY)>tol),'SeamSide','A crack boundary element has no unambiguous face.');
upperEdges=edges(seamBoundary(ownerY>0),:);
lowerEdges=edges(seamBoundary(ownerY<0),:);
need(isequal(sortrows(upperEdges),sortrows(sort([upper(1:end-1),upper(2:end)],2))) && ...
     isequal(sortrows(lowerEdges),sortrows(sort([lower(1:end-1),lower(2:end)],2))), ...
    'FaceTopology','Crack edges do not form the complete separate upper/lower chains.');
nodeSide=zeros(nX,1); nodeSide(upper)=1; nodeSide(lower)=-1;
for corner=1:3
    s=nodeSide(T(:,corner));
    need(all(s.*mean(reshape(X(T,2),nT,3),2)>=-tol), ...
        'FaceOwnership','An element uses a node belonging to the opposite crack face.');
end
alias=(1:nX).'; alias(lower)=upper;
remainingBoundary=sort(alias(edges(boundary(~seamEdge(boundary)),:)),2);
need(isequal(sortrows(remainingBoundary),parentEdges(parentEdgeCount==1,:)), ...
    'BoundaryPreservation','The original circular polygon boundaries changed, or a hanging edge exists.');
degree=accumarray(reshape(edges(boundary,:),[],1),1,[nX 1]);
need(all(degree(degree>0)==2),'BoundaryManifold','Boundary vertices must have degree two.');
need(nX-size(edges,1)+nT==1,'Topology','The cut annulus must have Euler characteristic one.');
% All cells must belong to one component even though the crack faces split.
adj=sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],true,nX,nX);
visited=false(nX,1); frontier=1; visited(1)=true;
while ~isempty(frontier)
    next=find(any(adj(:,frontier),2) & ~visited);
    visited(next)=true; frontier=next;
end
need(all(visited),'Connectivity','The cut mesh has disconnected components.');

[Q,minAngle]=quality(X,T,A);
near=touched(pid);
% The approved cut has Qmin approximately 0.7513. This modest margin rejects
% slivers while accommodating arithmetic variation across MATLAB versions.
need(min(Q(near))>=0.70 && min(minAngle(near))>=25,'CutQuality', ...
    'Near-cut quality must satisfy Q >= 0.70 and minimum angle >= 25 degrees.');
audit=struct('passed',true,'stage','T3','nParentNodes',nP,'nParentElements',nPT, ...
    'nT3Nodes',nX,'nT3Elements',nT,'nSplitParents',nnz(split), ...
    'nTouchedParents',nnz(touched),'nCrackPairs',numel(upper), ...
    'nBoundaryEdges',numel(boundary),'areaMin',min(A),'qualityMin',min(Q), ...
    'qualityCutMin',min(Q(near)),'minAngleDeg',min(minAngle), ...
    'minAngleCutDeg',min(minAngle(near)),'areaRelativeErrorMax',max(relativeAreaError), ...
    'qualityFloor',0.70,'minAngleFloorDeg',25);
audit.checks=struct('unchangedOutsideCut',true,'faceCoordinates',true, ...
    'distinctFaceIDs',true,'positiveAreas',true,'cutQuality',true, ...
    'noCrossing',true,'areaConservation',true,'childContainment',true, ...
    'conformingManifold',true,'preservedBoundary',true,'connectedCutAnnulus',true);

if isfield(mesh,'coord') || isfield(mesh,'connect')
    need(isfield(mesh,'coord') && isfield(mesh,'connect'),'T6Fields','Both T6 coord and connect are required.');
    X6=mesh.coord; T6=mesh.connect; check_arrays(X6,T6,6,'T6');
    need(size(T6,1)==nT && isequal(T6(:,1:3),T) && isequal(X6(1:nX,:),X), ...
        'T6Corners','T6 conversion changed the validated T3 vertices or connectivity.');
    mids=[T6(:,4);T6(:,5);T6(:,6)];
    minMid=accumarray(edgeGroup,mids,[],@min); maxMid=accumarray(edgeGroup,mids,[],@max);
    need(isequal(minMid,maxMid) && all(minMid>nX) && ...
        numel(unique(minMid))==size(edges,1) && size(X6,1)==nX+size(edges,1), ...
        'T6EdgeIDs','Each topological T3 edge must have exactly one distinct new midpoint ID.');
    expectedMid=0.5*(X(edges(:,1),:)+X(edges(:,2),:));
    need(isequal(X6(minMid,:),expectedMid),'T6Midpoints','A T6 midside node is not the exact edge midpoint.');
    need(isfield(mesh,'crackUpperT6IDs') && isfield(mesh,'crackLowerT6IDs'), ...
        'T6Faces','Complete T6 face-node lists are required.');
    upper6=mesh.crackUpperT6IDs(:); lower6=mesh.crackLowerT6IDs(:);
    check_faces(X6,upper6,lower6,tol,'T6');
    expectedUpper=[upper;minMid(seamBoundary(ownerY>0))];
    expectedLower=[lower;minMid(seamBoundary(ownerY<0))];
    need(isequal(sort(upper6),sort(expectedUpper)) && isequal(sort(lower6),sort(expectedLower)), ...
        'T6FaceCompleteness','T6 faces must include every corner and boundary midpoint, including endpoints.');
    audit.stage='T6'; audit.nT6Nodes=size(X6,1); audit.nT6CrackPairs=numel(upper6);
    audit.checks.t6EdgeMidpoints=true; audit.checks.t6SeparateFaces=true;
end
end

function check_arrays(X,T,ncol,label)
need(isnumeric(X) && size(X,2)==2 && all(isfinite(X(:))) && ...
    isnumeric(T) && size(T,2)==ncol && ~isempty(T) && valid_ids(T(:),size(X,1)), ...
    'Arrays',[label,' coordinate/connectivity arrays are invalid.']);
end

function yes=valid_ids(ids,n)
yes=all(isfinite(ids) & ids==fix(ids) & ids>=1 & ids<=n);
end

function check_faces(X,upper,lower,tol,label)
need(numel(upper)==numel(lower) && numel(upper)>=2 && ...
    valid_ids(upper,size(X,1)) && valid_ids(lower,size(X,1)) && ...
    numel(unique([upper;lower]))==2*numel(upper), ...
    'DistinctFaces',[label,' upper/lower face IDs must be distinct and disjoint.']);
need(all(X([upper;lower],2)==0) && all(X([upper;lower],1)<-tol) && ...
    isequal(X(upper,:),X(lower,:)) && all(diff(X(upper,1))>tol), ...
    'FaceCoordinates',[label,' matched faces must coincide on y=0 and be ordered by increasing x.']);
end

function A=areas(X,T)
a=X(T(:,1),:); b=X(T(:,2),:); c=X(T(:,3),:);
A=0.5*((b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(b(:,2)-a(:,2)).*(c(:,1)-a(:,1)));
end

function [touched,crossed]=ray_hits(X,T,tol)
y=reshape(X(T,2),size(T));
x=reshape(X(T,1),size(T));
touched=any(abs(y)<=tol & x<-tol,2);
for j=1:3
    k=mod(j,3)+1;
    opposite=(y(:,j)>tol & y(:,k)<-tol) | (y(:,j)<-tol & y(:,k)>tol);
    ids=find(opposite);
    hit=x(ids,j)-y(ids,j).*(x(ids,k)-x(ids,j))./(y(ids,k)-y(ids,j));
    touched(ids(hit<-tol))=true;
end
crossed=touched & any(y>tol,2) & any(y<-tol,2);
end

function [edges,count,group,owner]=edge_table(T)
allEdges=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
[edges,~,group]=unique(allEdges,'rows');
count=accumarray(group,1);
owner=accumarray(group,repmat((1:size(T,1)).',3,1),[],@min);
end

function [Q,minAngle]=quality(X,T,A)
l2=zeros(size(T,1),3);
for j=1:3
    d=X(T(:,j),:)-X(T(:,mod(j,3)+1),:); l2(:,j)=sum(d.^2,2);
end
Q=4*sqrt(3)*A./sum(l2,2);
angles=zeros(size(l2));
for j=1:3
    others=setdiff(1:3,j);
    c=(sum(l2(:,others),2)-l2(:,j))./(2*sqrt(prod(l2(:,others),2)));
    angles(:,j)=acosd(max(-1,min(1,c)));
end
minAngle=min(angles,[],2);
end

function need(condition,id,message)
if ~condition, error(['literalCrackCut:',id],'%s',message); end
end
