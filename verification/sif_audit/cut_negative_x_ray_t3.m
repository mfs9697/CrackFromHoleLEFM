function [mesh,cut] = cut_negative_x_ray_t3(parent)
%CUT_NEGATIVE_X_RAY_T3 Insert a traction-free negative-x seam into a T3 mesh.
% The input annulus is immutable. Original node IDs and original element rows
% are retained; extra children and nodes are appended. Only triangles whose
% interiors meet the ray are subdivided. Lower triangles touching the ray at
% a vertex merely substitute that vertex's lower-face ID.
%
% This cutter requires the ray origin to lie outside the material (as in the
% literal annulus). It does not insert a crack tip inside an element. Existing
% negative-axis vertices have only their sin(pi) roundoff canonicalized to
% exact zero in the output; every other original coordinate is unchanged.
X = parent.coord3;
T = parent.connect3;
n0 = size(X,1);
ne = size(T,1);
assert(size(X,2)==2 && size(T,2)==3 && all(isfinite(X(:))), ...
    'cut_negative_x_ray_t3:InvalidInput','Require finite 2D T3 data.');
tol = 64*eps(max(abs(X(:))));
originalCrack = find(X(:,1)<0 & abs(X(:,2))<=tol);
X(originalCrack,2) = 0;

% One intersection per original edge, shared by both incident triangles.
E = unique(sort([T(:,[1 2]); T(:,[2 3]); T(:,[3 1])],2),'rows');
y1 = X(E(:,1),2); y2 = X(E(:,2),2);
crossesAxis = (y1>tol & y2<-tol) | (y1<-tol & y2>tol);
Ec = E(crossesAxis,:);
fraction = -X(Ec(:,1),2)./(X(Ec(:,2),2)-X(Ec(:,1),2));
points = X(Ec(:,1),:) + fraction.*(X(Ec(:,2),:)-X(Ec(:,1),:));
cutEdge = points(:,1)<-tol;
sourceEdges = Ec(cutEdge,:);
points = points(cutEdge,:);
points(:,2) = 0;
inserted = (n0+1:n0+size(points,1)).';
X = [X; points];
edgeNode = sparse(sourceEdges(:,1),sourceEdges(:,2),inserted,n0,n0);

parentID = (1:ne).';
splitIDs = zeros(0,1);
touched = any(ismember(T,originalCrack),2);
for e = 1:ne
    ids = T(e,:);
    y = X(ids,2);
    if ~(any(y>tol) && any(y<-tol)), continue; end
    pairs = sort(ids([1 2;2 3;3 1]),2);
    mids = full(edgeNode(sub2ind([n0 n0],pairs(:,1),pairs(:,2))));
    nIntersections = nnz(mids) + nnz(ismember(ids,originalCrack));
    if nIntersections==0, continue; end % positive-x axis is not a crack
    assert(nIntersections==2, 'cut_negative_x_ray_t3:InteriorTip', ...
        'The ray origin must lie outside the material.');
    touched(e) = true;
    splitIDs(end+1,1) = e; %#ok<AGROW>
    upper = clip_half(ids,X,edgeNode,n0,+1);
    lower = clip_half(ids,X,edgeNode,n0,-1);
    children = [triangulate_polygon(upper,X); triangulate_polygon(lower,X)];
    T(e,:) = children(1,:);
    T = [T; children(2:end,:)]; %#ok<AGROW>
    parentID = [parentID; repmat(e,size(children,1)-1,1)]; %#ok<AGROW>
end

% Sorting establishes a one-to-one geometric face pairing from outer to inner.
upper = [originalCrack; inserted];
[~,order] = sort(X(upper,1));
upper = upper(order);
lower = (size(X,1)+1:size(X,1)+numel(upper)).';
X = [X; X(upper,:)];
lowerMap = (1:size(X,1)).';
lowerMap(upper) = lower;
isLower = any(reshape(X(T,2),size(T)) < -tol,2);
hasFace = any(ismember(T,upper),2);
T(isLower & hasFace,:) = lowerMap(T(isLower & hasFace,:));

mesh = struct('coord3',X,'connect3',T);
cut = struct('parentElementID',parentID, ...
    'splitParentIDs',splitIDs,'cutNeighborhoodParentIDs',find(touched), ...
    'originalCrackNodeIDs',originalCrack,'insertedEdgeNodeIDs',inserted, ...
    'insertedSourceEdges',sourceEdges, ...
    'crackUpperIDs',upper,'crackLowerIDs',lower,'tolerance',tol);
end

function polygon = clip_half(ids,X,edgeNode,n0,side)
% Sutherland-Hodgman clipping of a CCW triangle against the half-plane.
polygon = zeros(1,0);
for k = 1:3
    a = ids(k); b = ids(mod(k,3)+1);
    if side*X(a,2)>=0, polygon(end+1)=a; end %#ok<AGROW>
    if X(a,2)*X(b,2)<0
        edge = sort([a b]);
        node = full(edgeNode(sub2ind([n0 n0],edge(1),edge(2))));
        assert(node>0,'cut_negative_x_ray_t3:InteriorTip', ...
            'Cannot clip an element containing the ray origin.');
        polygon(end+1)=node; %#ok<AGROW>
    end
end
end

function triangles = triangulate_polygon(ids,X)
if numel(ids)==3
    triangles = ids;
elseif numel(ids)==4
    a = ids([1 2 3;1 3 4]);
    b = ids([1 2 4;2 3 4]);
    % Choose the diagonal giving the better worst child; no smoothing/flips
    % of original edges or any triangles outside this parent are allowed.
    if min(quality(a,X))>=min(quality(b,X)), triangles=a; else, triangles=b; end
else
    error('cut_negative_x_ray_t3:InvalidClip','Unexpected clipped polygon.');
end
end

function q = quality(T,X)
a = X(T(:,2),:)-X(T(:,1),:);
b = X(T(:,3),:)-X(T(:,2),:);
c = X(T(:,1),:)-X(T(:,3),:);
twiceArea = abs(a(:,1).*b(:,2)-a(:,2).*b(:,1));
q = 2*sqrt(3)*twiceArea./sum(a.^2+b.^2+c.^2,2);
end
