function [mesh,meta]=flip_literal_cut_mesh_lower_diagonals(baseMesh,baseInfo,fraction)
%FLIP_LITERAL_CUT_MESH_LOWER_DIAGONALS
% Connectivity-only asymmetry on the approved literal crack-cut mesh.
%
% T3 vertex coordinates, node IDs, crack faces, ring radii, and element count
% are unchanged. In the strict lower half only, selected natural annular-cell
% diagonals are replaced by the opposite diagonal of the same convex
% quadrilateral. The upper half and crack-neighborhood topology are untouched.
%
% fraction is the fraction (0..1) of eligible lower-half annular diagonals to
% flip. Selection is deterministic and nested through a low-discrepancy score.
%
% This is the clean A_conn experiment: geometry is mirror symmetric, but
% lower-half interpolation topology is not.

if nargin<3
    error('flip_literal_cut_mesh_lower_diagonals:Input', ...
        'baseMesh, baseInfo, and fraction are required.');
end
if ~(isnumeric(fraction)&&isscalar(fraction)&&isfinite(fraction)&& ...
        fraction>=0&&fraction<=1)
    error('flip_literal_cut_mesh_lower_diagonals:Fraction', ...
        'fraction must lie in [0,1].');
end

X=baseMesh.coord3;
T0=baseMesh.connect3;
T=T0;
tol=256*eps(max(1,max(abs(X(:)))));

% Build edge -> two owner triangles.
allE=sort([T0(:,[1 2]);T0(:,[2 3]);T0(:,[3 1])],2);
owner=repmat((1:size(T0,1)).',3,1);
[edges,~,grp]=unique(allE,'rows');
count=accumarray(grp,1);
first=accumarray(grp,owner,[],@min);
last=accumarray(grp,owner,[],@max);

cand=find(count==2);
eligible=zeros(0,7); % [edgeID t1 t2 a b c d]

for kk=1:numel(cand)
    ie=cand(kk);
    t1=first(ie); t2=last(ie);
    if t1==t2, continue; end

    verts=unique([T0(t1,:),T0(t2,:)]);
    if numel(verts)~=4, continue; end

    % Strict lower half: no seam or boundary-cut element is touched.
    if any(X(verts,2)>=-tol), continue; end

    a=edges(ie,1); b=edges(ie,2);
    c=setdiff(T0(t1,:),[a b]);
    d=setdiff(T0(t2,:),[a b]);
    if numel(c)~=1||numel(d)~=1||c==d, continue; end

    % Natural annular-cell diagonal: the four corners occupy exactly two
    % ring radii, two nodes on each ring, and both old/new diagonals connect
    % across the two radii.
    ids=[a b c d];
    rr=hypot(X(ids,1),X(ids,2));
    [rlo,rhi,ok]=two_radius_levels(rr);
    if ~ok, continue; end
    level=abs(rr-rlo) < abs(rr-rhi);
    if sum(level)~=2, continue; end

    ia=find(ids==a,1); ib=find(ids==b,1);
    ic=find(ids==c,1); id=find(ids==d,1);
    if level(ia)==level(ib), continue; end
    if level(ic)==level(id), continue; end

    % Trial flip must create two positive nondegenerate triangles.
    trial=[c d a; d c b];
    trial=orient_ccw(trial,X);
    A=signed_area(X,trial);
    if any(A<=100*eps(max(A))), continue; end

    % The opposite diagonal midpoint must lie inside the quadrilateral union.
    % Positive trial areas plus area conservation is a robust convexity gate.
    Aold=sum(signed_area(X,T0([t1 t2],:)));
    if abs(sum(A)-Aold)>1e-10*Aold, continue; end

    eligible(end+1,:)=[ie,t1,t2,a,b,c,d]; %#ok<AGROW>
end

if isempty(eligible)
    error('flip_literal_cut_mesh_lower_diagonals:NoEligibleEdges', ...
        'No lower-half annular-cell diagonals were found.');
end

% Stable radial/angle ordering, then a deterministic low-discrepancy score.
mid=0.5*(X(eligible(:,4),:)+X(eligible(:,5),:));
rm=hypot(mid(:,1),mid(:,2));
th=mod(atan2(mid(:,2),mid(:,1))+2*pi,2*pi);
[~,ord]=sortrows([rm,th],[1 2]);
eligible=eligible(ord,:);

nEligible=size(eligible,1);
phi=(sqrt(5)-1)/2;
score=mod((1:nEligible)'*phi,1);

if fraction==0
    take=false(nEligible,1);
elseif fraction==1
    take=true(nEligible,1);
else
    nTake=round(fraction*nEligible);
    [~,is]=sort(score);
    take=false(nEligible,1);
    take(is(1:nTake))=true;
end

chosen=eligible(take,:);

for k=1:size(chosen,1)
    t1=chosen(k,2); t2=chosen(k,3);
    a=chosen(k,4); b=chosen(k,5);
    c=chosen(k,6); d=chosen(k,7);
    trial=orient_ccw([c d a; d c b],X);
    T(t1,:)=trial(1,:);
    T(t2,:)=trial(2,:);
end

A=signed_area(X,T);
if any(A<=0)
    error('flip_literal_cut_mesh_lower_diagonals:InvertedElement', ...
        'A diagonal flip created a non-positive triangle.');
end

mesh=struct('coord3',X,'connect3',T);
[mesh.coord,mesh.connect]=T3toT6_fast(mesh.coord3,mesh.connect3);
mesh.crackUpperT6IDs=quadratic_face(mesh,baseInfo.cut.crackUpperIDs);
mesh.crackLowerT6IDs=quadratic_face(mesh,baseInfo.cut.crackLowerIDs);

% Hard isolation checks.
if ~isequal(mesh.coord3,baseMesh.coord3)
    error('flip_literal_cut_mesh_lower_diagonals:CoordinatesChanged', ...
        'Connectivity-only experiment must preserve every T3 coordinate.');
end
if size(mesh.connect3,1)~=size(baseMesh.connect3,1)
    error('flip_literal_cut_mesh_lower_diagonals:ElementCountChanged', ...
        'Connectivity-only experiment must preserve T3 element count.');
end
if ~isequal(mesh.connect3(~ismember((1:size(T,1))',reshape(chosen(:,2:3),[],1)),:), ...
            baseMesh.connect3(~ismember((1:size(T,1))',reshape(chosen(:,2:3),[],1)),:))
    error('flip_literal_cut_mesh_lower_diagonals:OutsideConnectivityChanged', ...
        'Only owner triangles of selected diagonals may change.');
end

% Upper-half T3 rows must be exactly unchanged.
cy0=mean(reshape(X(T0,2),size(T0)),2);
upperRows=cy0>tol;
if ~isequal(T(upperRows,:),T0(upperRows,:))
    error('flip_literal_cut_mesh_lower_diagonals:UpperChanged', ...
        'Upper-half connectivity must remain exactly frozen.');
end

% Crack-face corner chains and their geometry stay unchanged.
up3=baseInfo.cut.crackUpperIDs(:);
lo3=baseInfo.cut.crackLowerIDs(:);
if ~isequal(mesh.coord3([up3;lo3],:),baseMesh.coord3([up3;lo3],:))
    error('flip_literal_cut_mesh_lower_diagonals:CrackFaceMoved', ...
        'Crack-face corner coordinates changed.');
end

[Q,minAngle]=quality_metrics(X,T,A);

meta=struct();
meta.requestedFraction=fraction;
meta.nEligible=nEligible;
meta.nFlipped=size(chosen,1);
meta.actualFraction=meta.nFlipped/nEligible;
meta.eligible=eligible;
meta.chosen=chosen;
meta.qualityMin=min(Q);
meta.qualityP05=local_percentile(Q,5);
meta.qualityMedian=median(Q);
meta.minAngleDeg=min(minAngle);
meta.nChangedTriangles=numel(unique(reshape(chosen(:,2:3),[],1)));
meta.nT3Vertices=size(X,1);
meta.nT3Elements=size(T,1);
meta.geometryChanged=false;
meta.description='same T3 vertices; selected lower annular diagonals flipped';
end

function [lo,hi,ok]=two_radius_levels(r)
rs=sort(r(:));
scale=max(rs);
tolR=1e-10*max(1,scale);
lo=mean(rs(1:2)); hi=mean(rs(3:4));
ok=abs(rs(1)-rs(2))<=tolR && abs(rs(3)-rs(4))<=tolR && ...
   hi-lo>tolR;
end

function T=orient_ccw(T,X)
A=signed_area(X,T);
cw=A<0;
if any(cw)
    tmp=T(cw,2); T(cw,2)=T(cw,3); T(cw,3)=tmp;
end
end

function A=signed_area(X,T)
a=X(T(:,1),:); b=X(T(:,2),:); c=X(T(:,3),:);
A=0.5*((b(:,1)-a(:,1)).*(c(:,2)-a(:,2)) - ...
       (b(:,2)-a(:,2)).*(c(:,1)-a(:,1)));
end

function [Q,minAngle]=quality_metrics(X,T,A)
l2=zeros(size(T,1),3);
for j=1:3
    d=X(T(:,j),:)-X(T(:,mod(j,3)+1),:);
    l2(:,j)=sum(d.^2,2);
end
Q=4*sqrt(3)*A./sum(l2,2);
angles=zeros(size(T,1),3);
for j=1:3
    others=setdiff(1:3,j);
    c=(sum(l2(:,others),2)-l2(:,j))./(2*sqrt(prod(l2(:,others),2)));
    angles(:,j)=acosd(max(-1,min(1,c)));
end
minAngle=min(angles,[],2);
end

function ids=quadratic_face(mesh,corners)
edges=[mesh.connect(:,[1 2]);mesh.connect(:,[2 3]);mesh.connect(:,[3 1])];
mids=[mesh.connect(:,4);mesh.connect(:,5);mesh.connect(:,6)];
onFace=all(ismember(edges,corners),2);
ids=unique([corners(:);mids(onFace)]);
[~,order]=sort(mesh.coord(ids,1));
ids=ids(order);
end

function y=local_percentile(x,p)
x=sort(x(:));
if isempty(x), y=NaN; return; end
z=1+(numel(x)-1)*p/100;
i=floor(z); j=ceil(z);
if i==j
    y=x(i);
else
    y=x(i)+(z-i)*(x(j)-x(i));
end
end
