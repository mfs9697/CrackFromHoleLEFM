function [mesh,meta]=perturb_literal_cut_mesh_lower_angles(baseMesh,baseInfo,alpha)
%PERTURB_LITERAL_CUT_MESH_LOWER_ANGLES
% Controlled one-sided angular perturbation of the approved crack-cut mesh.
%
% Only T3 vertices strictly below x2=0 are moved. Their radius is preserved
% and their polar angle is changed by
%
%   theta_new = theta + alpha*dtheta*sin(theta),   -pi < theta < 0.
%
% Thus the perturbation is zero at theta=0 and theta=-pi, maximal near the
% lower vertical axis, and proportional to one nominal sector width dtheta.
% The upper half, both crack faces, all node IDs, and all T3 connectivity are
% frozen. T6 midside nodes are regenerated from the perturbed T3 geometry.
%
% alpha=0 reproduces the approved S0 mesh exactly.

if nargin<3
    error('perturb_literal_cut_mesh_lower_angles:Input', ...
        'baseMesh, baseInfo, and alpha are required.');
end
if ~(isnumeric(alpha)&&isscalar(alpha)&&isfinite(alpha)&&alpha>=0)
    error('perturb_literal_cut_mesh_lower_angles:Alpha', ...
        'alpha must be a finite nonnegative scalar.');
end
if ~isfield(baseInfo,'parent') || ~isfield(baseInfo.parent,'dtheta')
    error('perturb_literal_cut_mesh_lower_angles:Metadata', ...
        'baseInfo.parent.dtheta is required.');
end

X0=baseMesh.coord3;
T0=baseMesh.connect3;
tol=256*eps(max(1,max(abs(X0(:)))));
dtheta=baseInfo.parent.dtheta;

X=X0;
lower=find(X0(:,2)<-tol);
theta=atan2(X0(lower,2),X0(lower,1));
r=hypot(X0(lower,1),X0(lower,2));

delta=alpha*dtheta.*sin(theta);
thetaNew=theta+delta;

X(lower,1)=r.*cos(thetaNew);
X(lower,2)=r.*sin(thetaNew);

% Crack-face corner nodes are exactly frozen.
faceT3=unique([baseInfo.cut.crackUpperIDs(:);baseInfo.cut.crackLowerIDs(:)]);
X(faceT3,:)=X0(faceT3,:);

mesh=struct('coord3',X,'connect3',T0);
[mesh.coord,mesh.connect]=T3toT6_fast(mesh.coord3,mesh.connect3);

% Preserve the complete quadratic face definitions on the regenerated T6.
mesh.crackUpperT6IDs=quadratic_face(mesh,baseInfo.cut.crackUpperIDs);
mesh.crackLowerT6IDs=quadratic_face(mesh,baseInfo.cut.crackLowerIDs);

% Hard topology/geometry checks for the controlled family.
A=signed_area(X,T0);
if any(A<=0)
    error('perturb_literal_cut_mesh_lower_angles:InvertedElement', ...
        'alpha=%.6g creates a non-positive T3 area.',alpha);
end
if ~isequal(mesh.connect3,baseMesh.connect3)
    error('perturb_literal_cut_mesh_lower_angles:ConnectivityChanged', ...
        'T3 connectivity must remain exactly frozen.');
end
if ~isequal(mesh.connect,baseMesh.connect)
    error('perturb_literal_cut_mesh_lower_angles:T6ConnectivityChanged', ...
        'T6 connectivity must remain exactly frozen.');
end

upperFrozen=X0(:,2)>tol;
if ~isequal(X(upperFrozen,:),X0(upperFrozen,:))
    error('perturb_literal_cut_mesh_lower_angles:UpperMoved', ...
        'Upper-half T3 coordinates must remain exactly frozen.');
end

faceT6=unique([mesh.crackUpperT6IDs(:);mesh.crackLowerT6IDs(:)]);
if ~isequal(mesh.coord(faceT6,:),baseMesh.coord(faceT6,:))
    error('perturb_literal_cut_mesh_lower_angles:CrackFaceMoved', ...
        'Complete T6 crack-face coordinates must remain exactly frozen.');
end

% Radii of every moved corner are preserved to roundoff.
rNew=hypot(X(lower,1),X(lower,2));
radialRel=max(abs(rNew-r)./max(r,eps));

[Q,minAngle]=quality_metrics(X,T0,A);

% Independent dimensionless asymmetry measure. For every moved lower corner,
% compare its displacement from S0 with local polar spacing r*dtheta.
moveMag=sqrt(sum((X(lower,:)-X0(lower,:)).^2,2));
localH=r*dtheta;
aLocal=moveMag./localH;

% Pairwise mirror mismatch: match each strictly upper S0 node to its exact
% reflected lower S0 partner, then measure the perturbed mismatch.
[mirrorMedian,mirrorP95,mirrorMax,nPairs]=mirror_mismatch(X0,X,tol,dtheta);

meta=struct();
meta.alpha=alpha;
meta.nMovedLowerCorners=numel(lower);
meta.maxRadialRelativeChange=radialRel;
meta.qualityMin=min(Q);
meta.qualityP05=local_percentile(Q,5);
meta.qualityMedian=median(Q);
meta.minAngleDeg=min(minAngle);
meta.localAsymMedian=median(aLocal);
meta.localAsymP95=local_percentile(aLocal,95);
meta.localAsymMax=max(aLocal);
meta.mirrorMismatchMedian=mirrorMedian;
meta.mirrorMismatchP95=mirrorP95;
meta.mirrorMismatchMax=mirrorMax;
meta.nMirrorPairs=nPairs;
meta.dtheta=dtheta;
meta.formula='theta_new = theta + alpha*dtheta*sin(theta)';
end

function ids=quadratic_face(mesh,corners)
edges=[mesh.connect(:,[1 2]);mesh.connect(:,[2 3]);mesh.connect(:,[3 1])];
mids=[mesh.connect(:,4);mesh.connect(:,5);mesh.connect(:,6)];
onFace=all(ismember(edges,corners),2);
ids=unique([corners(:);mids(onFace)]);
[~,order]=sort(mesh.coord(ids,1));
ids=ids(order);
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

function [med,p95,mx,nPairs]=mirror_mismatch(X0,X,tol,dtheta)
up=find(X0(:,2)>tol);
lo=find(X0(:,2)<-tol);
vals=nan(numel(up),1);
keep=false(numel(up),1);

for k=1:numel(up)
    target=[X0(up(k),1),-X0(up(k),2)];
    d2=sum((X0(lo,:)-target).^2,2);
    [d2min,j]=min(d2);
    scale=max(hypot(target(1),target(2))*dtheta,eps);
    if sqrt(d2min)<=1e-10*max(1,hypot(target(1),target(2)))
        vals(k)=norm(X(lo(j),:)-target)/scale;
        keep(k)=true;
    end
end

vals=vals(keep);
nPairs=numel(vals);
if isempty(vals)
    med=NaN; p95=NaN; mx=NaN;
else
    med=median(vals);
    p95=local_percentile(vals,95);
    mx=max(vals);
end
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
