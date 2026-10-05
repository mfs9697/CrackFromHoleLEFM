function [Z,T,info,physicalIDs,physicalEdges]=build_stage3c_polyline_exterior(X,Told,cr,Zp,Tp,rp,design)
%BUILD_STAGE3C_POLYLINE_EXTERIOR Deterministic full-domain exterior for a
% retained kinked crack polyline in the current-tip frame.
% Stage III-C generalizes the qualified C03 exterior construction from a
% straight negative-axis crack to an explicit retained polyline. The last
% crack leg is aligned with local +x toward the current tip; the paired core
% intersects that leg at [-rCore,0]. Every earlier path vertex is retained
% exactly as a crack-face vertex.
% Adapted from the closed Step62/C03 SIF-audit exterior construction.
% All crack-local dimensional controls are supplied through `design` so the
% historical 8-mm audit can be transferred to the current trial length a0.
% Original physical boundary geometry is fixed. All interior old triangles
% are replaced. Ring seeds follow a monotone C1 exterior size law; fixed
% geometric boundary subdivisions receive a compatible boundary layer.
% Exterior constraint subdivision is allowed; original vertices never move.
E=sort([Told(:,[1 2]);Told(:,[2 3]);Told(:,[3 1])],2);
[E,~,g]=unique(E,'rows');E=E(accumarray(g,1)==1,:);
face=all(ismember(E,[cr.upperNodes(:);cr.tipNode]),2)| ...
    all(ismember(E,[cr.lowerNodes(:);cr.tipNode]),2);
physicalEdges=E(~face,:);physicalIDs=unique(physicalEdges(:));
path=cr.Pmid;
seg=diff(path,1,1);
segLen=vecnorm(seg,2,2);
assert(size(path,1)>=3&&all(segLen>1e-12),'stage3c:BadPath', ...
    'A nondegenerate polyline with at least two legs is required.');
pathLength=sum(segLen);
% The last leg must enter the paired core along the negative local x-axis.
lastDir=seg(end,:)/segLen(end);
assert(norm(lastDir-[1 0])<=1e-10,'stage3c:LastLegFrame', ...
    'Current-tip frame must align the last crack leg with +x.');
assert(segLen(end)>rp+1e-12,'stage3c:KinkInsideCore', ...
    'Previous kink must lie outside the paired core.');
cutPath=[path(1:end-1,:);[-rp,0]];
oldIDs=physicalIDs;
xy=X(oldIDs,:);[~,ix,group]=unique(round(xy/1e-12),'rows','stable');
Q=xy(ix,:);oldMap=zeros(size(X,1),1);oldMap(oldIDs)=group;
% Geometric inner polygon; face duplication is applied after triangulation.
PE=sort([Tp(:,[1 2]);Tp(:,[2 3]);Tp(:,[3 1])],2);
[PE,~,g]=unique(PE,'rows');PE=PE(accumarray(g,1)==1,:);
rr=reshape(vecnorm(Zp(PE(:),:),2,2),size(PE));PE=PE(all(rr>rp-1e-12,2),:);
inner=Zp(unique(PE(:)),:);[~,ix]=unique(round(inner/1e-12),'rows','stable');inner=inner(ix,:);
[~,ix]=sort(atan2(inner(:,2),inner(:,1)));inner=inner(ix,:);
nOld=size(Q,1);nI=size(inner,1);Q=[Q;inner];
C=[oldMap(physicalEdges);nOld+(1:nI)',nOld+[2:nI,1]'];
% Retain and resample the complete crack polyline outside the paired core.
% Every original path vertex is inserted exactly. The final retained point
% is the rear core intersection [-rp,0].
mouthTarget=path(1,:);
dM=vecnorm(Q-mouthTarget,2,2);
[dmin,mouth]=min(dM);
assert(dmin<=1e-10,'stage3c:MouthNotFound', ...
    'Could not recover the crack mouth from the source physical boundary.');

cutInner=nOld+find(vecnorm(inner-[-rp,0],2,2)<=1e-10,1);
assert(~isempty(cutInner),'stage3c:CoreIntersection', ...
    'Could not recover the rear core-axis vertex.');

hCut=design.scale*(design.hBase_m+design.slope*rp);
cutPts=[];
for kk=1:size(cutPath,1)-1
    Aseg=cutPath(kk,:);Bseg=cutPath(kk+1,:);
    ns=max(1,ceil(norm(Bseg-Aseg)/hCut));
    Pseg=Aseg+(0:ns)'/ns.*(Bseg-Aseg);
    if kk>1,Pseg=Pseg(2:end,:);end
    cutPts=[cutPts;Pseg]; %#ok<AGROW>
end

% Use existing mouth/core vertices and add only interior retained-path points.
midPts=cutPts(2:end-1,:);
cut=[mouth;(size(Q,1)+(1:size(midPts,1)))';cutInner];
Q=[Q;midPts];
C=[C;cut(1:end-1),cut(2:end)];
C=unique(sort(C,2),'rows');
assert(all(C(:,1)~=C(:,2)),'stage3c:ZeroConstraint','Zero constraint.');
sourceTri=triangulation(Told,X);
maxR=max(vecnorm(X(physicalIDs,:),2,2));
scale=design.scale;hRp=scale*(design.hBase_m+design.slope*rp);
cal=exterior_calibration(design);
% Ring widths start consistently with the structured patch, then increase.
radial=rp;seedRows=[];ringRows=[];
while radial<maxR
    hh=exterior_h(radial,rp,hRp,design.farCap_m,scale,design.slope,cal);
    dr=design.radialFactor*hh;
    radial=radial+dr;
    n=max(6,6*ceil(2*pi*radial/hh/6));
    th=(0:n-1)'*2*pi/n;xy=radial*[cos(th),sin(th)];
    keep=~isnan(pointLocation(sourceTri,xy));
    keep=keep&~inpolygon(xy(:,1),xy(:,2),inner(:,1),inner(:,2));
    xy=xy(keep,:);
    safe=false(size(xy,1),1);
    for j=1:size(xy,1)
        safe(j)=segment_distance(xy(j,:),X,physicalEdges)>.40*hh&& ...
            point_polyline_distance(xy(j,:),cutPath)>.40*hh;
    end
    seedRows=[seedRows;xy(safe,:)]; %#ok<AGROW>
    ringRows(end+1,:)=[radial,hh,dr,n,nnz(safe)]; %#ok<AGROW>
end
% Respect the immutable physical segment spacing near boundaries. The
% interior normal comes from the saved adjacent triangle, not a nominal hole.
layer=[];
for j=1:size(physicalEdges,1)
    edge=physicalEdges(j,:);a=X(edge(1),:);b=X(edge(2),:);mid=.5*(a+b);
    hit=find(sum(ismember(Told,edge),2)==2,1);
    third=setdiff(Told(hit,:),edge);toward=X(third,:)-mid;
    v=b-a;normal=[-v(2),v(1)]/norm(v);
    if dot(normal,toward)<0,normal=-normal;end
    hh=exterior_h(norm(mid),rp,hRp,design.farCap_m,scale,design.slope,cal);
    z=mid+design.radialFactor*min(scale*norm(v),hh)*normal;
    if ~isnan(pointLocation(sourceTri,z))&&norm(z)>rp+1e-12&& ...
            point_polyline_distance(z,cutPath)>.1*min(norm(v),hh)
        layer(end+1,:)=z; %#ok<AGROW>
    end
end
% Keep ring seeds apart from immutable-geometry boundary-layer seeds.
keep=true(size(seedRows,1),1);
for j=1:size(seedRows,1)
    hh=exterior_h(norm(seedRows(j,:)),rp,hRp,design.farCap_m,scale,design.slope,cal);
    keep(j)=isempty(layer)||min(vecnorm(layer-seedRows(j,:),2,2))>.45*hh;
end
Q=[Q;layer;seedRows(keep,:)];
% Coincident seed removal never changes the first fixed boundary vertices.
[~,ix,group]=unique(round(Q/1e-12),'rows','stable');Q=Q(ix,:);C=group(C);
C=unique(sort(C,2),'rows');
fixed=unique(C(:));
[Q,T,excluded]=exterior_cells(Q,C,sourceTri,inner);
% Deterministic bounded smoothing, with fixed geometry and C1 size-law
% seeds preserved as a reference. Restore triangulation after each move.
for it=1:cal.smoothingSteps
    ED=unique(sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2),'rows');
    A=sparse([ED(:,1);ED(:,2)],[ED(:,2);ED(:,1)],1,size(Q,1),size(Q,1));
    avg=(A*Q)./max(sum(A,2),1);trial=Q;move=setdiff(unique(T(:)),fixed);
    trial(move,:)=.75*Q(move,:)+.25*avg(move,:);
    a=trial(T(:,2),:)-trial(T(:,1),:);b=trial(T(:,3),:)-trial(T(:,1),:);
    if all(a(:,1).*b(:,2)-a(:,2).*b(:,1)>0)
        [Q,T,~]=exterior_cells(trial,C,sourceTri,inner);
    end
end
% Constrained Delaunay refinement. Encroached exterior constraints may be
% subdivided, retaining the exact saved polygon (and every original vertex).
% Structured inner edges are protected. Candidate centers are spaced before
% batch insertion, preventing coincident-center refinement cascades.
refinementRows=[];nSplits=0;nPhysicalSplits=0;
for it=1:cal.refinementMaxPasses
    [minAngle,longest,L]=triangle_quality(Q,T);
    cent=(Q(T(:,1),:)+Q(T(:,2),:)+Q(T(:,3),:))/3;
    hr=exterior_h(vecnorm(cent,2,2),rp,hRp,design.farCap_m,scale,design.slope,cal);
    heff=min(hr,scale*boundary_metric(cent,X,physicalEdges,cal.boundaryMetricGrowth));
    bad=minAngle<cal.refinementMinAngle_deg-1e-8 | ...
        longest>cal.refinementLongestFactor*heff;
    [ratio,badNeighbor]=neighbor_ratio(T,longest,cal.neighborRatioTarget);bad=bad|badNeighbor;
    refinementRows(end+1,:)=[it,size(T,1),nnz(bad),min(minAngle),ratio]; %#ok<AGROW>
    if cal.verbose
        fprintf('Exterior refinement %d: %d cells, %d flagged, angle %.6g, ratio %.6g\n', ...
            it,size(T,1),nnz(bad),min(minAngle),ratio);
    end
    if ~any(bad),break,end
    ca=Q(C(:,1),:);cb=Q(C(:,2),:);cm=.5*(ca+cb);
    cR2=.25*sum((ca-cb).^2,2);
    protected=vecnorm(ca,2,2)>rp-1e-12&vecnorm(cb,2,2)>rp-1e-12 & ...
        vecnorm(ca,2,2)<rp+1e-12&vecnorm(cb,2,2)<rp+1e-12;
    split=false(size(C,1),1);add=[];addRadius=[];
    rows=find(bad);[~,order]=sortrows([minAngle(rows),-longest(rows),rows],[1 2 3]);
    edgeCols=[2 3;1 3;1 2];
    for j=rows(order).'
        aa=Q(T(j,2),:)-Q(T(j,1),:);bb=Q(T(j,3),:)-Q(T(j,1),:);
        den=2*(aa(1)*bb(2)-aa(2)*bb(1));
        cc=Q(T(j,1),:)+ ...
            [sum(aa.^2)*bb(2)-sum(bb.^2)*aa(2), ...
             aa(1)*sum(bb.^2)-bb(1)*sum(aa.^2)]/den;
        valid=~isnan(pointLocation(sourceTri,cc))&& ...
            ~inpolygon(cc(1),cc(2),inner(:,1),inner(:,2));
        if ~valid
            [~,k]=max(L(j,:));edge=sort(T(j,edgeCols(k,:)));
            [found,where]=ismember(edge,C,'rows');
            if found
                if ~protected(where),split(where)=true;end
                continue
            end
            cc=.5*(Q(edge(1),:)+Q(edge(2),:));
        end
        encroached=sum((cm-cc).^2,2)<(1-1e-10)*cR2;
        if any(encroached)
            split=split|(encroached&~protected);continue
        end
        radius=min(vecnorm(Q-cc,2,2));
        if radius<.30*min(L(j,:)),continue,end
        if isempty(add)||all(vecnorm(add-cc,2,2)>.70*min(addRadius,radius))
            add(end+1,:)=cc;addRadius(end+1,1)=radius; %#ok<AGROW>
        end
    end
    ns=nnz(split);nSplits=nSplits+ns;
    if ns>0
        for j=find(split).'
            nPhysicalSplits=nPhysicalSplits+ ...
                (segment_distance(cm(j,:),X,physicalEdges)<1e-12);
        end
        ids=size(Q,1)+(1:ns)';Q=[Q;cm(split,:)];
        C=[C(~split,:);C(split,1),ids;ids,C(split,2)];
    end
    if isempty(add)&&ns==0,break,end
    Q=[Q;add];[~,ix,group]=unique(round(Q/1e-12),'rows','stable');Q=Q(ix,:);C=group(C);
    [Q,T,~]=exterior_cells(Q,C,sourceTri,inner);
end

% Split the entire retained crack polyline, including all historical
% vertices and the rear-core intersection. The original copy is upper/left;
% duplicated IDs are assigned to triangles on the lower/right side.
dCut=point_polyline_distance(Q,cutPath);
cut=find(dCut<1e-11);
lower=(1:size(Q,1))';
lower(cut)=size(Q,1)+(1:numel(cut))';

cent=(Q(T(:,1),:)+Q(T(:,2),:)+Q(T(:,3),:))/3;
side=polyline_signed_side(cent,cutPath);
bottom=side<0;
T(bottom,:)=lower(T(bottom,:));

Z=[Q;Q(cut,:)];
nodeSide=ones(size(Z,1),1);
nodeSide(lower(cut))=-1;
crackNodeMask=false(size(Z,1),1);
crackNodeMask(cut)=true;
crackNodeMask(lower(cut))=true;
info=struct('innerPolygon',inner,'nodeSide',nodeSide, ...
    'originalPhysicalIDs',oldIDs,'excludedZeroInteriorHullSlivers',excluded, ...
    'allExteriorInteriorTrianglesReplaced',true,'usesRandomness',false, ...
    'hMax_m',design.farCap_m,'hRp_m',hRp,'farSlope',cal.farSlope, ...
    'slopeTransitionLength_m',cal.transitionLength_m, ...
    'boundaryLayerSeeds',size(layer,1), ...
    'ringSeedCount',nnz(keep),'smoothingSteps',cal.smoothingSteps, ...
    'boundaryMetricGrowth',cal.boundaryMetricGrowth, ...
    'qualityRefinementMaxPasses',cal.refinementMaxPasses, ...
    'exteriorConstraintSubdivisions',nSplits, ...
    'physicalBoundarySubdivisions',nPhysicalSplits, ...
    'originalPhysicalSegments',size(physicalEdges,1), ...
    'retainedCrackPath',cutPath,'crackNodeMask',crackNodeMask, ...
    'refinementMinimumAngle_deg',cal.refinementMinAngle_deg, ...
    'refinementLongestFactor',cal.refinementLongestFactor, ...
    'neighborRatioTarget',cal.neighborRatioTarget, ...
    'calibration',cal);
info.refinementTable=array2table(refinementRows,'VariableNames', ...
    {'pass','exteriorCells','flaggedCells','minimumAngle_deg','maximumNeighborRatio'});
info.ringTable=array2table(ringRows,'VariableNames', ...
    {'radius_m','targetH_m','bandWidth_m','angularIntervals','keptSeeds'});
end

function d=point_polyline_distance(X,P)
if size(X,1)==1
    X=reshape(X,1,2);
end
d=inf(size(X,1),1);
for k=1:size(P,1)-1
    A=P(k,:);B=P(k+1,:);v=B-A;L2=max(dot(v,v),1e-30);
    t=((X-A)*v.')/L2;t=max(0,min(1,t));
    Q=A+t.*v;
    d=min(d,vecnorm(X-Q,2,2));
end
end

function s=polyline_signed_side(X,P)
% Signed distance numerator relative to the nearest retained crack segment.
% Positive means left/upper side of the oriented mouth->tip path.
s=zeros(size(X,1),1);
for i=1:size(X,1)
    best=inf;sgn=0;
    for k=1:size(P,1)-1
        A=P(k,:);B=P(k+1,:);v=B-A;L2=max(dot(v,v),1e-30);
        t=dot(X(i,:)-A,v)/L2;t=max(0,min(1,t));
        Q=A+t*v;w=X(i,:)-Q;d=norm(w);
        if d<best
            best=d;
            sgn=v(1)*w(2)-v(2)*w(1);
        end
    end
    s(i)=sgn;
end
end

function cal=exterior_calibration(design)
assert(isfield(design,'transitionLength_m')&&design.transitionLength_m>0, ...
    'stage3c:MissingTransitionLength','design.transitionLength_m is required.');
assert(isfield(design,'farCap_m')&&design.farCap_m>0, ...
    'stage3c:MissingFarCap','design.farCap_m is required.');
% Closed Step62B C03 calibration transferred nondimensionally:
% transition/a0=1, far cap/a0=0.625, far slope=.10,
% boundary metric growth=.25, neighbor target=1.8.
cal=struct('transitionLength_m',design.transitionLength_m,'farSlope',.10, ...
    'boundaryMetricGrowth',.25,'smoothingSteps',6, ...
    'refinementMaxPasses',140,'refinementMinAngle_deg',25, ...
    'refinementLongestFactor',1.65,'neighborRatioTarget',1.8, ...
    'verbose',true);
if isfield(design,'exteriorCalibration')&&~isempty(design.exteriorCalibration)
    u=design.exteriorCalibration; names=fieldnames(u);
    for k=1:numel(names)
        assert(isfield(cal,names{k}),'stage3c:UnknownExteriorCalibration', ...
            'Unknown exterior calibration field %s.',names{k});
        cal.(names{k})=u.(names{k});
    end
end
assert(design.farCap_m>design.scale*(design.hBase_m+design.slope*design.rCore_m), ...
    'stage3c:FarCapTooSmall','Far-field cap must exceed the core-boundary target size.');
assert(cal.transitionLength_m>0&&cal.farSlope>0&& ...
    cal.boundaryMetricGrowth>=0&&cal.smoothingSteps>=0&& ...
    cal.refinementMaxPasses>=1&&cal.refinementMinAngle_deg>0&& ...
    cal.refinementLongestFactor>1&&cal.neighborRatioTarget>1);
end

function h=exterior_h(r,rp,hRp,farCap,scale,nearSlope,cal)
t=max(0,r-rp);L=cal.transitionLength_m;cap=farCap-hRp;
z=t/L;logcosh=z+log1p(exp(-2*z))-log(2);
increment=scale*(nearSlope*t+(cal.farSlope-nearSlope)*L*logcosh);
h=hRp+cap*tanh(increment/cap);
end
function [angle,longest,L]=triangle_quality(P,T)
a=P(T(:,1),:);b=P(T(:,2),:);c=P(T(:,3),:);
L=[vecnorm(b-c,2,2),vecnorm(a-c,2,2),vecnorm(a-b,2,2)];
angle=Inf(size(T,1),1);
for k=1:3,j=mod(k,3)+1;l=mod(k+1,3)+1;
    v=acosd(max(-1,min(1,(L(:,j).^2+L(:,l).^2-L(:,k).^2)./(2*L(:,j).*L(:,l)))));
    angle=min(angle,v);
end
longest=max(L,[],2);
end
function h=boundary_metric(Z,P,E,growth)
a=P(E(:,1),:);b=P(E(:,2),:);v=b-a;ell=sum(v.^2,2);
h=zeros(size(Z,1),1);
for first=1:500:size(Z,1)
    rows=first:min(first+499,size(Z,1));
    sx=Z(rows,1)-a(:,1).';sy=Z(rows,2)-a(:,2).';
    t=(sx.*v(:,1).'+sy.*v(:,2).')./ell.';t=max(0,min(1,t));
    d2=(sx-t.*v(:,1).').^2+(sy-t.*v(:,2).').^2;
    h(rows)=min(1.2*sqrt(ell).'+growth*sqrt(d2),[],2);
end
end
function [ratio,bad]=neighbor_ratio(T,L,target)
n=size(T,1);E=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
which=repmat((1:n)',3,1);[~,~,g]=unique(E,'rows');
lo=accumarray(g,L(which),[],@min);hi=accumarray(g,L(which),[],@max);
% Refine only the larger neighbor. Refining the small side as well would
% amplify the size jump and create a refinement cascade.
ratio=max(hi./lo);bad=accumarray(which,L(which)>target*lo(g),[n,1],@max)>0;
end
function d=segment_distance(z,P,E)
a=P(E(:,1),:);b=P(E(:,2),:);v=b-a;
t=sum((z-a).*v,2)./sum(v.^2,2);t=max(0,min(1,t));
d=min(vecnorm(a+t.*v-z,2,2));
end
function [Q,T,excluded]=exterior_cells(Q,C,sourceTri,inner)
dt=delaunayTriangulation(Q,C);Q=dt.Points;T=dt.ConnectivityList;
cent=(Q(T(:,1),:)+Q(T(:,2),:)+Q(T(:,3),:))/3;
use=~isnan(pointLocation(sourceTri,cent)) & ...
    ~inpolygon(cent(:,1),cent(:,2),inner(:,1),inner(:,2));
a=Q(T(:,2),:)-Q(T(:,1),:);b=Q(T(:,3),:)-Q(T(:,1),:);
area=.5*(a(:,1).*b(:,2)-a(:,2).*b(:,1));
edgeScale=max([sum(a.^2,2),sum(b.^2,2),sum((a-b).^2,2)],[],2);
degenerate=abs(area)<=1e-12*edgeScale;excluded=nnz(use&degenerate);
T=T(use&~degenerate,:);area=area(use&~degenerate);
T(area<0,[2 3])=T(area<0,[3 2]);
end
