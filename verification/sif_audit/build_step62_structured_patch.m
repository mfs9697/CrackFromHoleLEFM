function [Z,T,mirror,nUpper,rings,axisIDs,design]=build_step62_structured_patch(rp,hTip,level)
%BUILD_STEP62_STRUCTURED_PATCH Explicit deterministic polar T3 topology.
% No Delaunay, random points, smoothing, or historical refinement pattern.
% Three 60-degree upper sectors define six equilateral tip triangles.
% Ring subdivisions are multiples of three. Adjacent sectors are joined
% by an ordered zipper with a shortest-diagonal rule and fixed tie breaks.
% One upper connectivity is reflected; the intact positive axis is shared.
if nargin<3,level=0;end
validateattributes(rp,{'double'},{'scalar','positive','finite'});
validateattributes(hTip,{'double'},{'scalar','positive','finite'});
validateattributes(level,{'double'},{'scalar','integer','nonnegative','<=',4});
scale=2^(-level);rFirst=scale*hTip;
assert(rp>4*rFirst,'step62:PatchSize','Insufficient room for graded rings.');
slope=.028; radialFactor=sqrt(3)/2;
% h(r)=scale*(hTip+slope*r), increasing everywhere, including beyond ro.
% Equal increments in integral dr/h give strictly increasing band widths.
% Fit the outer boundary exactly, without a short arbitrary last band.
metricSpan=log((hTip+slope*rp)/(hTip+slope*rFirst))/slope;
nBands=ceil(metricSpan/(radialFactor*scale));ds=metricSpan/nBands;
k=(0:nBands)';
rings=((hTip+slope*rFirst)*exp(slope*ds*k)-hTip)/slope;
rings(1)=rFirst;rings(end)=rp;
h=scale*(hTip+slope*rings);
nTheta=3*ceil(pi*rings./h/3);nTheta(1)=3;
assert(all(diff(nTheta)>=0)&&all(mod(nTheta,3)==0));
assert(all(diff(diff(rings))>=-1e-15),'step62:NonmonotoneRadial','Ring widths decrease.');
U=[0 0]; ringIDs=cell(numel(rings),1);nodeRing=0;nodeAngle=0;
for j=1:numel(rings)
    th=(0:nTheta(j))'*pi/nTheta(j);
    xy=rings(j)*[cos(th),sin(th)];xy([1 end],2)=0;
    ringIDs{j}=(size(U,1)+(1:size(xy,1)))';
    U=[U;xy];nodeRing=[nodeRing;repmat(j,size(xy,1),1)]; ...
        nodeAngle=[nodeAngle;th]; %#ok<AGROW>
end
first=ringIDs{1};Tu=[ones(3,1),first(1:3),first(2:4)];
bandRange=zeros(nBands+1,2);bandRange(1,:)=[1 3];
for j=2:numel(rings)
    start=size(Tu,1)+1;
    A=ringIDs{j-1};B=ringIDs{j};na=nTheta(j-1)/3;nb=nTheta(j)/3;
    for sector=0:2
        aa=A(sector*na+(1:na+1));bb=B(sector*nb+(1:nb+1));
        i=1;l=1;
        while i<=na||l<=nb
            if i>na,advanceInner=false;
            elseif l>nb,advanceInner=true;
            else
                di=sum((U(aa(i+1),:)-U(bb(l),:)).^2);
                dout=sum((U(aa(i),:)-U(bb(l+1),:)).^2);
                tol=64*eps(max(di,dout));
                if abs(di-dout)<=tol
                    advanceInner=mod(j+sector+i+l,2)==0;
                else,advanceInner=di<dout;end
            end
            if advanceInner
                Tu(end+1,:)=[aa(i),bb(l),aa(i+1)];i=i+1; %#ok<AGROW>
            else
                Tu(end+1,:)=[aa(i),bb(l),bb(l+1)];l=l+1; %#ok<AGROW>
            end
        end
    end
    bandRange(j,:)=[start,size(Tu,1)];
end
a=U(Tu(:,2),:)-U(Tu(:,1),:);b=U(Tu(:,3),:)-U(Tu(:,1),:);
assert(all(a(:,1).*b(:,2)-a(:,2).*b(:,1)>0),'step62:ExplicitOrientation','Invalid strip.');
assert(size(Tu,1)==3+sum(nTheta(1:end-1)+nTheta(2:end)));
nUpper=size(U,1);axisIDs=find(U(:,2)==0);
shared=axisIDs(U(axisIDs,1)>=0);nonShared=setdiff((1:nUpper)',shared);
mirror=(1:nUpper)';mirror(nonShared)=nUpper+(1:numel(nonShared))';
Z=[U;U(nonShared,1),-U(nonShared,2)];
T=[Tu;mirror(Tu(:,[1 3 2]))];
design=struct('family','explicit three-sector radial zipper', ...
    'level',level,'scale',scale,'hTipTarget_m',rFirst,'hBase_m',hTip, ...
    'slope',slope,'radialFactor',radialFactor,'metricStep',ds, ...
    'upperRingIDs',{ringIDs},'upperNodeRing',nodeRing,'upperNodeAngle',nodeAngle, ...
    'upperBandElementRange',bandRange,'angularIntervalsUpper',nTheta, ...
    'tipTriangles',6,'tipIncidentTopologicalEdges',7, ...
    'hLaw','h(r)=2^(-level)*(hTip+0.028*r)', ...
    'connectivityRule','three sector strips; shortest new bridge; fixed parity ties', ...
    'usesDelaunayInPairedRegion',false,'usesRandomness',false, ...
    'usesSmoothingInPairedRegion',false,'allBandWidthsIncreasing',true);
width=[NaN;diff(rings)];
design.ringTable=table((1:numel(rings))',rings*1e3,h*1e3,width*1e3, ...
    nTheta,pi*1e3*rings./nTheta, ...
    'VariableNames',{'ring','radius_mm','targetH_mm','bandWidth_mm', ...
    'upperAngularIntervals','arcSpacing_mm'});
end
