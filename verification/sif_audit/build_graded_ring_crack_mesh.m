function [mesh,info]=build_graded_ring_crack_mesh(varargin)
%BUILD_GRADED_RING_CRACK_MESH
% Graded concentric-ring crack-tip mesh with uniform subdivision of every
% ring and near-equilateral T3 triangles.
%
% IMPORTANT: every half-ring is divided into equal angular segments. To
% interlace adjacent rings without creating half-size seam segments, the
% number of half-ring intervals alternates M, M+1, M, M+1, ... . Thus both
% x-axes are nodes of every ring and the negative-x crack seam remains clean.
% The radial ratio is chosen from the equilateral-altitude estimate
%   q_target = 1 + sqrt(3)*sin(dtheta/2),
% then adjusted slightly so the final ring lands exactly at r1.
%
% Variants:
%   mirror_reflected : exact reflected lower half (canonical symmetric case)
%   lower_shift      : same density; lower interior angles shifted smoothly
%   lower_coarse     : lower half uses fewer angular sectors
%
% Name-value:
%   r0                 0.005
%   r1                 0.20
%   Ntheta             64 full-circle nominal angular divisions
%   Variant            'mirror_reflected'
%   LowerShiftFraction 0.20, measured in upper-sector widths
%   LowerCoarsenFactor 1.5
%
% The positive x-axis is shared by upper/lower halves. The negative x-axis
% is the crack: upper/lower face nodes are geometrically coincident but have
% distinct node IDs.

ip=inputParser;
addParameter(ip,'r0',0.005,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'r1',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'Ntheta',64,@(x)isnumeric(x)&&isscalar(x)&&x>=16);
addParameter(ip,'Variant','mirror_reflected',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'LowerShiftFraction',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<0.45);
addParameter(ip,'LowerCoarsenFactor',1.5,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
parse(ip,varargin{:}); S=ip.Results;

if S.r1<=S.r0, error('Require r1>r0.'); end
Ntheta=round(S.Ntheta);
if mod(Ntheta,2)~=0, error('Ntheta must be even.'); end
Nh=Ntheta/2;
dth=pi/Nh;

% Target radial growth from the altitude of an equilateral triangle built
% on one inner-ring chord. This is a shape target, not a claim that every
% annular triangle can be exactly equilateral.
qtarget=1 + sqrt(3)*sin(dth/2);
Nr=max(2,round(log(S.r1/S.r0)/log(qtarget)));
q=(S.r1/S.r0)^(1/Nr);
rv=S.r0*q.^(0:Nr);
rv(end)=S.r1;

% Upper half.
[coordU,ringsU,thetaU]=make_half(rv,Nh,+1,0);
Tup=connect_half(ringsU,thetaU);

variant=lower(char(S.Variant));

if strcmp(variant,'mirror_reflected')
    coord3=coordU;
    mirror=zeros(size(coordU,1),1);

    % theta=0 nodes are shared.
    for j=1:numel(ringsU)
        mirror(ringsU{j}(1))=ringsU{j}(1);
    end

    n0=size(coord3,1);
    next=n0;
    for j=1:numel(ringsU)
        ids=ringsU{j};
        for k=2:numel(ids)
            iu=ids(k);
            next=next+1;
            coord3(next,:)=[coordU(iu,1),-coordU(iu,2)]; %#ok<AGROW>
            mirror(iu)=next;
        end
    end

    Tlo=mirror(Tup);
    Tlo=Tlo(:,[1 3 2]);
    connect3=[Tup;Tlo];

    crackUpper=cellfun(@(x)x(end),ringsU).';
    crackLower=mirror(crackUpper);

    actualLowerFactor=1;
    lowerShiftFraction=0;
else
    if strcmp(variant,'lower_shift')
        NhL=Nh;
        shift=S.LowerShiftFraction;
        actualLowerFactor=1;
        lowerShiftFraction=shift;
    elseif strcmp(variant,'lower_coarse')
        NhL=max(4,round(Nh/S.LowerCoarsenFactor));
        shift=0;
        actualLowerFactor=Nh/NhL;
        lowerShiftFraction=0;
    else
        error('Unknown Variant: %s',variant);
    end

    % Build lower half independently, then merge its theta=0 nodes into the
    % already-existing upper positive-x radial line.
    [coordL,ringsL,thetaL]=make_half(rv,NhL,-1,shift);
    TloLocal=connect_half(ringsL,thetaL);

    map=zeros(size(coordL,1),1);
    coord3=coordU;
    next=size(coord3,1);

    for j=1:numel(ringsL)
        ids=ringsL{j};
        map(ids(1))=ringsU{j}(1); % share theta=0
        for k=2:numel(ids)
            next=next+1;
            coord3(next,:)=coordL(ids(k),:); %#ok<AGROW>
            map(ids(k))=next;
        end
    end

    Tlo=map(TloLocal);
    connect3=[Tup;Tlo];

    crackUpper=cellfun(@(x)x(end),ringsU).';
    crackLower=zeros(numel(ringsL),1);
    for j=1:numel(ringsL), crackLower(j)=map(ringsL{j}(end)); end
end

% Enforce CCW and reject invalid elements.
A=tri_area_signed(connect3,coord3);
cw=A<0;
if any(cw)
    tmp=connect3(cw,2);
    connect3(cw,2)=connect3(cw,3);
    connect3(cw,3)=tmp;
end
A=tri_area_signed(connect3,coord3);
if any(A<=0), error('Non-positive T3 area in graded-ring mesh.'); end

[coord6,connect6]=T3toT6_fast(coord3,connect3);

% Triangle quality: 1 for equilateral.
Q=triangle_quality(connect3,coord3);

mesh=struct('coord3',coord3,'connect3',connect3, ...
            'coord',coord6,'connect',connect6);

info=struct();
info.variant=variant;
info.r0=S.r0; info.r1=S.r1;
info.Ntheta=Ntheta; info.NhalfUpper=Nh;
info.Nr=Nr; info.radii=rv;
info.qTargetNearEquilateral=qtarget; info.qActual=q;
info.relativeQMismatch=(q-qtarget)/qtarget;
% Backward-compatible alias for older audit scripts; do not interpret as
% exact equilateral geometry.
info.qEquilateral=qtarget;
info.lowerAngularFactor=actualLowerFactor;
info.lowerShiftFraction=lowerShiftFraction;
info.nT3Vertices=size(coord3,1);
info.nT3Elements=size(connect3,1);
info.nT6Nodes=size(coord6,1);
info.nT6Elements=size(connect6,1);
info.nDOF=2*size(coord6,1);
info.qualityMin=min(Q);
info.qualityMedian=median(Q);
info.qualityP05=local_percentile(Q,5);
info.qualityP95=local_percentile(Q,95);
info.upperRingSegmentRelSpreadMax=ring_segment_rel_spread(coord3,ringsU);
if exist('ringsL','var')
    info.lowerRingSegmentRelSpreadMax=ring_segment_rel_spread(coordL,ringsL);
else
    info.lowerRingSegmentRelSpreadMax=info.upperRingSegmentRelSpreadMax;
end
info.nQualityBelow08=nnz(Q<0.8);
info.crackUpperIDs=crackUpper;
info.crackLowerIDs=crackLower;
info.crackFacesDistinct=all(crackUpper~=crackLower);

if strcmp(variant,'mirror_reflected')
    nU=size(coordU,1);
    ids=(1:nU).';
    mapped=mirror(ids);
    xr=[coord3(ids,1),-coord3(ids,2)];
    info.maxMirrorCoordError=max(vecnorm(coord3(mapped,:)-xr,2,2));
else
    info.maxMirrorCoordError=NaN;
end
end

function [coord,rings,thetaR]=make_half(rv,Nh,sgn,shiftFraction)
% Every ring is uniformly divided. Adjacent rings alternate between Nh and
% Nh+1 half-ring intervals so their nodes interlace without forcing a
% half-sector at theta=0 or theta=pi.

dth0=pi/Nh;
rings=cell(numel(rv),1);
thetaR=cell(numel(rv),1);
coord=zeros(0,2);

for j=1:numel(rv)
    nseg=Nh + mod(j-1,2);
    th=linspace(0,pi,nseg+1);

    if shiftFraction>0
        % Deliberate asymmetry variant only. Boundaries stay fixed.
        th=th + shiftFraction*dth0*sin(th);
        th(1)=0; th(end)=pi;
    end

    ids=zeros(1,numel(th));
    for k=1:numel(th)
        ids(k)=size(coord,1)+1;
        coord(ids(k),:)=[rv(j)*cos(th(k)),sgn*rv(j)*sin(th(k))];
    end
    rings{j}=ids;
    thetaR{j}=th;
end
end

function T=connect_half(rings,thetaR)
T=zeros(0,3);
for j=1:numel(rings)-1
    in=rings{j}; out=rings{j+1};
    a=thetaR{j}; b=thetaR{j+1};
    i=1; k=1;
    while i<numel(in) || k<numel(out)
        if i==numel(in)
            T(end+1,:)=[in(i),out(k),out(k+1)]; %#ok<AGROW>
            k=k+1;
        elseif k==numel(out)
            T(end+1,:)=[in(i),in(i+1),out(k)]; %#ok<AGROW>
            i=i+1;
        elseif a(i+1) <= b(k+1)
            T(end+1,:)=[in(i),in(i+1),out(k)]; %#ok<AGROW>
            i=i+1;
        else
            T(end+1,:)=[in(i),out(k),out(k+1)]; %#ok<AGROW>
            k=k+1;
        end
    end
end
end

function A=tri_area_signed(T,X)
v1=X(T(:,1),:); v2=X(T(:,2),:); v3=X(T(:,3),:);
A=0.5*((v2(:,1)-v1(:,1)).*(v3(:,2)-v1(:,2)) - ...
       (v2(:,2)-v1(:,2)).*(v3(:,1)-v1(:,1)));
end

function Q=triangle_quality(T,X)
a=sqrt(sum((X(T(:,2),:)-X(T(:,1),:)).^2,2));
b=sqrt(sum((X(T(:,3),:)-X(T(:,2),:)).^2,2));
c=sqrt(sum((X(T(:,1),:)-X(T(:,3),:)).^2,2));
A=abs(tri_area_signed(T,X));
Q=4*sqrt(3)*A./(a.^2+b.^2+c.^2);
end


function smax=ring_segment_rel_spread(X,rings)
% Maximum within-ring chord-length spread, normalized by ring mean.
smax=0;
for j=1:numel(rings)
    ids=rings{j};
    P=X(ids,:);
    L=sqrt(sum(diff(P,1,1).^2,2));
    if isempty(L), continue; end
    sm=mean(L);
    if sm>0
        smax=max(smax,(max(L)-min(L))/sm);
    end
end
end

function y=local_percentile(x,p)
x=sort(x(:));
if isempty(x), y=NaN; return; end
z=1+(numel(x)-1)*p/100;
i=floor(z); j=ceil(z);
if i==j, y=x(i); else, y=x(i)+(z-i)*(x(j)-x(i)); end
end
