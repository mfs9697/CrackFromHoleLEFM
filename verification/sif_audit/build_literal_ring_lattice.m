function [mesh,info]=build_literal_ring_lattice(varargin)
%BUILD_LITERAL_RING_LATTICE
% Full annular lattice matching the original illustrated mesh concept.
%
% Every circular ring has exactly the same number N of equal chord segments.
% Successive rings alternate angular phase by half a sector:
%   phase_j = 0, dtheta/2, 0, dtheta/2, ...
% The radial growth ratio is chosen from the outward equilateral-triangle
% construction
%   q_eq = cos(dtheta/2) + sqrt(3)*sin(dtheta/2),
% and adjusted slightly so the final ring lands exactly at r1.
%
% This function builds the CLOSED annular lattice only. It does not yet cut
% the negative-x crack seam. The seam must be inserted afterwards so the
% attractive regular ring pattern is not distorted globally.

ip=inputParser;
addParameter(ip,'r0',0.005,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'r1',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'Ntheta',64,@(x)isnumeric(x)&&isscalar(x)&&x>=16);
parse(ip,varargin{:}); S=ip.Results;

N=round(S.Ntheta);
if mod(N,2)~=0, error('Ntheta must be even.'); end
dth=2*pi/N;

qeq=cos(dth/2)+sqrt(3)*sin(dth/2);
Nr=max(2,round(log(S.r1/S.r0)/log(qeq)));
q=(S.r1/S.r0)^(1/Nr);
rv=S.r0*q.^(0:Nr);
rv(end)=S.r1;

coord=zeros((Nr+1)*N,2);
rings=cell(Nr+1,1);
next=0;
for j=1:Nr+1
    phase=mod(j-1,2)*dth/2;
    ids=zeros(1,N);
    for k=0:N-1
        next=next+1;
        th=phase+k*dth;
        coord(next,:)=rv(j)*[cos(th),sin(th)];
        ids(k+1)=next;
    end
    rings{j}=ids;
end

connect=zeros(2*Nr*N,3);
e=0;
for j=1:Nr
    in=rings{j}; out=rings{j+1};
    if mod(j-1,2)==0
        for k=1:N
            kp=mod(k,N)+1;
            e=e+1; connect(e,:)=[in(k),in(kp),out(k)];
            e=e+1; connect(e,:)=[in(kp),out(kp),out(k)];
        end
    else
        for k=1:N
            kp=mod(k,N)+1;
            e=e+1; connect(e,:)=[in(k),in(kp),out(kp)];
            e=e+1; connect(e,:)=[in(k),out(kp),out(k)];
        end
    end
end

A=tri_area_signed(connect,coord);
cw=A<0;
if any(cw)
    tmp=connect(cw,2);
    connect(cw,2)=connect(cw,3);
    connect(cw,3)=tmp;
end

Q=triangle_quality(connect,coord);
segSpread=ring_segment_rel_spread(coord,rings);

mesh=struct('coord3',coord,'connect3',connect,'rings',{rings});
info=struct();
info.r0=S.r0; info.r1=S.r1; info.Ntheta=N; info.Nr=Nr;
info.dtheta=dth; info.qEquilateral=qeq; info.qActual=q;
info.radii=rv;
info.ringSegmentRelSpreadMax=segSpread;
info.qualityMin=min(Q);
info.qualityMedian=median(Q);
info.qualityP05=local_percentile(Q,5);
info.nT3Vertices=size(coord,1);
info.nT3Elements=size(connect,1);
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
smax=0;
for j=1:numel(rings)
    ids=rings{j};
    P=X([ids,ids(1)],:);
    L=sqrt(sum(diff(P,1,1).^2,2));
    sm=mean(L);
    if sm>0
        smax=max(smax,(max(L)-min(L))/sm);
    end
end
end

function y=local_percentile(x,p)
x=sort(x(:));
z=1+(numel(x)-1)*p/100;
i=floor(z); j=ceil(z);
if i==j, y=x(i); else, y=x(i)+(z-i)*(x(j)-x(i)); end
end
