function [r,app,diag]=native_COD_audit(mesh,U,mat,crack,minPts,includeFace)
%NATIVE_COD_AUDIT Shared EXACT Step-38 native crack-face COD extractor.
% Extracted without changing the numerical algorithm. For Step-44 synthetic
% replay, includeFace=true optionally exposes topology-based face labels;
% the default remains compact to avoid saving N-by-1 labels into O38.
if nargin<5 || isempty(minPts),minPts=8;end
if nargin<6 || isempty(includeFace),includeFace=false;end
X=mesh.coord;T=mesh.connect;
n=size(X,1);
tip=crack.Pmid(end,:);
vec=crack.Pmid(end,:)-crack.Pmid(1,:);
a0=norm(vec);vec=vec/a0;
perp=[-vec(2),vec(1)];
R=[vec(:),perp(:)];
xl=(X-tip)*R;
face=zeros(n,1);
tipID=crack.tipNode;
up=unique(crack.upperNodes(:));
lo=unique(crack.lowerNodes(:));
face(setdiff(up,[lo;tipID]))=1;
face(setdiff(lo,[up;tipID]))=-1;
emap=[1 2 4;2 3 5;3 1 6];
for j=1:3
    edge=T(:,emap(j,:));
    v1=edge(:,1);v2=edge(:,2);
    mU=(face(v1)==1&(face(v2)==1|v2==tipID)) | ...
       (face(v2)==1&(face(v1)==1|v1==tipID));
    mL=(face(v1)==-1&(face(v2)==-1|v2==tipID)) | ...
       (face(v2)==-1&(face(v1)==-1|v1==tipID));
    idsU=unique(edge(mU,3));idsL=unique(edge(mL,3));
    if any(face(idsU)==-1) || any(face(idsL)==1)
        error('step38:FaceConflict','Crack-face midside sets conflict.');
    end
    face(idsU)=+1;face(idsL)=-1;
end
tol=max(1e-12,1e-8*a0);
onFace=xl(:,1)<-tol & -xl(:,1)<=a0+tol & abs(xl(:,2))<tol;
up=find(onFace&face==1);lo=find(onFace&face==-1);
if numel(up)<minPts||numel(lo)<minPts
    error('step38:FaceNodes','Too few classified crack-face T6 nodes.');
end
[rU,iu]=sort(-xl(up,1));[rL,il]=sort(-xl(lo,1));
up=up(iu);lo=lo(il);
[rU,~,gU]=unique(rU);[rL,~,gL]=unique(rL);
u=reshape(U,2,[]).'*R;
Uu=zeros(numel(rU),2);Ul=zeros(numel(rL),2);
for k=1:2
    Uu(:,k)=accumarray(gU,u(up,k),[],@mean);
    Ul(:,k)=accumarray(gL,u(lo,k),[],@mean);
end
mask=rU>=min(rL)&rU<=max(rL);
r=rU(mask);
if isempty(r),error('step38:FaceOverlap','No face abscissa overlap.');end
jump=Uu(mask,:) - interp1(rL,Ul,r,'pchip');
if any(~isfinite(jump(:)))
    error('step38:NonfiniteCOD','Native COD interpolation failed.');
end
mu=mat.E/(2*(1+mat.nu));
if mat.ps==1
    kappa=3-4*mat.nu;
else
    kappa=(3-mat.nu)/(1+mat.nu);
end
scale=mu/(kappa+1)*sqrt(2*pi./r);
app=bsxfun(@times,[jump(:,2),jump(:,1)],scale);
mismatch=NaN;
if numel(rU)==numel(rL)
    mismatch=max(abs(rU-rL));
end
diag=struct('nUpper',numel(rU),'nLower',numel(rL), ...
    'gridMismatch',mismatch);
if includeFace,diag.faceSide=face;end
fprintf('  COD new mesh: native nodes upper/lower=%d/%d; ', ...
    diag.nUpper,diag.nLower);
fprintf('abscissa mismatch %.5e m\n',mismatch);
end
