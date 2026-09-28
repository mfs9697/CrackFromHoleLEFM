function U=exact_williams_displacement_audit(coord,KI,KII,E,nu,ps,varargin)
%EXACT_WILLIAMS_DISPLACEMENT_AUDIT
% Leading-order isotropic Williams displacement field for a crack on x<0.
% Tip is at the origin and local crack direction is +x.
% For a mesh with duplicated coincident crack-face nodes on x<0, pass
%   'UpperFaceIDs', idsUpper, 'LowerFaceIDs', idsLower
% so the branch angle is forced to +pi on the upper face and -pi on the
% lower face. Coordinates alone cannot distinguish coincident face nodes.

if nargin<6 || isempty(ps), ps=1; end
ip=inputParser;
addParameter(ip,'UpperFaceIDs',[],@(x)isempty(x)||(isnumeric(x)&&isvector(x)));
addParameter(ip,'LowerFaceIDs',[],@(x)isempty(x)||(isnumeric(x)&&isvector(x)));
parse(ip,varargin{:});
upperFace=unique(ip.Results.UpperFaceIDs(:));
lowerFace=unique(ip.Results.LowerFaceIDs(:));

n=size(coord,1);
if ~isempty([upperFace;lowerFace])
    ids=[upperFace;lowerFace];
    if any(ids<1|ids>n|ids~=fix(ids)) || numel(unique(ids))~=numel(ids)
        error('exact_williams_displacement_audit:BadFaceIDs', ...
            'Upper/lower crack-face IDs must be valid and disjoint.');
    end
    tol=128*eps(max(1,max(abs(coord(:)))));
    rr=hypot(coord(ids,1),coord(ids,2));
    badFace=abs(coord(ids,2))>tol | ((coord(ids,1)>=-tol) & (rr>tol));
    if any(badFace)
        error('exact_williams_displacement_audit:BadFaceGeometry', ...
            'Crack-face IDs must lie on x<=0, y=0; the origin is allowed.');
    end
end
mu=E/(2*(1+nu));
if ps==1
    kappa=3-4*nu;
else
    kappa=(3-nu)/(1+nu);
end

U=zeros(2*n,1);

isUpper=false(n,1); isLower=false(n,1);
isUpper(upperFace)=true; isLower(lowerFace)=true;

for i=1:n
    x=coord(i,1); y=coord(i,2);
    r=hypot(x,y);
    if isUpper(i)
        th=pi;
    elseif isLower(i)
        th=-pi;
    else
        th=atan2(y,x);
    end

    fac=sqrt(r/(2*pi))/(2*mu);
    c=cos(th/2); s=sin(th/2);

    u1I=KI*fac*c*(kappa-1+2*s^2);
    u2I=KI*fac*s*(kappa+1-2*c^2);

    u1II=KII*fac*s*(kappa+1+2*c^2);
    u2II=-KII*fac*c*(kappa-1-2*s^2);

    U(2*i-1)=u1I+u1II;
    U(2*i)=u2I+u2II;
end
end
