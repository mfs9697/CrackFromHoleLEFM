function U=exact_williams_displacement_audit(coord,KI,KII,E,nu,ps)
%EXACT_WILLIAMS_DISPLACEMENT_AUDIT
% Leading-order isotropic Williams displacement field for a crack on x<0.
% Tip is at the origin and local crack direction is +x.

if nargin<6, ps=1; end
mu=E/(2*(1+nu));
if ps==1
    kappa=3-4*nu;
else
    kappa=(3-nu)/(1+nu);
end

n=size(coord,1);
U=zeros(2*n,1);

for i=1:n
    x=coord(i,1); y=coord(i,2);
    r=hypot(x,y);
    th=atan2(y,x);

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
