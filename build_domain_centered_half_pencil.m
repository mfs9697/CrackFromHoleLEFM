function D=build_domain_centered_half_pencil(Pmid,C,w)
hole=C.hole;
c=hole.center(:).';
R=hole.r;
Aplate=C.A;
B=C.B;
xSym=c(1);
A0=[xSym+R,c(2)];
d=Pmid(2,:)-Pmid(1,:);
a0=norm(d);
edir=d/a0;
theta=atan2(edir(2),edir(1));
[~,Gapp]=build_appended_hole_loop(hole,A0,[1,0],[0,1],theta,a0,w, ...
    'epsMode','arclength','nArc',max(80,hole.npoly),'orientation','cw');
Mup=Gapp.Mup;
Mlo=Gapp.Mlo;
xtip=Gapp.xtip;
if Mup(2)<Mlo(2)
    tmp=Mup; Mup=Mlo; Mlo=tmp;
end
Gapp.Mup=Mup;
Gapp.Mlo=Mlo;
Gapp.face_upper=[Mup;xtip];
Gapp.face_lower=[xtip;Mlo];

phiUp=atan2(Mup(2)-c(2),Mup(1)-c(1));
phiLo=atan2(Mlo(2)-c(2),Mlo(1)-c(1));
dphi=2*pi/hole.npoly;

arcTop=local_arc(c,R,pi/2,phiUp,dphi);
arcBot=local_arc(c,R,phiLo,-pi/2,dphi);

outerPoly=[xSym,-B;Aplate,-B;Aplate,B;xSym,B;arcTop;xtip;arcBot];
outerPoly=local_dedup(outerPoly);
if local_area(outerPoly)<0
    outerPoly=flipud(outerPoly);
end

D=struct();
D.outerPoly=outerPoly;
D.holeLoops={};
D.channelPoly=[];
D.Pmid=Pmid;
D.A=Aplate;
D.B=B;
D.w=w;
D.holes={};
D.xSym=xSym;
D.channelGeom=struct('mode','centered_half','append',Gapp);
D.topology=struct('mode','centered_half','xSym',xSym,'A0',A0,'theta',theta,'a0',a0);
end

function P=local_arc(c,R,p1,p2,dpt)
n=max(1,ceil(abs(p2-p1)/dpt));
q=linspace(p1,p2,n+1).';
P=c+R*[cos(q),sin(q)];
end

function P=local_dedup(P)
keep=true(size(P,1),1);
for k=2:size(P,1)
    keep(k)=norm(P(k,:)-P(k-1,:),inf)>1e-13;
end
P=P(keep,:);
end

function A=local_area(P)
x=P(:,1); y=P(:,2);
A=0.5*sum(x.*[y(2:end);y(1)]-[x(2:end);x(1)].*y);
end
