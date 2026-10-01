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
D=struct();
end
