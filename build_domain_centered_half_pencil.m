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
D=struct();
end
