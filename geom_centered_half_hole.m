function G=geom_centered_half_hole(C)
%GEOM_CENTERED_HALF_HOLE  Right-half plate with a semicircular hole cutout.
%
% The geometry is the exact symmetry reduction of the centered full plate:
%
%   x in [x_sym,A], y in [-B,B]
%
% with the right semicircle of the centered circular hole removed from the
% left symmetry boundary.  The semicircle is represented by straight
% polygonal segments using half of C.hole.npoly full-circle segments.

must(C,'A'); must(C,'B'); must(C,'hole');
must(C.hole,'center'); must(C.hole,'r'); must(C.hole,'npoly');

A=C.A; B=C.B;
c=C.hole.center(:).';
R=C.hole.r;
Np=C.hole.npoly;

xSym=c(1);

if abs(c(2))>1e-14*max([1,A,B,R])
    error('geom_centered_half_hole:NotCenteredY', ...
        'This benchmark requires hole center y=0.');
end
if abs(xSym-0.5*A)>1e-12*max(1,A)
    error('geom_centered_half_hole:NotCenteredX', ...
        'This benchmark requires hole center x=A/2.');
end
if mod(Np,2)~=0
    error('geom_centered_half_hole:OddNpoly', ...
        'C.hole.npoly must be even so the semicircle uses Npoly/2 segments.');
end
if ~(R>0 && R<B && xSym+R<A)
    error('geom_centered_half_hole:BadGeometry', ...
        'Semicircle must fit strictly inside the right half-plate.');
end

% Arc from top mouth to bottom mouth along the material/cavity boundary.
nArcSeg=Np/2;
thetaArc=linspace(pi/2,-pi/2,nArcSeg+1).';
arc=[c(1)+R*cos(thetaArc), c(2)+R*sin(thetaArc)];

% Counterclockwise simple polygon:
% bottom symmetry corner -> right-bottom -> right-top -> top symmetry
% corner -> top mouth -> clockwise right semicircle -> bottom mouth ->
% implicit close along lower symmetry segment.
outerPoly=[ ...
    xSym,-B;
    A,-B;
    A, B;
    xSym, B;
    arc];

[mdl,dl,bt,gd,ns,sf]=build_pde_geometry_with_holes(outerPoly,{});

mesh1=getf(C,'mesh1',struct());
Hmin=getf(mesh1,'hmin',2*pi*R/Np);
Hmax=getf(mesh1,'hmax',20*Hmin);
Hgrad=getf(mesh1,'hgrad',1.2);

msh=generateMesh(mdl, ...
    'Hmin',Hmin, ...
    'Hmax',Hmax, ...
    'Hgrad',max(1.01,Hgrad), ...
    'GeometricOrder','linear');

p=msh.Nodes.';
t=msh.Elements.';

scale=max([1,A,2*B]);
tol=1e-8*scale;

x=p(:,1); y=p(:,2);
rr=hypot(x-c(1),y-c(2));

edgeSets=struct();
edgeSets.symmetry=find(abs(x-xSym)<tol);
edgeSets.right=find(abs(x-A)<tol);
edgeSets.top=find(abs(y-B)<tol);
edgeSets.bottom=find(abs(y+B)<tol);
edgeSets.semicircle=find(abs(rr-R)<5*tol & x>=xSym-5*tol);
edgeSets.corners=struct();
edgeSets.corners.symmetry_bottom=nearest_node(p,[xSym,-B]);
edgeSets.corners.symmetry_top=nearest_node(p,[xSym,B]);
edgeSets.corners.right_bottom=nearest_node(p,[A,-B]);
edgeSets.corners.right_top=nearest_node(p,[A,B]);

showMesh=false;
if isfield(C,'plot')&&isfield(C.plot,'show_mesh1')
    showMesh=logical(C.plot.show_mesh1);
end

if showMesh
    figure('Name','Centered half-domain Stage I mesh','Color','w'); clf
    hold on; axis equal; box on
    triplot(t,p(:,1),p(:,2),'Color',[0.75 0.75 0.75]);
    plot([outerPoly(:,1);outerPoly(1,1)], ...
         [outerPoly(:,2);outerPoly(1,2)],'k-','LineWidth',1.2);
    plot(arc(:,1),arc(:,2),'r-','LineWidth',1.4);
    plot(p(edgeSets.symmetry,1),p(edgeSets.symmetry,2),'bo','MarkerSize',3);
    xlabel('x'); ylabel('y');
    title('Centered-hole right half-domain');
    xlim([xSym,A]); ylim([-B,B]);
end

G=struct();
G.p=p;
G.t=t;
G.geom=mdl;
G.meshobj=msh;
G.outerPoly=outerPoly;
G.holeLoops={};
G.hole=C.hole;
G.holeArc=arc;
G.holeArcTheta=thetaArc;
G.edgeSets=edgeSets;

G.meta=struct();
G.meta.A=A;
G.meta.B=B;
G.meta.x_sym=xSym;
G.meta.domain_mode='centered_right_half';
G.meta.hole=C.hole;
G.meta.mesh1=mesh1;
G.meta.Hmin=Hmin;
G.meta.Hmax=Hmax;
G.meta.Hgrad=Hgrad;
G.meta.decsg=struct('dl',dl,'bt',bt,'gd',gd,'ns',ns,'sf',sf);
end


function idx=nearest_node(p,xq)
[~,idx]=min(sum((p-xq).^2,2));
end

function must(S,f)
if ~isfield(S,f)||isempty(S.(f))
    error('geom_centered_half_hole:MissingField', ...
        'Required field "%s" is missing or empty.',f);
end
end

function v=getf(S,f,d)
if isstruct(S)&&isfield(S,f)&&~isempty(S.(f))
    v=S.(f);
else
    v=d;
end
end
