function B=sample_hole_boundary_stress_direct_U_audit(C,G,S1,varargin)
%SAMPLE_HOLE_BOUNDARY_STRESS_DIRECT_U_AUDIT
% Evaluate Stage-I stress directly from the T6 displacement field:
%     eps(xq) = B(xq) * Ue
%     sig(xq) = D * eps(xq)
%
% Query points are on the material side of the circular hole:
%     xq = xb + ShiftFraction * hhole * n_mat.
%
% No GP-stress extrapolation, nodal stress recovery, nodal averaging, or
% stress-field interpolation is used.

ip=inputParser;
addParameter(ip,'ShiftFraction',0.25, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
parse(ip,varargin{:});
O=ip.Results;

hole=G.hole;
nphi=C.stage1.nphi;

c=hole.center(:).';
R=hole.r;
phi=linspace(0,2*pi,nphi+1).';
phi(end)=[];

n_mat=[cos(phi),sin(phi)];
t_hat=[-sin(phi),cos(phi)];
xb=c+R*n_mat;

hhole=[];
if isfield(C,'mesh1')&&isfield(C.mesh1,'hhole')
    hhole=C.mesh1.hhole;
end
if isempty(hhole)
    hhole=C.mesh1.hmax;
end

eps_shift=O.ShiftFraction*hhole;
eps_shift=min(eps_shift,0.10*R);
eps_shift=max(eps_shift,1e-8*max(R,1));

xq=xb+eps_shift*n_mat;

TR=triangulation(G.t,G.p);
[elem,bary]=pointLocation(TR,xq(:,1),xq(:,2));
inside=isfinite(elem);

if any(~inside)
    error('sample_hole_boundary_stress_direct_U_audit:Outside', ...
        '%d/%d material-side query points lie outside the FEM domain.', ...
        nnz(~inside),numel(inside));
end

coord=S1.mesh.coord;
conn=S1.mesh.connect;
U=S1.U;
D=S1.mat.D;

nq=numel(phi);
sig_xx=zeros(nq,1);
sig_yy=zeros(nq,1);
sig_xy=zeros(nq,1);
detJ=zeros(nq,1);

for k=1:nq
    e=elem(k);
    nodes=conn(e,:);
    X=coord(nodes,:);

    % pointLocation returns barycentric coordinates associated with the T3
    % corner nodes; BN_local uses the first two and reconstructs the third.
    xi0=bary(k,1:2).';
    [Bmat,detJ(k)]=BN_local_direct(xi0,X);

    eldof=zeros(12,1);
    for a=1:6
        eldof(2*a-1:2*a)=[2*nodes(a)-1;2*nodes(a)];
    end

    epsv=Bmat*U(eldof);
    sigv=D*epsv;

    sig_xx(k)=sigv(1);
    sig_yy(k)=sigv(2);
    sig_xy(k)=sigv(3);
end

sig_tt=zeros(nq,1);
sig_nn=zeros(nq,1);
sig_nt=zeros(nq,1);

for k=1:nq
    Sm=[sig_xx(k),sig_xy(k);sig_xy(k),sig_yy(k)];
    n=n_mat(k,:).';
    tt=t_hat(k,:).';

    sig_nn(k)=n.'*Sm*n;
    sig_tt(k)=tt.'*Sm*tt;
    sig_nt(k)=n.'*Sm*tt;
end

B=struct();
B.phi=phi;
B.x=xb;
B.xq=xq;
B.n_out=n_mat;
B.t_hat=t_hat;
B.sig_tt=sig_tt;
B.sig_nn=sig_nn;
B.sig_nt=sig_nt;
B.sig_tt_eff=sig_tt;
B.hole=hole;
B.eps_shift=eps_shift;
B.querySide='material';
B.estimator='direct_U';
B.fractionInsideActualFEMDomain=mean(inside);
B.parentElement=elem;
B.parentBary=bary;
B.detJ=detJ;
end


function [B,Det]=BN_local_direct(xi0,X)
% Exact copy of the T6 strain-displacement construction used in StressExt.
xi=[xi0;1-sum(xi0)];

Nap=[ ...
    4*xi(1)-1, 0,           1-4*xi(3), 4*xi(2),          -4*xi(2),         4*xi(3)-4*xi(1);
    0,           4*xi(2)-1, 1-4*xi(3), 4*xi(1),           4*xi(3)-4*xi(2), -4*xi(1)];

dxdxi=Nap*X;
Det=det(dxdxi);

N1=[ dxdxi(2,2),-dxdxi(1,2);
    -dxdxi(2,1), dxdxi(1,1)]/Det*Nap;

eldf=12;
inx=(2:2:eldf)'-1;
iny=inx+1;

B=zeros(3,eldf);
B(1,inx)=N1(1,:);
B(2,iny)=N1(2,:);
B(3,inx)=N1(2,:);
B(3,iny)=N1(1,:);
end
