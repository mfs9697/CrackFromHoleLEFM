function B=sample_hole_boundary_stress_audit(C,G,S1,varargin)
%SAMPLE_HOLE_BOUNDARY_STRESS_AUDIT
% Audit-only variant of Stage-I hole-boundary stress sampling.
%
% It keeps the existing recovered nodal stresses from StressExt and lets us
% isolate two later choices:
%   QuerySide:
%     'legacy_cavity' : xq = xb - eps*n_mat  (current production behavior)
%     'material'      : xq = xb + eps*n_mat  (query inside the solid)
%   Interpolator:
%     'scattered'     : current scatteredInterpolant over all T6 nodes
%     't6'            : topology-respecting T6 interpolation of the same
%                       recovered nodal stresses using the actual FEM parent
%                       T3 element containing xq
%
% No stress recovery is changed here; only query location/interpolation are
% varied for diagnosis.

ip=inputParser;
addParameter(ip,'QuerySide','legacy_cavity', ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Interpolator','scattered', ...
    @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'ShiftFraction',0.25, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
parse(ip,varargin{:});
O=ip.Results;

querySide=lower(strtrim(char(O.QuerySide)));
interpMode=lower(strtrim(char(O.Interpolator)));

hole=G.hole;
nphi=C.stage1.nphi;
coord=S1.mesh.coord;
sNodal=S1.stress;

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

switch querySide
    case 'legacy_cavity'
        xq=xb-eps_shift*n_mat;
    case 'material'
        xq=xb+eps_shift*n_mat;
    otherwise
        error('sample_hole_boundary_stress_audit:BadQuerySide', ...
            'Unknown QuerySide %s.',querySide);
end

% Actual FEM-domain membership, independent of interpolation choice.
TR=triangulation(G.t,G.p);
[elem,bary]=pointLocation(TR,xq(:,1),xq(:,2));
inside=isfinite(elem);

switch interpMode
    case 'scattered'
        Fx=scatteredInterpolant(coord(:,1),coord(:,2),sNodal(:,1),'linear','nearest');
        Fy=scatteredInterpolant(coord(:,1),coord(:,2),sNodal(:,2),'linear','nearest');
        Fxy=scatteredInterpolant(coord(:,1),coord(:,2),sNodal(:,3),'linear','nearest');
        sig_xx=Fx(xq(:,1),xq(:,2));
        sig_yy=Fy(xq(:,1),xq(:,2));
        sig_xy=Fxy(xq(:,1),xq(:,2));

    case 't6'
        if any(~inside)
            error('sample_hole_boundary_stress_audit:T6Outside', ...
                ['Topology-respecting T6 interpolation requested, but %d/%d ', ...
                 'query points lie outside the actual FEM domain.'], ...
                 nnz(~inside),numel(inside));
        end
        sig_xx=zeros(nphi,1);
        sig_yy=zeros(nphi,1);
        sig_xy=zeros(nphi,1);
        for k=1:nphi
            e=elem(k);
            lam=bary(k,:);
            N=[ ...
                lam(1)*(2*lam(1)-1), ...
                lam(2)*(2*lam(2)-1), ...
                lam(3)*(2*lam(3)-1), ...
                4*lam(1)*lam(2), ...
                4*lam(2)*lam(3), ...
                4*lam(3)*lam(1)];
            nodes=S1.mesh.connect(e,:);
            sv=N*sNodal(nodes,:);
            sig_xx(k)=sv(1);
            sig_yy(k)=sv(2);
            sig_xy(k)=sv(3);
        end

    otherwise
        error('sample_hole_boundary_stress_audit:BadInterpolator', ...
            'Unknown Interpolator %s.',interpMode);
end

sig_tt=zeros(nphi,1);
sig_nn=zeros(nphi,1);
sig_nt=zeros(nphi,1);
for k=1:nphi
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
B.n_out=n_mat; % same field name as current sampler
B.t_hat=t_hat;
B.sig_tt=sig_tt;
B.sig_nn=sig_nn;
B.sig_nt=sig_nt;
B.sig_tt_eff=sig_tt;
B.hole=hole;
B.eps_shift=eps_shift;
B.querySide=querySide;
B.interpolator=interpMode;
B.insideActualFEMDomain=inside;
B.fractionInsideActualFEMDomain=mean(inside);
B.parentElement=elem;
B.parentBary=bary;
end
