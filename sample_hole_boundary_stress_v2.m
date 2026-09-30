function B=sample_hole_boundary_stress_v2(C,G,S1)
%SAMPLE_HOLE_BOUNDARY_STRESS_V2
% Production Stage-I boundary-stress sampler for a circular hole.
%
% Redesign after SIF-audit Steps 12--14:
%   1) sample strictly on the MATERIAL side of the hole;
%   2) respect the actual FEM topology when interpolating recovered stresses;
%   3) use several offsets eps_j = alpha_j*h_hole;
%   4) linearly extrapolate sigma(phi,eps) to eps -> 0+.
%
% The stress estimator retained for production is the recovered nodal stress
% field from StressExt, interpolated through the actual containing T6 element.
% Direct-from-U remains an audit/reference estimator.
%
% Required C.stage1 fields (defaults shown):
%   nphi             = 1440
%   shift_fractions  = [0.05 0.10 0.25]
%   radial_fit_order = 1
%
% Output B.sig_tt_eff is the extrapolated boundary-limit tangential stress.

must(C,'stage1');
must(S1,'mesh');
must(S1,'stress');

hole=get_single_circular_hole(C,G);

nphi=getf(C.stage1,'nphi',1440);
alpha=getf(C.stage1,'shift_fractions',[0.05 0.10 0.25]);
fitOrder=getf(C.stage1,'radial_fit_order',1);

alpha=alpha(:).';
if numel(alpha)<2
    error('sample_hole_boundary_stress_v2:NeedOffsets', ...
        'At least two material-side offsets are required.');
end
if any(~isfinite(alpha)) || any(alpha<=0) || any(diff(alpha)<=0)
    error('sample_hole_boundary_stress_v2:BadOffsets', ...
        'C.stage1.shift_fractions must be finite, positive, strictly increasing.');
end
if fitOrder~=1
    error('sample_hole_boundary_stress_v2:FitOrder', ...
        'Current production implementation supports radial_fit_order=1 only.');
end

coord=S1.mesh.coord;
connect=S1.mesh.connect;
sNodal=S1.stress;

if size(connect,2)~=6
    error('sample_hole_boundary_stress_v2:NotT6', ...
        'S1.mesh.connect must be T6 connectivity.');
end
if size(sNodal,2)~=3
    error('sample_hole_boundary_stress_v2:BadStress', ...
        'S1.stress must have columns [sxx syy sxy].');
end

c=hole.center(:).';
R=hole.r;

phi=linspace(0,2*pi,nphi+1).';
phi(end)=[];

n_mat=[cos(phi),sin(phi)];
t_hat=[-sin(phi),cos(phi)];
xb=c+R*n_mat;

hhole=get_hhole(C,R);
epsList=alpha*hhole;
nA=numel(alpha);

sig_xx=zeros(nphi,nA);
sig_yy=zeros(nphi,nA);
sig_xy=zeros(nphi,nA);
sig_tt=zeros(nphi,nA);
sig_nn=zeros(nphi,nA);
sig_nt=zeros(nphi,nA);
insideFrac=zeros(1,nA);
parentElem=cell(1,nA);
parentBary=cell(1,nA);
xq=cell(1,nA);

TR=triangulation(G.t,G.p);

for ia=1:nA
    epsShift=epsList(ia);
    Xq=xb+epsShift*n_mat;
    xq{ia}=Xq;

    [elem,bary]=pointLocation(TR,Xq(:,1),Xq(:,2));
    inside=isfinite(elem);
    insideFrac(ia)=mean(inside);

    if any(~inside)
        error('sample_hole_boundary_stress_v2:QueryOutside', ...
            ['Material-side offset alpha=%.6g has %d/%d query points ', ...
             'outside the actual FEM domain.'], ...
             alpha(ia),nnz(~inside),nphi);
    end

    parentElem{ia}=elem;
    parentBary{ia}=bary;

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

        nodes=connect(e,:);
        sv=N*sNodal(nodes,:);

        sig_xx(k,ia)=sv(1);
        sig_yy(k,ia)=sv(2);
        sig_xy(k,ia)=sv(3);

        Sm=[sv(1),sv(3);sv(3),sv(2)];
        n=n_mat(k,:).';
        tt=t_hat(k,:).';

        sig_nn(k,ia)=n.'*Sm*n;
        sig_tt(k,ia)=tt.'*Sm*tt;
        sig_nt(k,ia)=n.'*Sm*tt;
    end
end

% Linear boundary extrapolation in alpha=eps/hhole. Scaling the independent
% variable avoids conditioning problems without changing the intercept.
[sig_xx_b,fit_xx]=linear_intercept(alpha,sig_xx);
[sig_yy_b,fit_yy]=linear_intercept(alpha,sig_yy);
[sig_xy_b,fit_xy]=linear_intercept(alpha,sig_xy);
[sig_tt_b,fit_tt]=linear_intercept(alpha,sig_tt);
[sig_nn_b,fit_nn]=linear_intercept(alpha,sig_nn);
[sig_nt_b,fit_nt]=linear_intercept(alpha,sig_nt);

B=struct();
B.phi=phi;
B.x=xb;
B.n_out=n_mat;       % kept for compatibility
B.n_mat=n_mat;
B.t_hat=t_hat;

B.sig_tt=sig_tt_b;
B.sig_nn=sig_nn_b;
B.sig_nt=sig_nt_b;
B.sig_tt_eff=sig_tt_b;

B.sig_xx=sig_xx_b;
B.sig_yy=sig_yy_b;
B.sig_xy=sig_xy_b;

B.offset=struct();
B.offset.shift_fractions=alpha;
B.offset.eps=epsList;
B.offset.hhole=hhole;
B.offset.xq={xq{:}};
B.offset.sig_xx=sig_xx;
B.offset.sig_yy=sig_yy;
B.offset.sig_xy=sig_xy;
B.offset.sig_tt=sig_tt;
B.offset.sig_nn=sig_nn;
B.offset.sig_nt=sig_nt;
B.offset.fraction_inside_actual_FEM_domain=insideFrac;
B.offset.parentElement={parentElem{:}};
B.offset.parentBary={parentBary{:}};

B.fit=struct();
B.fit.order=fitOrder;
B.fit.sig_xx=fit_xx;
B.fit.sig_yy=fit_yy;
B.fit.sig_xy=fit_xy;
B.fit.sig_tt=fit_tt;
B.fit.sig_nn=fit_nn;
B.fit.sig_nt=fit_nt;

B.hole=hole;
B.method='boundary_extrapolated_recovered_T6';
B.is_boundary_extrapolated=true;
end


function [intercept,F]=linear_intercept(alpha,Y)
% Y is nPoint x nOffset. Fit Y = slope*alpha + intercept at every point.
X=[alpha(:),ones(numel(alpha),1)];
coef=X\Y.'; % 2 x nPoint
slope=coef(1,:).';
intercept=coef(2,:).';

Yfit=(X*coef).';
res=Y-Yfit;

F=struct();
F.slope_per_h=slope;
F.intercept=intercept;
F.rmse=sqrt(mean(res.^2,2));
F.max_abs_residual=max(abs(res),[],2);
end


function h=get_hhole(C,R)
h=[];
if isfield(C,'mesh1')&&isfield(C.mesh1,'hhole')&&~isempty(C.mesh1.hhole)
    h=C.mesh1.hhole;
elseif isfield(C,'mesh1')&&isfield(C.mesh1,'hmin')&&~isempty(C.mesh1.hmin)
    h=C.mesh1.hmin;
end
if isempty(h)
    h=1e-3*R;
end
if ~(isscalar(h)&&isfinite(h)&&h>0)
    error('sample_hole_boundary_stress_v2:BadHhole','Invalid hole mesh scale.');
end
end


function hole=get_single_circular_hole(C,G)
if isfield(G,'hole')&&~isempty(G.hole)
    hole=G.hole;
elseif isfield(C,'hole')&&~isempty(C.hole)
    hole=C.hole;
else
    error('sample_hole_boundary_stress_v2:HoleSpec', ...
        'Exactly one circular hole is required.');
end
if ~isfield(hole,'type')||~strcmpi(strtrim(hole.type),'circle')
    error('sample_hole_boundary_stress_v2:HoleType', ...
        'Current production Stage-I v2 supports a circular hole.');
end
must(hole,'center');
must(hole,'r');
end


function must(S,f)
if ~isfield(S,f)||isempty(S.(f))
    error('sample_hole_boundary_stress_v2:MissingField', ...
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
