function I=find_hole_initiation_point_v2(C,B)
%FIND_HOLE_INITIATION_POINT_V2
% Determine Stage-I initiation from the extrapolated boundary stress field.
%
% Redesign features:
%   - uses B.sig_tt_eff from material-side eps->0 extrapolation;
%   - finds the discrete tensile maximum;
%   - optionally refines its angular position by a local periodic quadratic fit;
%   - evaluates the circular-hole point/frame at the fitted angle.
%
% Config defaults:
%   C.stage1.angular_fit_enable = true
%   C.stage1.angular_fit_points = 5  (odd integer >=3)

must(C,'sig_c');
must(C,'stage1');
must(B,'phi');
must(B,'sig_tt_eff');
must(B,'hole');

phi=B.phi(:);
sig=B.sig_tt_eff(:);
sigPos=max(sig,0);

[sigDisc,idxDisc]=max(sigPos);
if ~(isfinite(sigDisc)&&sigDisc>0)
    error('find_hole_initiation_point_v2:NoTension', ...
        'No positive finite boundary tangential stress was found.');
end

doFit=logical(getf(C.stage1,'angular_fit_enable',true));
nFit=round(getf(C.stage1,'angular_fit_points',5));
if mod(nFit,2)==0 || nFit<3
    error('find_hole_initiation_point_v2:AngularFitPoints', ...
        'C.stage1.angular_fit_points must be an odd integer >=3.');
end

phiDisc=phi(idxDisc);
phiStar=phiDisc;
sigStar=sigDisc;
fitAccepted=false;
fitInfo=struct('accepted',false,'reason','disabled', ...
    'indices',[],'dphi',[],'values',[],'coefficients',nan(1,3), ...
    'curvature',NaN,'vertex_offset',NaN,'rmse',NaN);

if doFit
    useScaledWindow=isfield(C.stage1,'angular_fit_halfwidth_factor') && ...
        ~isempty(C.stage1.angular_fit_halfwidth_factor);

    if useScaledWindow
        factor=C.stage1.angular_fit_halfwidth_factor;
        if ~(isscalar(factor)&&isfinite(factor)&&factor>0)
            error('find_hole_initiation_point_v2:AngularHalfwidthFactor', ...
                'C.stage1.angular_fit_halfwidth_factor must be a positive scalar.');
        end
        if ~isfield(B,'offset') || ~isfield(B.offset,'hhole') || ...
                ~isfield(B,'hole') || ~isfield(B.hole,'r')
            error('find_hole_initiation_point_v2:MissingMeshScale', ...
                'Mesh-scaled angular fitting requires B.offset.hhole and B.hole.r.');
        end
        halfWidth=factor*B.offset.hhole/B.hole.r;
        isPeriodic=logical(getf(B,'periodic',true));
        [phiFit,sigFit,fitInfo]=window_quadratic_peak( ...
            phi,sigPos,idxDisc,halfWidth,isPeriodic);
        fitInfo.halfwidth_factor=factor;
        fitInfo.halfwidth_rad=halfWidth;
    else
        [phiFit,sigFit,fitInfo]=periodic_quadratic_peak(phi,sigPos,idxDisc,nFit);
        fitInfo.halfwidth_factor=NaN;
        fitInfo.halfwidth_rad=NaN;
    end

    if fitInfo.accepted && isfinite(sigFit) && sigFit>0
        phiStar=phiFit;
        sigStar=sigFit;
        fitAccepted=true;
    end
end

hole=B.hole;
c=hole.center(:).';
R=hole.r;

nMat=[cos(phiStar),sin(phiStar)];
tHat=[-sin(phiStar),cos(phiStar)];
xStar=c+R*nMat;

lambdaIni=C.sig_c/sigStar;
sig0=1.0;
if isfield(C,'load')&&isfield(C.load,'sig0')&&~isempty(C.load.sig0)
    sig0=C.load.sig0;
end

% Retain a conventional near-tie list for diagnostics only.
tol=getf(C.stage1,'max_tie_rel_tol',1e-8)*max(1,abs(sigDisc));
allMaxIdx=find(abs(sigPos-sigDisc)<=tol);

I=struct();
I.idx_star=idxDisc;
I.idx_discrete=idxDisc;
I.phi_discrete=phiDisc;
if logical(getf(B,'periodic',true))
    I.phi_star=mod(phiStar,2*pi);
else
    I.phi_star=phiStar;
end
I.x_star=xStar;
I.n_mat_star=nMat;
I.n_hole_star=-nMat;
I.t_hat_star=tHat;

I.sig_tt_discrete_unit=sigDisc;
I.sig_tt_unit=sigStar;
I.sig_tt_pos_unit=sigStar;

I.lambda_ini=lambdaIni;
I.sig_applied_ini=lambdaIni*sig0;

I.all_max_idx=allMaxIdx;
if isfield(C.stage1,'angular_fit_halfwidth_factor') && ...
        ~isempty(C.stage1.angular_fit_halfwidth_factor)
    I.selection_rule='boundary_extrapolation_plus_mesh_scaled_quadratic_peak';
else
    I.selection_rule='boundary_extrapolation_plus_local_quadratic_peak';
end
I.angular_fit=fitInfo;
I.angular_fit.accepted=fitAccepted;
I.boundary_stress_method=getf(B,'method','unknown');
end


function [phiStar,sigStar,F]=periodic_quadratic_peak(phi,sig,idx0,nFit)
n=numel(phi);
m=(nFit-1)/2;
off=(-m:m).';

idx=1+mod((idx0-1)+off,n);
phi0=phi(idx0);

dphi=atan2(sin(phi(idx)-phi0),cos(phi(idx)-phi0));
y=sig(idx);

[dphi,ord]=sort(dphi);
y=y(ord);
idx=idx(ord);

p=polyfit(dphi,y,2);

F=struct();
F.accepted=false;
F.reason='';
F.indices=idx;
F.dphi=dphi;
F.values=y;
F.coefficients=p;
F.curvature=p(1);
F.vertex_offset=NaN;
F.rmse=sqrt(mean((y-polyval(p,dphi)).^2));

phiStar=phi0;
sigStar=sig(idx0);

if ~(isfinite(p(1))&&isfinite(p(2))&&isfinite(p(3)))
    F.reason='nonfinite_coefficients';
    return;
end
if p(1)>=0
    F.reason='nonnegative_curvature';
    return;
end

dv=-p(2)/(2*p(1));
span=max(abs(dphi));
F.vertex_offset=dv;

if ~isfinite(dv) || abs(dv)>span
    F.reason='vertex_outside_fit_window';
    return;
end

sf=polyval(p,dv);
if ~isfinite(sf) || sf<=0
    F.reason='invalid_fitted_peak';
    return;
end

phiStar=mod(phi0+dv,2*pi);
sigStar=sf;
F.accepted=true;
F.reason='accepted';
end


function [phiStar,sigStar,F]=window_quadratic_peak(phi,sig,idx0,halfWidth,isPeriodic)
phi0=phi(idx0);

if isPeriodic
    dphi=atan2(sin(phi-phi0),cos(phi-phi0));
else
    dphi=phi-phi0;
end

keep=abs(dphi)<=halfWidth+100*eps;
idx=find(keep);
x=dphi(keep);
y=sig(keep);

[x,ord]=sort(x);
y=y(ord);
idx=idx(ord);

F=struct();
F.accepted=false;
F.reason='';
F.indices=idx;
F.dphi=x;
F.values=y;
F.coefficients=nan(1,3);
F.curvature=NaN;
F.vertex_offset=NaN;
F.rmse=NaN;

phiStar=phi0;
sigStar=sig(idx0);

if numel(x)<5
    F.reason='too_few_points';
    return;
end

p=polyfit(x,y,2);
F.coefficients=p;
F.curvature=p(1);
F.rmse=sqrt(mean((y-polyval(p,x)).^2));

if any(~isfinite(p))
    F.reason='nonfinite_coefficients';
    return;
end
if p(1)>=0
    F.reason='nonnegative_curvature';
    return;
end

dv=-p(2)/(2*p(1));
F.vertex_offset=dv;
if ~isfinite(dv) || abs(dv)>halfWidth
    F.reason='vertex_outside_fit_window';
    return;
end

sf=polyval(p,dv);
if ~isfinite(sf)||sf<=0
    F.reason='invalid_fitted_peak';
    return;
end

phiStar=phi0+dv;
if isPeriodic
    phiStar=mod(phiStar,2*pi);
end
sigStar=sf;
F.accepted=true;
F.reason='accepted';
end


function must(S,f)
if ~isfield(S,f)||isempty(S.(f))
    error('find_hole_initiation_point_v2:MissingField', ...
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
