function Out=main_step16_stage1_angular_window_audit(varargin)
%MAIN_STEP16_STAGE1_ANGULAR_WINDOW_AUDIT
% Audit the local angular regression window used to refine the Stage-I peak.
%
% Motivation:
% Step 15 showed that the redesigned radial stress estimator is excellent,
% but a 5-point quadratic fit locks onto a tiny mesh-scale ripple and returns
% about -0.49 deg instead of the exact centered-hole symmetry direction.
%
% This audit keeps the new boundary-limit stress field fixed and varies the
% angular half-width of a local quadratic regression according to
%
%     halfWidth = factor * (h_hole/R).
%
% Thus the smoothing/regression scale is tied to the boundary mesh and shrinks
% consistently under refinement.
%
% Defaults:
%   NpolyList = [240 480]
%   HalfWidthFactors = [0.5 1 1.5 2 3]
%   Nphi = 1440
%
% No production angular-fit settings are changed by this audit.

ip=inputParser;
addParameter(ip,'NpolyList',[240 480], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=32));
addParameter(ip,'HalfWidthFactors',[0.5 1 1.5 2 3], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
addParameter(ip,'Nphi',1440,@(x)isnumeric(x)&&isscalar(x)&&x>=360);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 16: STAGE-I ANGULAR WINDOW SENSITIVITY\n');
fprintf('============================================================\n');

Nlist=round(O.NpolyList(:));
fac=O.HalfWidthFactors(:).';

rows=[];
Curves=cell(numel(Nlist),1);

for il=1:numel(Nlist)
    C=cfg_hole_initiation();
    C.solver.verbose=0;
    C.plot.show_mesh1=false;
    C.stage1.nphi=O.Nphi;

    Np=Nlist(il);
    C.hole.npoly=Np;
    C.holes={C.hole};
    hArc=2*pi*C.hole.r/Np;
    C.mesh1.hmin=hArc;
    C.mesh1.hhole=hArc;
    C.mesh1.hmax=20*hArc;

    G=geom_hole_only(C);
    S=solve_hole_only(C,G,'lambda',1.0);
    B=sample_hole_boundary_stress_v2(C,G,S);

    Curves{il}=B;

    sig=max(B.sig_tt_eff(:),0);
    phi=B.phi(:);
    [sigDisc,idxDisc]=max(sig);
    phiDisc=phi(idxDisc);

    meshAngle=C.mesh1.hhole/C.hole.r;

    fprintf('\n--- Npoly=%d | h/R=%.8g rad = %.6f deg ---\n', ...
        Np,meshAngle,rad2deg(meshAngle));
    fprintf('  discrete maximum: phi=%+.6f deg, sigma=%.10e\n', ...
        rad2deg(local_wrap(phiDisc)),sigDisc);

    for jf=1:numel(fac)
        f=fac(jf);
        hw=f*meshAngle;

        F=local_quadratic_window(phi,sig,phiDisc,hw);

        [~,i0]=min(abs(local_wrap(phi)));
        sig0=sig(i0);

        if F.accepted
            phiFit=local_wrap(F.phiStar);
            sigFit=F.sigStar;
            [phiRef,localOffset]=local_nearest_centered_hole_peak(phiFit);
            dist0=abs(rad2deg(localOffset));

            [~,iRef]=min(abs(local_wrap(phi-phiRef)));
            sigRef=sig(iRef);
            peakExcess=(sigFit-sigRef)/max(abs(sigFit),eps);
        else
            phiFit=NaN;
            sigFit=NaN;
            dist0=NaN;
            peakExcess=NaN;
            localOffset=NaN;
        end

        rows(end+1,:)=[ ... %#ok<AGROW>
            il,Np,meshAngle,f,hw,rad2deg(hw), ...
            numel(F.indices),double(F.accepted), ...
            rad2deg(local_wrap(phiDisc)),sigDisc, ...
            rad2deg(phiFit),sigFit,dist0,rad2deg(localOffset), ...
            peakExcess,F.rmse,F.rmse/max(abs(sigFit),eps), ...
            F.curvature,F.vertexOffsetDeg];
        
        fprintf(['  factor=%4.1f | halfwidth=%6.3f deg | n=%2d | ', ...
                 'phi_fit=%+9.5f deg | excess0=% .3e | relRMSE=% .3e\n'], ...
            f,rad2deg(hw),numel(F.indices),rad2deg(phiFit), ...
            peakExcess,F.rmse/max(abs(sigFit),eps));
    end
end

T=array2table(rows,'VariableNames',{ ...
    'level','Npoly','mesh_angle_rad','halfwidth_factor', ...
    'halfwidth_rad','halfwidth_deg','nfit','accepted', ...
    'phi_discrete_deg','sig_discrete', ...
    'phi_fit_deg','sig_fit','distance_to_nearest_exact_peak_deg', ...
    'local_peak_offset_deg','peak_excess_over_exact_peak_rel', ...
    'fit_rmse','fit_rmse_rel', ...
    'curvature','vertex_offset_from_discrete_deg'});

fprintf('\nANGULAR WINDOW AUDIT TABLE\n');
disp(T);

% Cross-mesh comparison by factor.
Srows=nan(numel(fac),7);
for jf=1:numel(fac)
    Q=T(abs(T.halfwidth_factor-fac(jf))<1e-12,:);
    if height(Q)>=2
        p240=Q.local_peak_offset_deg(1);
        p480=Q.local_peak_offset_deg(end);
        spread=max(Q.local_peak_offset_deg)-min(Q.local_peak_offset_deg);
        d0=max(abs(Q.local_peak_offset_deg));
        maxEx=max(abs(Q.peak_excess_over_exact_peak_rel));
        maxRMSE=max(Q.fit_rmse_rel);
        minN=min(Q.nfit);
    else
        p240=NaN;p480=NaN;spread=NaN;d0=NaN;maxEx=NaN;maxRMSE=NaN;minN=NaN;
    end
    Srows(jf,:)=[fac(jf),p240,p480,spread,d0,maxEx,maxRMSE];
end

S=array2table(Srows,'VariableNames',{ ...
    'halfwidth_factor','local_offset_first_deg','local_offset_last_deg', ...
    'cross_mesh_spread_deg','max_abs_local_offset_deg', ...
    'max_peak_excess_rel','max_fit_rmse_rel'});

fprintf('\nCROSS-MESH WINDOW SUMMARY\n');
disp(S);

if logical(O.Plot)
    local_plot_phi(T,fac);
    local_plot_fit_quality(T,fac);
    local_plot_curves(Curves,Nlist,T,fac);
end

Out=struct();
Out.table=T;
Out.summary=S;
Out.curves=Curves;
Out.settings=O;

fprintf('\nSTEP 16 completed.\n');
fprintf(['Preferred window should drive the centered-hole fitted angle toward ', ...
    'the nearest exact symmetry-equivalent peak (0 or 180 deg), remain stable ', ...
    'under refinement, and retain a ', ...
    'small regression residual without using a window so broad that real ', ...
    'asymmetry would be washed out.\n']);
end


function F=local_quadratic_window(phi,sig,phiCenter,halfWidth)
d=atan2(sin(phi-phiCenter),cos(phi-phiCenter));
keep=abs(d)<=halfWidth+100*eps;

x=d(keep);
y=sig(keep);
idx=find(keep);

[x,ord]=sort(x);
y=y(ord);
idx=idx(ord);

F=struct();
F.accepted=false;
F.indices=idx;
F.rmse=NaN;
F.curvature=NaN;
F.vertexOffsetDeg=NaN;
F.phiStar=NaN;
F.sigStar=NaN;

if numel(x)<5
    return;
end

p=polyfit(x,y,2);
yf=polyval(p,x);
F.rmse=sqrt(mean((y-yf).^2));
F.curvature=p(1);

if any(~isfinite(p)) || p(1)>=0
    return;
end

dv=-p(2)/(2*p(1));
F.vertexOffsetDeg=rad2deg(dv);

if ~isfinite(dv) || abs(dv)>halfWidth
    return;
end

sf=polyval(p,dv);
if ~isfinite(sf)||sf<=0
    return;
end

F.accepted=true;
F.phiStar=mod(phiCenter+dv,2*pi);
F.sigStar=sf;
end


function [phiRef,offset]=local_nearest_centered_hole_peak(phi)
% The centered-hole benchmark has two physically equivalent exact maxima:
% phi = 0 and phi = pi. Return the nearest representative and signed offset.
d0=local_wrap(phi);
dpi=local_wrap(phi-pi);
if abs(d0)<=abs(dpi)
    phiRef=0;
    offset=d0;
else
    phiRef=pi;
    offset=dpi;
end
end


function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end


function local_plot_phi(T,fac)
figure('Name','Step 16: fitted peak versus window','Color','w');
clf; hold on; box on; grid on;
N=unique(T.Npoly);
for k=1:numel(N)
    Q=T(T.Npoly==N(k),:);
    plot(Q.halfwidth_factor,Q.phi_fit_deg,'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('Npoly=%d',N(k)));
end
yline(0,'k--','HandleVisibility','off');
xlabel('half-width / (h_{hole}/R)');
ylabel('\phi_* [deg]');
legend('Location','best');
title('Angular regression window sensitivity');
end


function local_plot_fit_quality(T,fac) %#ok<INUSD>
figure('Name','Step 16: fit residual','Color','w');
clf; hold on; box on; grid on;
N=unique(T.Npoly);
for k=1:numel(N)
    Q=T(T.Npoly==N(k),:);
    semilogy(Q.halfwidth_factor,Q.fit_rmse_rel,'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('Npoly=%d',N(k)));
end
xlabel('half-width / (h_{hole}/R)');
ylabel('relative quadratic-fit RMSE');
legend('Location','best');
title('Angular fit quality versus window');
end


function local_plot_curves(Curves,Nlist,T,fac)
for il=1:numel(Nlist)
    B=Curves{il};
    phi=rad2deg(local_wrap(B.phi));
    keep=abs(phi)<=8;
    [x,ord]=sort(phi(keep));
    y=B.sig_tt_eff(keep); y=y(ord);

    figure('Name',sprintf('Step 16: Npoly=%d',Nlist(il)),'Color','w');
    clf; hold on; box on; grid on;
    plot(x,y,'-o','LineWidth',1.0,'MarkerSize',3,'DisplayName','boundary stress');

    Q=T(T.Npoly==Nlist(il),:);
    for jf=1:numel(fac)
        qq=Q(abs(Q.halfwidth_factor-fac(jf))<1e-12,:);
        if ~isempty(qq)&&qq.accepted
            xline(qq.phi_fit_deg(1),':', ...
                'DisplayName',sprintf('f=%.1f',fac(jf)));
        end
    end

    xline(0,'k--','HandleVisibility','off');
    xlabel('\phi [deg]');
    ylabel('\sigma_{tt}^{boundary}');
    legend('Location','best');
    title(sprintf('Angular-window audit, Npoly=%d',Nlist(il)));
end
end
