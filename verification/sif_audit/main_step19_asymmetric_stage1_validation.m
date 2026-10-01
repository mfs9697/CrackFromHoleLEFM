function Out=main_step19_asymmetric_stage1_validation(varargin)
%MAIN_STEP19_ASYMMETRIC_STAGE1_VALIDATION
% Validate redesigned Stage I on a genuinely asymmetric full-domain hole.
%
% The benchmark uses cfg_asymmetric_full_domain() and checks:
%   1) a unique dominant boundary-stress maximum exists;
%   2) the fitted initiation angle is stable under mesh refinement;
%   3) the mesh-scaled angular window c*h_hole/R preserves the displaced
%      physical peak rather than smoothing it back toward a symmetry axis;
%   4) radial boundary extrapolation remains well conditioned.
%
% Defaults:
%   NpolyList = [240 480]
%   WindowFactors = [1 1.5 2 3]
%   Plot = true

ip=inputParser;
addParameter(ip,'NpolyList',[240 480], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=32));
addParameter(ip,'WindowFactors',[1 1.5 2 3], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 19: ASYMMETRIC FULL-DOMAIN STAGE-I VALIDATION\n');
fprintf('============================================================\n');

Nlist=round(O.NpolyList(:));
fac=O.WindowFactors(:).';
nL=numel(Nlist);
nF=numel(fac);

Runs=cell(nL,1);
Fitted=cell(nL,1);
Peaks=cell(nL,1);

rows=nan(nL*nF,16);
rr=0;

for il=1:nL
    C=cfg_asymmetric_full_domain();

    Np=Nlist(il);
    C.hole.npoly=Np;
    C.holes={C.hole};

    hArc=2*pi*C.hole.r/Np;
    C.mesh1.hmin=hArc;
    C.mesh1.hhole=hArc;
    C.mesh1.hmax=20*hArc;

    C.mesh2.hmax=C.mesh1.hmax;
    C.mesh2.hhole=C.mesh1.hmin;
    C.mesh2.hcrack=C.mesh1.hmin;

    C.solver.verbose=0;
    C.plot.show_mesh1=false;

    % Solve once per mesh.  The default c=3 finder result is stored in R,
    % but all window-factor comparisons below reuse the same B field.
    R=run_stage1_hole_initiation(C);
    Runs{il}=R;

    B=R.B;
    sig=max(B.sig_tt_eff(:),0);
    phi=B.phi(:);

    P=local_rank_peaks(B,3.0,6);
    Peaks{il}=P;

    fprintf('\n--- Npoly=%d | h/R=%.8g rad = %.6f deg ---\n', ...
        Np,hArc/C.hole.r,rad2deg(hArc/C.hole.r));
    fprintf('  hole center = [%.6f, %.6f] m, R=%.6f m\n', ...
        C.hole.center(1),C.hole.center(2),C.hole.r);
    fprintf('  all radial query fractions inside FEM = %s\n', ...
        mat2str(B.offset.fraction_inside_actual_FEM_domain,5));
    fprintf('  max relative radial-fit RMSE(sigma_tt) = %.6e\n', ...
        max(B.fit.sig_tt.rmse)/max(max(abs(B.sig_tt_eff)),eps));

    if height(P)>=1
        fprintf('  dominant local peak: phi=%+.8f deg, sigma=%.10e\n', ...
            P.phi_fit_deg(1),P.sig_fit(1));
    end
    if height(P)>=2
        fprintf(['  second local peak  : phi=%+.8f deg, sigma=%.10e | ', ...
            'top-second gap=% .6e\n'], ...
            P.phi_fit_deg(2),P.sig_fit(2),P.top_minus_this_rel(2));
        fprintf('  angular separation of top two = %.6f deg\n', ...
            P.separation_from_top_deg(2));
    end

    Fits=cell(nF,1);

    for jf=1:nF
        Ct=C;
        Ct.stage1.angular_fit_halfwidth_factor=fac(jf);

        I=find_hole_initiation_point_v2(Ct,B);
        Fits{jf}=I;

        rr=rr+1;

        fitRMSE=I.angular_fit.rmse/max(abs(I.sig_tt_pos_unit),eps);
        halfWidth=I.angular_fit.halfwidth_rad;

        rightScale=max(abs(B.sig_tt_eff));
        maxNn=max(abs(B.sig_nn))/max(rightScale,eps);
        maxNt=max(abs(B.sig_nt))/max(rightScale,eps);

        if height(P)>=2
            gap=P.top_minus_this_rel(2);
            sep=P.separation_from_top_deg(2);
        else
            gap=NaN;
            sep=NaN;
        end

        rows(rr,:)=[ ...
            Np,hArc/C.hole.r,fac(jf),rad2deg(halfWidth), ...
            rad2deg(local_wrap(I.phi_discrete)), ...
            rad2deg(local_wrap(I.phi_star)), ...
            I.sig_tt_discrete_unit,I.sig_tt_pos_unit,I.lambda_ini, ...
            double(I.angular_fit.accepted),fitRMSE, ...
            gap,sep, ...
            max(B.fit.sig_tt.rmse)/max(max(abs(B.sig_tt_eff)),eps), ...
            maxNn,maxNt];

        fprintf(['  c=%3.1f | halfwidth=%6.3f deg | discrete=%+9.4f deg | ', ...
            'fit=%+10.6f deg | sigma=%.8e | relRMSE=% .3e\n'], ...
            fac(jf),rad2deg(halfWidth), ...
            rad2deg(local_wrap(I.phi_discrete)), ...
            rad2deg(local_wrap(I.phi_star)), ...
            I.sig_tt_pos_unit,fitRMSE);
    end

    Fitted{il}=Fits;
end

T=array2table(rows,'VariableNames',{ ...
    'Npoly','h_over_R','window_factor','halfwidth_deg', ...
    'phi_discrete_deg','phi_fit_deg', ...
    'sig_discrete','sig_fit','lambda_ini', ...
    'fit_accepted','fit_rmse_rel', ...
    'top_second_gap_rel','top_two_separation_deg', ...
    'radial_fit_rmse_relmax','sig_nn_relmax','sig_nt_relmax'});

fprintf('\nASYMMETRIC STAGE-I WINDOW TABLE\n');
disp(T);

% Cross-mesh comparison by window factor, using wrapped angular differences.
Srows=nan(nF,9);
for jf=1:nF
    Q=T(abs(T.window_factor-fac(jf))<1e-12,:);
    p1=deg2rad(Q.phi_fit_deg(1));
    p2=deg2rad(Q.phi_fit_deg(end));
    dphi=rad2deg(local_wrap(p2-p1));

    ds=(Q.sig_fit(end)-Q.sig_fit(1))/max(abs(Q.sig_fit(end)),eps);
    dl=(Q.lambda_ini(end)-Q.lambda_ini(1))/max(abs(Q.lambda_ini(end)),eps);

    Srows(jf,:)=[ ...
        fac(jf),Q.phi_fit_deg(1),Q.phi_fit_deg(end),dphi, ...
        abs(dphi),ds,dl, ...
        max(Q.fit_rmse_rel),min(Q.top_second_gap_rel)];
end

S=array2table(Srows,'VariableNames',{ ...
    'window_factor','phi_fit_first_deg','phi_fit_last_deg', ...
    'phi_change_signed_deg','phi_change_abs_deg', ...
    'sig_fit_rel_change','lambda_rel_change', ...
    'max_fit_rmse_rel','min_top_second_gap_rel'});

fprintf('\nCROSS-MESH ASYMMETRIC SUMMARY\n');
disp(S);

if logical(O.Plot)
    local_plot_boundary(Runs,Nlist);
    local_plot_window(T,Nlist);
    local_plot_peak_zoom(Runs,Fitted,Nlist,fac);
    local_plot_peak_ranking(Peaks,Nlist);
end

Out=struct();
Out.settings=O;
Out.runs=Runs;
Out.fits=Fitted;
Out.peaks=Peaks;
Out.table=T;
Out.summary=S;

fprintf('\nSTEP 19 completed.\n');
fprintf(['Gate: the asymmetric case should have a clearly preferred peak, ', ...
    'its fitted angle should remain nonzero and stable from Npoly=240 to 480, ', ...
    'and c=3 should agree with the narrower mesh-scaled windows without ', ...
    'erasing the physical displacement of the maximum.\n']);
end


function P=local_rank_peaks(B,factor,nKeep)
phi=B.phi(:);
sig=max(B.sig_tt_eff(:),0);
n=numel(sig);

prev=sig(1+mod((0:n-1)-1,n));
next=sig(1+mod((0:n-1)+1,n));
idx=find(sig>=prev(:) & sig>=next(:) & sig>0);

halfWidth=factor*B.offset.hhole/B.hole.r;

vals=[];
for k=1:numel(idx)
    F=local_refine(phi,sig,idx(k),halfWidth);
    if F.accepted
        ph=F.phiStar;
        sf=F.sigStar;
    else
        ph=phi(idx(k));
        sf=sig(idx(k));
    end
    vals(end+1,:)=[idx(k),rad2deg(local_wrap(phi(idx(k)))), ...
        rad2deg(local_wrap(ph)),sf,double(F.accepted),F.rmse]; %#ok<AGROW>
end

if isempty(vals)
    P=array2table(zeros(0,8),'VariableNames',{ ...
        'idx','phi_discrete_deg','phi_fit_deg','sig_fit','fit_accepted', ...
        'fit_rmse','top_minus_this_rel','separation_from_top_deg'});
    return;
end

[~,ord]=sort(vals(:,4),'descend');
vals=vals(ord,:);
vals=vals(1:min(nKeep,size(vals,1)),:);

topSig=vals(1,4);
topPhi=deg2rad(vals(1,3));

gap=(topSig-vals(:,4))/max(abs(topSig),eps);
sep=abs(rad2deg(local_wrap(deg2rad(vals(:,3))-topPhi)));

P=array2table([vals,gap,sep],'VariableNames',{ ...
    'idx','phi_discrete_deg','phi_fit_deg','sig_fit','fit_accepted', ...
    'fit_rmse','top_minus_this_rel','separation_from_top_deg'});
end


function F=local_refine(phi,sig,idx0,halfWidth)
phi0=phi(idx0);
d=atan2(sin(phi-phi0),cos(phi-phi0));
keep=abs(d)<=halfWidth+100*eps;

x=d(keep);
y=sig(keep);
[x,ord]=sort(x);
y=y(ord);

F=struct('accepted',false,'phiStar',phi0,'sigStar',sig(idx0), ...
    'rmse',NaN,'curvature',NaN,'vertexOffset',NaN);

if numel(x)<5
    return;
end

p=polyfit(x,y,2);
F.rmse=sqrt(mean((y-polyval(p,x)).^2));
F.curvature=p(1);

if any(~isfinite(p)) || p(1)>=0
    return;
end

dv=-p(2)/(2*p(1));
F.vertexOffset=dv;

if ~isfinite(dv) || abs(dv)>halfWidth
    return;
end

sf=polyval(p,dv);
if ~isfinite(sf) || sf<=0
    return;
end

F.accepted=true;
F.phiStar=mod(phi0+dv,2*pi);
F.sigStar=sf;
end


function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end


function local_plot_boundary(Runs,Nlist)
figure('Name','Step 19: asymmetric boundary stress','Color','w');
clf; hold on; box on; grid on

for il=1:numel(Nlist)
    B=Runs{il}.B;
    x=rad2deg(local_wrap(B.phi));
    [x,ord]=sort(x);
    y=B.sig_tt_eff(ord);
    plot(x,y,'LineWidth',1.2,'DisplayName',sprintf('Npoly=%d',Nlist(il)));
end

xlabel('\phi [deg]');
ylabel('\sigma_{tt}^{boundary}');
title('Asymmetric full-domain boundary-limit stress');
legend('Location','best');
end


function local_plot_window(T,Nlist)
figure('Name','Step 19: asymmetric fitted angle','Color','w');
clf; hold on; box on; grid on

for il=1:numel(Nlist)
    Q=T(T.Npoly==Nlist(il),:);
    plot(Q.window_factor,Q.phi_fit_deg,'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('Npoly=%d',Nlist(il)));
end

xlabel('angular half-width factor c');
ylabel('\phi_* [deg]');
title('Asymmetric initiation angle versus mesh-scaled window');
legend('Location','best');
end


function local_plot_peak_zoom(Runs,Fitted,Nlist,fac)
figure('Name','Step 19: asymmetric peak zoom','Color','w');
clf; tiledlayout(numel(Nlist),1,'Padding','compact','TileSpacing','compact');

for il=1:numel(Nlist)
    B=Runs{il}.B;
    I=Fitted{il}{end};
    center=local_wrap(I.phi_star);

    d=local_wrap(B.phi-center);
    keep=abs(d)<=deg2rad(12);
    x=rad2deg(d(keep));
    y=B.sig_tt_eff(keep);
    [x,ord]=sort(x); y=y(ord);

    nexttile; hold on; box on; grid on
    plot(x,y,'-o','MarkerSize',3,'LineWidth',1.0,'DisplayName','boundary stress');

    for jf=1:numel(fac)
        J=Fitted{il}{jf};
        dx=rad2deg(local_wrap(J.phi_star-center));
        xline(dx,':','DisplayName',sprintf('c=%.1f',fac(jf)));
    end

    xline(0,'k--','HandleVisibility','off');
    xlabel('\phi-\phi_*^{c=3} [deg]');
    ylabel('\sigma_{tt}^{boundary}');
    title(sprintf('Dominant asymmetric peak, Npoly=%d',Nlist(il)));
    legend('Location','best');
end
end


function local_plot_peak_ranking(Peaks,Nlist)
figure('Name','Step 19: local peak ranking','Color','w');
clf; hold on; box on; grid on

for il=1:numel(Nlist)
    P=Peaks{il};
    if isempty(P), continue; end
    k=(1:height(P)).';
    plot(k,P.sig_fit/P.sig_fit(1),'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('Npoly=%d',Nlist(il)));
end

xlabel('ranked local maximum');
ylabel('\sigma_{peak}/\sigma_{peak,1}');
title('Asymmetric boundary-stress peak hierarchy');
legend('Location','best');
end
