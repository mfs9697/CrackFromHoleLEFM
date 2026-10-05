function Out=main_step11_production_symmetry_parity(varargin)
%MAIN_STEP11_PRODUCTION_SYMMETRY_PARITY
% Symmetry-enforced production benchmark for the centered-hole problem.
%
% The current Stage-I numerical initiation detector selects phi=359 deg on
% the nominally symmetric centered-hole problem. This driver removes that
% Stage-I sampling/mesh bias by imposing the exact rightmost initiation point
% phi=0, n_mat=[1,0], t_hat=[0,1], then runs matched +/-theta Stage-II
% production meshes.
%
% The primary EDI parity checks are:
%   KII(0) ~ 0,
%   KII(-theta) ~ -KII(+theta),
%   KI(-theta) ~ KI(+theta).
%
% The historical decomposed-J branch is treated only as a magnitude check:
%   |KII(-theta)| ~ |KII(+theta)|.
%
% Defaults:
%   thetaDegList    = -3:0.5:3
%   radiusFractions = [0.50 0.65 0.80]

ip=inputParser;
addParameter(ip,'thetaDegList',-3:0.5:3, ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x)));
addParameter(ip,'radiusFractions',[0.50 0.65 0.80], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0)&&all(x<1));
addParameter(ip,'OldNtheta',240,@(x)isnumeric(x)&&isscalar(x)&&x>=40);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 11: PRODUCTION SYMMETRY / PARITY BENCHMARK\n');
fprintf('============================================================\n');

C=cfg_hole_initiation();

% Numerical Stage-I solution retained only to quantify its bias.
G=geom_hole_only(C);
S1=solve_hole_only(C,G,'lambda',1.0);
B=sample_hole_boundary_stress(C,G,S1);
Inum=find_hole_initiation_point(C,B);

% Exact symmetry-enforced rightmost initiation point.
ctr=C.hole.center(:).';
R=C.hole.r;
I=Inum;
I.phi_star=0;
I.x_star=ctr+[R,0];
I.n_mat_star=[1,0];
I.n_hole_star=[-1,0];
I.t_hat_star=[0,1];
I.selection_rule='exact_rightmost_symmetry_point';

% Compare sampled Stage-I stress at phi=0 with numerical selected maximum.
[~,i0]=min(abs(local_wrap_to_pi(B.phi)));
sig0=B.sig_tt_eff(i0);
sigmax=max(B.sig_tt_eff);
fprintf('Numerical Stage-I selected phi = %.6f deg\n',rad2deg(Inum.phi_star));
fprintf('Exact symmetry point phi       = 0 deg\n');
fprintf('sig_tt(phi=0)                  = %.10e\n',sig0);
fprintf('sig_tt numerical max           = %.10e\n',sigmax);
fprintf('relative Stage-I stress bias   = %.6e\n',(sigmax-sig0)/max(abs(sigmax),eps));

thetaDeg=O.thetaDegList(:);
if ~any(abs(thetaDeg)<1e-12)
    thetaDeg=sort([thetaDeg;0]);
end

% Require paired +/- angles for parity reporting.
tol=1e-12;
for k=1:numel(thetaDeg)
    if ~any(abs(thetaDeg+thetaDeg(k))<tol)
        error('step11:UnpairedTheta', ...
            'thetaDegList must contain +/- pairs and zero.');
    end
end

theta=deg2rad(thetaDeg);
rf=O.radiusFractions(:).';
nT=numel(thetaDeg);
nR=numel(rf);

KIold=nan(nT,nR);
KIIold=nan(nT,nR);
KIedi=nan(nT,nR);
KIIedi=nan(nT,nR);
rInner=nan(nT,nR);
hTip=nan(nT,1);
valid=false(nT,1);

for it=1:nT
    fprintf('\n--- theta %+7.3f deg ---\n',thetaDeg(it));
    try
        [G2,~,~,Mc]=build_stage2_cracked_mesh_for_theta( ...
            C,I,theta(it), ...
            'PlotGeom',false,'PlotMesh',false,'PlotCollapsed',false);
        S2=solve_cracked_LEFM(C,Mc);
        Llast=norm(Mc.crack.Pmid(end,:)-Mc.crack.Pmid(end-1,:));

        for ir=1:nR
            rOut=rf(ir)*Llast;
            Rcmp=compute_SIF_for_stage2_compare( ...
                C,G2,Mc,S2, ...
                'OldRadius',rOut, ...
                'EDIOuterRadius',rOut, ...
                'OldNtheta',O.OldNtheta);

            KIold(it,ir)=Rcmp.KI_old;
            KIIold(it,ir)=Rcmp.KII_old;
            KIedi(it,ir)=Rcmp.KI_EDI;
            KIIedi(it,ir)=Rcmp.KII_EDI;
            rInner(it,ir)=Rcmp.domain_EDI.r_inner;
            hTip(it)=Rcmp.tipMeshScale.median;

            fprintf(['  r/a0=%.2f | old |KII|~%.6e | ', ...
                     'EDI KI=% .6e KII=% .6e | rin/a0=%.3f\n'], ...
                rf(ir),abs(Rcmp.KII_old),Rcmp.KI_EDI,Rcmp.KII_EDI, ...
                Rcmp.domain_EDI.r_inner/Llast);
        end
        valid(it)=true;
    catch ME
        fprintf('  FAILED: %s\n',ME.message);
    end
end

% Parity diagnostics per radius.
P=nan(nR,10);
for ir=1:nR
    [k0,zeroFound]=find_zero_index(thetaDeg);
    if zeroFound
        KII0=KIIedi(k0,ir);
    else
        KII0=NaN;
    end

    oddErr=[];
    evenKIErr=[];
    oldMagEvenErr=[];

    pos=find(thetaDeg>tol & valid);
    for jp=1:numel(pos)
        ip=pos(jp);
        im=find(abs(thetaDeg+thetaDeg(ip))<tol,1);
        if isempty(im) || ~valid(im), continue; end

        oddErr(end+1,1)=abs(KIIedi(ip,ir)+KIIedi(im,ir)); %#ok<AGROW>
        evenKIErr(end+1,1)=abs(KIedi(ip,ir)-KIedi(im,ir)); %#ok<AGROW>
        oldMagEvenErr(end+1,1)=abs(abs(KIIold(ip,ir))-abs(KIIold(im,ir))); %#ok<AGROW>
    end

    scaleKII=max(abs(KIIedi(:,ir)),[],'omitnan');
    scaleKI=max(abs(KIedi(:,ir)),[],'omitnan');
    scaleOld=max(abs(KIIold(:,ir)),[],'omitnan');

    Root=local_signed_root(thetaDeg,KIIedi(:,ir),valid);

    P(ir,:)=[ ...
        rf(ir),KII0, ...
        local_max_or_nan(oddErr), ...
        local_max_or_nan(oddErr)/max(scaleKII,eps), ...
        local_max_or_nan(evenKIErr), ...
        local_max_or_nan(evenKIErr)/max(scaleKI,eps), ...
        local_max_or_nan(oldMagEvenErr), ...
        local_max_or_nan(oldMagEvenErr)/max(scaleOld,eps), ...
        Root.rootDeg,Root.nSignChanges];
end

Tpar=array2table(P,'VariableNames',{ ...
    'r_over_a0','EDI_KII_at_zero', ...
    'EDI_KII_odd_absmax','EDI_KII_odd_relmax', ...
    'EDI_KI_even_absmax','EDI_KI_even_relmax', ...
    'old_absKII_even_absmax','old_absKII_even_relmax', ...
    'EDI_root_theta_deg','EDI_n_sign_changes'});

fprintf('\nSYMMETRY / PARITY SUMMARY\n');
disp(Tpar);

% Long table.
nRows=nT*nR;
A=nan(nRows,9); row=0;
for it=1:nT
    for ir=1:nR
        row=row+1;
        A(row,:)=[thetaDeg(it),rf(ir), ...
            KIold(it,ir),KIIold(it,ir), ...
            KIedi(it,ir),KIIedi(it,ir), ...
            rInner(it,ir),hTip(it),double(valid(it))];
    end
end
T=array2table(A,'VariableNames',{ ...
    'thetaDeg','r_over_a0','KI_old','KII_old_legacy', ...
    'KI_EDI','KII_EDI','EDI_r_inner','h_tip_median','valid'});

if logical(O.Plot)
    local_plot_edi(thetaDeg,KIIedi,rf);
    local_plot_parity(thetaDeg,KIedi,KIIedi,rf);
end

Out=struct();
Out.config=C;
Out.stage1Numeric=Inum;
Out.stage1Symmetric=I;
Out.stage1Boundary=B;
Out.thetaDeg=thetaDeg;
Out.radiusFractions=rf;
Out.KI_old=KIold;
Out.KII_old=KIIold;
Out.KI_EDI=KIedi;
Out.KII_EDI=KIIedi;
Out.parity=Tpar;
Out.table=T;
Out.valid=valid;
Out.settings=O;

fprintf('\nSTEP 11 completed.\n');
fprintf(['Acceptance target: EDI should show one root near 0 deg, small KII(0), ', ...
    'approximately odd KII(theta), and approximately even KI(theta).\n']);
end


function a=local_wrap_to_pi(a)
a=mod(a+pi,2*pi)-pi;
end


function [idx,ok]=find_zero_index(x)
[amin,idx]=min(abs(x));
ok=amin<1e-12;
end


function y=local_max_or_nan(x)
if isempty(x), y=NaN; else, y=max(x); end
end


function R=local_signed_root(x,y,valid)
good=valid(:)&isfinite(x(:))&isfinite(y(:));
x=x(good); y=y(good);
R=struct('rootDeg',NaN,'nSignChanges',0);
if numel(x)<2, return; end

[amin,iz]=min(abs(y));
scale=max(max(abs(y)),eps);
if amin<=100*eps(scale)
    R.rootDeg=x(iz);
    R.nSignChanges=1;
    return;
end

idx=find(y(1:end-1).*y(2:end)<0);
R.nSignChanges=numel(idx);
if isempty(idx), return; end

score=abs(y(idx))+abs(y(idx+1));
[~,j]=min(score);
i=idx(j);
R.rootDeg=x(i)-y(i)*(x(i+1)-x(i))/(y(i+1)-y(i));
end


function local_plot_edi(thetaDeg,KII,rf)
figure('Name','Step 11: EDI symmetry benchmark','Color','w');
clf; hold on; box on; grid on;
for j=1:numel(rf)
    plot(thetaDeg,KII(:,j),'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('r/a_0=%.2f',rf(j)));
end
yline(0,'k--','HandleVisibility','off');
xline(0,'k:','HandleVisibility','off');
xlabel('\theta [deg]');
ylabel('K_{II}^{EDI}');
legend('Location','best');
title('Symmetry-enforced production benchmark');
end


function local_plot_parity(thetaDeg,KI,KII,rf)
[~,ir]=min(abs(rf-0.65));
figure('Name','Step 11: parity check','Color','w');
clf; hold on; box on; grid on;
plot(thetaDeg,KII(:,ir),'-o','LineWidth',1.1,'DisplayName','K_{II}^{EDI}(\theta)');
plot(thetaDeg,-flipud(KII(:,ir)),'--s','LineWidth',1.1, ...
    'DisplayName','-K_{II}^{EDI}(-\theta)');
xlabel('\theta [deg]');
ylabel('K_{II}');
legend('Location','best');
title(sprintf('Odd-parity check, r/a_0=%.2f',rf(ir)));

figure('Name','Step 11: KI evenness','Color','w');
clf; hold on; box on; grid on;
plot(thetaDeg,KI(:,ir),'-o','LineWidth',1.1,'DisplayName','K_I^{EDI}(\theta)');
plot(thetaDeg,flipud(KI(:,ir)),'--s','LineWidth',1.1, ...
    'DisplayName','K_I^{EDI}(-\theta)');
xlabel('\theta [deg]');
ylabel('K_I');
legend('Location','best');
title(sprintf('Even-parity check, r/a_0=%.2f',rf(ir)));
end
