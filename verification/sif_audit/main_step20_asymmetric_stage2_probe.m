function Out=main_step20_asymmetric_stage2_probe(varargin)
%MAIN_STEP20_ASYMMETRIC_STAGE2_PROBE
% First signed-EDI Stage-II direction probe for the asymmetric full domain.
%
% This is intentionally a small computational step.  It uses only the
% Npoly=240 asymmetric benchmark, evaluates a modest set of trial crack
% directions, and checks whether KII(theta) has a clean sign change.
%
% theta is measured relative to the Stage-I material-side normal at the
% fitted initiation point.
%
% Defaults:
%   Npoly = 240
%   ThetaDeg = [-4 -2 -1 0 1 2 4]
%   ROuterOverA0 = [0.50 0.65 0.80]
%
% No root-refinement mesh is generated in this step.  Root estimates are
% obtained only by interpolation of the signed-EDI probe values.

ip=inputParser;
addParameter(ip,'Npoly',240,@(x)isnumeric(x)&&isscalar(x)&&x>=32);
addParameter(ip,'ThetaDeg',[-4 -2 -1 0 1 2 4], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x)));
addParameter(ip,'ROuterOverA0',[0.50 0.65 0.80], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 20: ASYMMETRIC FULL-DOMAIN STAGE-II PROBE\n');
fprintf('============================================================\n');

C=cfg_asymmetric_full_domain();
C.hole.npoly=round(O.Npoly);
C.holes={C.hole};

hArc=2*pi*C.hole.r/C.hole.npoly;
C.mesh1.hmin=hArc;
C.mesh1.hhole=hArc;
C.mesh1.hmax=20*hArc;
C.mesh2.hmax=C.mesh1.hmax;
C.mesh2.hhole=C.mesh1.hmin;
C.mesh2.hcrack=C.mesh1.hmin;
C.solver.verbose=0;
C.plot.show_mesh1=false;

R1=run_stage1_hole_initiation(C);
I=R1.I;

fprintf('\nStage I:\n');
fprintf('  phi_* = %+12.8f deg\n',rad2deg(local_wrap(I.phi_star)));
fprintf('  x_*   = [%.10e, %.10e]\n',I.x_star(1),I.x_star(2));
fprintf('  n_mat = [%.10e, %.10e]\n',I.n_mat_star(1),I.n_mat_star(2));
fprintf('  sigma_tt,max(unit) = %.10e\n',I.sig_tt_pos_unit);
fprintf('  lambda_ini = %.10e\n',I.lambda_ini);

thetaDeg=O.ThetaDeg(:).';
thetaRad=deg2rad(thetaDeg);
rRat=O.ROuterOverA0(:).';

nT=numel(thetaDeg);
nR=numel(rRat);

KI=nan(nT,nR);
KII=nan(nT,nR);
rIn=nan(nT,nR);
hTip=nan(nT,1);
Cases=cell(nT,1);

for it=1:nT
    th=thetaRad(it);
    fprintf('\n--- theta=%+7.3f deg ---\n',thetaDeg(it));

    doPlotCase=logical(O.Plot)&&abs(thetaDeg(it))<1e-12;

    [G2,D,M,Mc]=build_stage2_cracked_mesh_for_theta( ...
        C,I,th, ...
        'PlotGeom',false, ...
        'PlotMesh',false, ...
        'PlotCollapsed',doPlotCase);

    S2=solve_cracked_LEFM(C,Mc,'lambda',1.0);

    H=local_tip_mesh_scale(S2.mesh,Mc.crack.Pmid(end,:));
    hTip(it)=H.median;

    mat=S2.mat;
    if ~isfield(mat,'Dmat')
        mat.Dmat=mat.D;
    end

    V=Mc.crack.Pmid;

    for ir=1:nR
        rOut=rRat(ir)*C.a0;
        rin=max(0.10*rOut,2.0*H.median);

        if rin>=rOut
            error('step20:BadEDIDomain', ...
                'Invalid EDI domain at theta=%g deg and r/a0=%g.', ...
                thetaDeg(it),rRat(ir));
        end

        dom=struct('r_inner',rin,'r_outer',rOut);

        [KI(it,ir),KII(it,ir)]=SIF_LEFM_interaction_EDI( ...
            S2.mesh,S2.U,V,mat,dom, ...
            'UsePlaneStrain',mat.ps==1, ...
            'Verbose',false, ...
            'WeightFunction','fe_nodal');

        rIn(it,ir)=rin;

        fprintf('  r_o/a0=%.2f | KI=%.8e | KII=%+.8e | KII/KI=%+.6e\n', ...
            rRat(ir),KI(it,ir),KII(it,ir),KII(it,ir)/KI(it,ir));
    end

    Cases{it}=struct('G2',G2,'D',D,'M',M,'Mc',Mc,'S2',S2,'H',H);
end

% Interpolated signed roots for each EDI outer radius.
rootRows=nan(nR,9);

for ir=1:nR
    q=KII(:,ir);
    k=KI(:,ir);

    rootBracket=local_sign_change_root(thetaRad,q);
    rootLinear=local_linear_root(thetaRad,q);

    if isfinite(rootBracket)
        KIroot=interp1(thetaRad,k,rootBracket,'linear');
        slope=local_local_slope(thetaRad,q,rootBracket);
    else
        KIroot=NaN;
        slope=NaN;
    end

    [~,i0]=min(abs(thetaDeg));
    q0=q(i0);
    k0=k(i0);

    rootRows(ir,:)=[ ...
        rRat(ir),q0,k0,q0/k0, ...
        rad2deg(rootBracket),rad2deg(rootLinear), ...
        KIroot,slope,max(abs(q./k))];
end

Roots=array2table(rootRows,'VariableNames',{ ...
    'r_outer_over_a0','KII_theta0','KI_theta0','KII_over_KI_theta0', ...
    'theta_root_bracket_deg','theta_root_linear_deg', ...
    'KI_at_bracket_root','dKII_dtheta_near_root','max_abs_KII_over_KI_probe'});

fprintf('\nSIGNED-EDI ROOT ESTIMATES\n');
disp(Roots);

% Cross-domain root spread.
validRoots=Roots.theta_root_bracket_deg(isfinite(Roots.theta_root_bracket_deg));
if isempty(validRoots)
    rootSpread=NaN;
    rootMean=NaN;
    rootMaxDev=NaN;
else
    rootSpread=max(validRoots)-min(validRoots);
    rootMean=mean(validRoots);
    rootMaxDev=max(abs(validRoots-rootMean));
end

fprintf('\nPROBE SUMMARY\n');
fprintf('  bracket roots available = %d / %d EDI domains\n',numel(validRoots),nR);
fprintf('  mean bracket root       = %+12.8f deg\n',rootMean);
fprintf('  root domain spread      = %.8g deg\n',rootSpread);
fprintf('  max deviation from mean = %.8g deg\n',rootMaxDev);

if logical(O.Plot)
    local_plot_signed(thetaDeg,KI,KII,rRat);
    local_plot_absolute(thetaDeg,KII,rRat);
end

Out=struct();
Out.C=C;
Out.Stage1=R1;
Out.thetaDeg=thetaDeg;
Out.rOuterOverA0=rRat;
Out.KI=KI;
Out.KII=KII;
Out.rInner=rIn;
Out.hTip=hTip;
Out.Cases=Cases;
Out.roots=Roots;
Out.rootMeanDeg=rootMean;
Out.rootSpreadDeg=rootSpread;
Out.rootMaxDeviationDeg=rootMaxDev;
Out.settings=O;

fprintf('\nSTEP 20 completed.\n');
fprintf(['Gate: all EDI domains should show the same KII sign trend and a ', ...
    'consistent sign-change root within the probe interval. Only then should ', ...
    'we generate a dedicated root-refinement mesh.\n']);
end


function root=local_sign_change_root(theta,y)
root=NaN;
best=inf;
for k=1:numel(theta)-1
    if y(k)==0
        cand=theta(k);
    elseif y(k)*y(k+1)<0
        cand=theta(k)-y(k)*(theta(k+1)-theta(k))/(y(k+1)-y(k));
    else
        continue;
    end

    if abs(cand)<best
        best=abs(cand);
        root=cand;
    end
end
end


function root=local_linear_root(theta,y)
p=polyfit(theta,y,1);
if abs(p(1))<=eps
    root=NaN;
else
    root=-p(2)/p(1);
end
end


function slope=local_local_slope(theta,y,root)
slope=NaN;
for k=1:numel(theta)-1
    lo=min(theta(k),theta(k+1));
    hi=max(theta(k),theta(k+1));
    if root>=lo && root<=hi && theta(k+1)~=theta(k)
        slope=(y(k+1)-y(k))/(theta(k+1)-theta(k));
        return;
    end
end
end


function H=local_tip_mesh_scale(mesh,tip)
X=mesh.coord3;
T=mesh.connect3;
d=sqrt(sum((X-tip).^2,2));
dmin=min(d);
tol=max(1e-12,1e-8*max(1,max(abs(X(:)))));
tipNodes=find(d<=dmin+tol);

hit=any(ismember(T,tipNodes),2);
Te=T(hit,:);

L=[];
for k=1:size(Te,1)
    P=X(Te(k,:),:);
    L=[L,norm(P(2,:)-P(1,:)),norm(P(3,:)-P(2,:)),norm(P(1,:)-P(3,:))]; %#ok<AGROW>
end

L=L(isfinite(L)&L>tol);
if isempty(L)
    error('step20:NoTipEdges','Could not determine the crack-tip mesh scale.');
end

H=struct('min',min(L),'median',median(L),'max',max(L), ...
    'nTipElements',size(Te,1),'tipNodeDistance',dmin);
end


function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end


function local_plot_signed(thetaDeg,KI,KII,rRat)
figure('Name','Step 20: asymmetric signed Stage-II response','Color','w');
clf; hold on; box on; grid on

for ir=1:numel(rRat)
    plot(thetaDeg,KII(:,ir)./KI(:,ir),'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('r_o/a_0=%.2f',rRat(ir)));
end

xline(0,'k--','HandleVisibility','off');
yline(0,'k:','HandleVisibility','off');
xlabel('\theta relative to Stage-I normal [deg]');
ylabel('K_{II}/K_I');
title('Asymmetric Stage-II signed interaction EDI');
legend('Location','best');
end


function local_plot_absolute(thetaDeg,KII,rRat)
figure('Name','Step 20: asymmetric KII','Color','w');
clf; hold on; box on; grid on

for ir=1:numel(rRat)
    plot(thetaDeg,KII(:,ir),'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('r_o/a_0=%.2f',rRat(ir)));
end

xline(0,'k--','HandleVisibility','off');
yline(0,'k:','HandleVisibility','off');
xlabel('\theta relative to Stage-I normal [deg]');
ylabel('K_{II}');
title('Asymmetric Stage-II local-symmetry crossing');
legend('Location','best');
end
