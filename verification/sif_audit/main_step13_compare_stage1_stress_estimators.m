function Out=main_step13_compare_stage1_stress_estimators(varargin)
%MAIN_STEP13_COMPARE_STAGE1_STRESS_ESTIMATORS
% Compare two material-side Stage-I stress estimators on identical FEM solves:
%
%   RECOVERED-T6:
%     StressExt GP stresses -> extrapolated/area-averaged nodal stresses ->
%     topology-respecting T6 interpolation.
%
%   DIRECT-U:
%     eps(xq)=B(xq)Ue -> sig(xq)=D eps(xq) directly in the containing T6.
%
% Both use exactly the same material-side query points and the same
% ShiftFraction*h_hole offset. No production code is changed.
%
% Default refinement family matches Step 12:
%   Npoly = [120 180 240 360 480], Nphi=1440.
%
% Diagnostics:
%   - sigma_tt(0), right-peak angle/value/bias,
%   - local +/-phi symmetry,
%   - estimator-to-estimator difference,
%   - sigma_nn and sigma_nt residuals near the right-hole peak,
%   - simple angular roughness indicator of sigma_tt.

ip=inputParser;
addParameter(ip,'NpolyList',[120 180 240 360 480], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=32));
addParameter(ip,'Nphi',1440, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>=360);
addParameter(ip,'WindowDeg',12, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<45);
addParameter(ip,'ShiftFraction',0.25, ...
    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 13: STAGE-I STRESS ESTIMATOR COMPARISON\n');
fprintf('============================================================\n');
fprintf('ShiftFraction = %.4g of h_hole\n',O.ShiftFraction);

Nlist=round(O.NpolyList(:));
nL=numel(Nlist);

rows=nan(nL,30);
Recovered=cell(nL,1);
Direct=cell(nL,1);
Solutions=cell(nL,1);

for il=1:nL
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

    fprintf('\n--- refinement %d/%d : Npoly=%d, h_arc=%.8g ---\n', ...
        il,nL,Np,hArc);

    G=geom_hole_only(C);
    S=solve_hole_only(C,G,'lambda',1.0);

    Br=sample_hole_boundary_stress_audit(C,G,S, ...
        'QuerySide','material','Interpolator','t6', ...
        'ShiftFraction',O.ShiftFraction);

    Bd=sample_hole_boundary_stress_direct_U_audit(C,G,S, ...
        'ShiftFraction',O.ShiftFraction);

    Mr=local_metrics(Br,O.WindowDeg);
    Md=local_metrics(Bd,O.WindowDeg);
    Cmp=local_compare(Br,Bd,O.WindowDeg);

    Recovered{il}=Br;
    Direct{il}=Bd;
    Solutions{il}=struct('C',C,'G',G,'S',S);

    rows(il,:)=[ ...
        Np,hArc,hArc/C.hole.r, ...
        size(G.p,1),size(G.t,1),size(S.mesh.coord,1),numel(S.U), ...
        Br.eps_shift,Br.eps_shift/C.hole.r, ...
        Mr.rightPeakPhiDeg,Mr.sigAtZero,Mr.rightPeakSig,Mr.rightPeakBiasRel, ...
        Mr.localSymRelMax,Mr.sigNN_relmax,Mr.sigNT_relmax,Mr.roughnessRel, ...
        Md.rightPeakPhiDeg,Md.sigAtZero,Md.rightPeakSig,Md.rightPeakBiasRel, ...
        Md.localSymRelMax,Md.sigNN_relmax,Md.sigNT_relmax,Md.roughnessRel, ...
        Cmp.sigTT_at0_relDiff,Cmp.sigTT_window_relLinf,Cmp.sigTT_window_relL2, ...
        Cmp.peakSig_relDiff,Cmp.peakAngle_absDiffDeg];

    fprintf(['  recovered-T6: phi_R=%+6.3f deg, sig0=%.8f, ', ...
             'bias=% .3e, sym=% .3e, nn=% .3e, nt=% .3e\n'], ...
        Mr.rightPeakPhiDeg,Mr.sigAtZero,Mr.rightPeakBiasRel, ...
        Mr.localSymRelMax,Mr.sigNN_relmax,Mr.sigNT_relmax);

    fprintf(['  direct-U    : phi_R=%+6.3f deg, sig0=%.8f, ', ...
             'bias=% .3e, sym=% .3e, nn=% .3e, nt=% .3e\n'], ...
        Md.rightPeakPhiDeg,Md.sigAtZero,Md.rightPeakBiasRel, ...
        Md.localSymRelMax,Md.sigNN_relmax,Md.sigNT_relmax);

    fprintf(['  estimator diff: sig0=% .3e rel, window Linf=% .3e, ', ...
             'window L2=% .3e\n'], ...
        Cmp.sigTT_at0_relDiff,Cmp.sigTT_window_relLinf, ...
        Cmp.sigTT_window_relL2);
end

T=array2table(rows,'VariableNames',{ ...
    'Npoly','h_arc','h_arc_over_R', ...
    'T3_nodes','T3_elements','T6_nodes','DOF', ...
    'eps_shift','eps_shift_over_R', ...
    'rec_peak_phi_deg','rec_sig_tt_phi0','rec_peak_sig_tt', ...
    'rec_peak_bias_rel','rec_sym_relmax','rec_sig_nn_relmax', ...
    'rec_sig_nt_relmax','rec_roughness_rel', ...
    'dir_peak_phi_deg','dir_sig_tt_phi0','dir_peak_sig_tt', ...
    'dir_peak_bias_rel','dir_sym_relmax','dir_sig_nn_relmax', ...
    'dir_sig_nt_relmax','dir_roughness_rel', ...
    'sig_tt_phi0_relDiff','sig_tt_window_relLinf','sig_tt_window_relL2', ...
    'peak_sig_relDiff','peak_angle_absDiff_deg'});

fprintf('\nSTRESS-ESTIMATOR COMPARISON TABLE\n');
disp(T);

fprintf('\nCOMPACT CONVERGENCE VIEW\n');
disp(T(:,{ ...
    'Npoly','h_arc_over_R', ...
    'rec_peak_phi_deg','dir_peak_phi_deg', ...
    'rec_sig_tt_phi0','dir_sig_tt_phi0','sig_tt_phi0_relDiff', ...
    'sig_tt_window_relLinf','sig_tt_window_relL2', ...
    'rec_sym_relmax','dir_sym_relmax', ...
    'rec_sig_nn_relmax','dir_sig_nn_relmax', ...
    'rec_sig_nt_relmax','dir_sig_nt_relmax', ...
    'rec_roughness_rel','dir_roughness_rel'}));

if logical(O.Plot)
    local_plot_sig0(T);
    local_plot_estimator_difference(T);
    local_plot_traction_residuals(T);
    local_plot_curves(Recovered,Direct,Nlist,O.WindowDeg);
end

Out=struct();
Out.table=T;
Out.recovered=Recovered;
Out.direct=Direct;
Out.solutions=Solutions;
Out.NpolyList=Nlist;
Out.settings=O;

fprintf('\nSTEP 13 completed.\n');
fprintf(['Primary interpretation: convergence of recovered-T6 and direct-U ', ...
    'toward the same material-side stress field validates stress recovery; ', ...
    'systematic disagreement would implicate StressExt/recovery.\n']);
end


function M=local_metrics(B,windowDeg)
phi=local_wrap(B.phi);
sig=B.sig_tt_eff(:);
nn=B.sig_nn(:);
nt=B.sig_nt(:);
win=deg2rad(windowDeg);

right=abs(phi)<=win;
[sp,ipLocal]=max(sig(right));
ids=find(right);
ip=ids(ipLocal);

[~,i0]=min(abs(phi));
sig0=sig(i0);

pos=find(phi>0 & phi<=win);
diffs=nan(numel(pos),1);
for j=1:numel(pos)
    [~,im]=min(abs(phi+phi(pos(j))));
    diffs(j)=abs(sig(pos(j))-sig(im));
end

scale=max(abs(sig(right)));
symAbs=max(diffs);
symRel=symAbs/max(scale,eps);

% Traction-free boundary diagnostics evaluated at the offset material-side
% query ring. They need not vanish exactly, but should decrease as the query
% approaches the boundary with refinement.
nnRel=max(abs(nn(right)))/max(scale,eps);
ntRel=max(abs(nt(right)))/max(scale,eps);

% Simple roughness measure based on second differences. For a smooth angular
% curve this decays rapidly; element-crossing noise tends to increase it.
ids=find(right);
[phis,ord]=sort(phi(ids));
vals=sig(ids(ord)); %#ok<ASGLU>
if numel(vals)>=3
    sec=diff(vals,2);
    rough=max(abs(sec))/max(scale,eps);
else
    rough=NaN;
end

M=struct();
M.rightPeakPhiDeg=rad2deg(phi(ip));
M.sigAtZero=sig0;
M.rightPeakSig=sp;
M.rightPeakBiasRel=(sp-sig0)/max(abs(sp),eps);
M.localSymRelMax=symRel;
M.sigNN_relmax=nnRel;
M.sigNT_relmax=ntRel;
M.roughnessRel=rough;
end


function Cmp=local_compare(A,B,windowDeg)
phi=local_wrap(A.phi);
if numel(phi)~=numel(B.phi) || max(abs(local_wrap(A.phi)-local_wrap(B.phi)))>1e-12
    error('step13:SamplingMismatch','Estimator sampling grids do not match.');
end

win=abs(phi)<=deg2rad(windowDeg);
sa=A.sig_tt_eff(:);
sb=B.sig_tt_eff(:);

[~,i0]=min(abs(phi));
scale=max([abs(sa(win));abs(sb(win));eps]);

d=sa(win)-sb(win);

Cmp=struct();
Cmp.sigTT_at0_relDiff=(sa(i0)-sb(i0))/max([abs(sa(i0)),abs(sb(i0)),eps]);
Cmp.sigTT_window_relLinf=max(abs(d))/scale;
Cmp.sigTT_window_relL2=norm(d)/max(norm([sa(win);sb(win)])/sqrt(2),eps);

spa=max(sa(win));
spb=max(sb(win));
Cmp.peakSig_relDiff=(spa-spb)/max([abs(spa),abs(spb),eps]);

ids=find(win);
[~,ia]=max(sa(win)); ia=ids(ia);
[~,ib]=max(sb(win)); ib=ids(ib);
Cmp.peakAngle_absDiffDeg=abs(rad2deg(phi(ia)-phi(ib)));
end


function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end


function local_plot_sig0(T)
figure('Name','Step 13: sigma_tt(phi=0) convergence','Color','w');
clf; hold on; box on; grid on;
plot(T.h_arc_over_R,T.rec_sig_tt_phi0,'-o','LineWidth',1.1, ...
    'DisplayName','recovered-T6');
plot(T.h_arc_over_R,T.dir_sig_tt_phi0,'-s','LineWidth',1.1, ...
    'DisplayName','direct-U');
set(gca,'XDir','reverse');
xlabel('h_{arc}/R (refinement -> right)');
ylabel('\sigma_{tt}(0)');
legend('Location','best');
title('Stage-I stress estimator convergence at symmetry point');
end


function local_plot_estimator_difference(T)
figure('Name','Step 13: estimator difference','Color','w');
clf; hold on; box on; grid on;
semilogy(T.h_arc_over_R,abs(T.sig_tt_phi0_relDiff),'-o','LineWidth',1.1, ...
    'DisplayName','|difference at \phi=0|');
semilogy(T.h_arc_over_R,T.sig_tt_window_relLinf,'-s','LineWidth',1.1, ...
    'DisplayName','window L_\infty difference');
semilogy(T.h_arc_over_R,T.sig_tt_window_relL2,'-^','LineWidth',1.1, ...
    'DisplayName','window relative L_2 difference');
set(gca,'XDir','reverse');
xlabel('h_{arc}/R (refinement -> right)');
ylabel('relative difference');
legend('Location','best');
title('Recovered-T6 versus direct-U');
end


function local_plot_traction_residuals(T)
figure('Name','Step 13: material-side traction residuals','Color','w');
clf; hold on; box on; grid on;
semilogy(T.h_arc_over_R,T.rec_sig_nn_relmax,'-o','LineWidth',1.1, ...
    'DisplayName','recovered: |\sigma_{nn}|');
semilogy(T.h_arc_over_R,T.dir_sig_nn_relmax,'--o','LineWidth',1.1, ...
    'DisplayName','direct: |\sigma_{nn}|');
semilogy(T.h_arc_over_R,T.rec_sig_nt_relmax,'-s','LineWidth',1.1, ...
    'DisplayName','recovered: |\sigma_{nt}|');
semilogy(T.h_arc_over_R,T.dir_sig_nt_relmax,'--s','LineWidth',1.1, ...
    'DisplayName','direct: |\sigma_{nt}|');
set(gca,'XDir','reverse');
xlabel('h_{arc}/R (refinement -> right)');
ylabel('max residual / \sigma_{tt} scale');
legend('Location','best');
title('Approach to traction-free hole boundary');
end


function local_plot_curves(R,D,Nlist,windowDeg)
figure('Name','Step 13: finest stress curves','Color','w');
clf; hold on; box on; grid on;
il=numel(Nlist);
for k=1:2
    if k==1
        B=R{il}; name='recovered-T6';
    else
        B=D{il}; name='direct-U';
    end
    phi=rad2deg(local_wrap(B.phi));
    keep=abs(phi)<=windowDeg;
    [x,ord]=sort(phi(keep));
    y=B.sig_tt_eff(keep); y=y(ord);
    plot(x,y,'LineWidth',1.2,'DisplayName',name);
end
xline(0,'k--','HandleVisibility','off');
xlabel('\phi [deg]');
ylabel('\sigma_{tt}');
legend('Location','best');
title(sprintf('Finest mesh comparison, Npoly=%d',Nlist(il)));
end
