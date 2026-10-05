function Out=main_step14_stage1_offset_sensitivity(varargin)
%MAIN_STEP14_STAGE1_OFFSET_SENSITIVITY
% Small material-side query-offset experiment for Stage-I hole initiation.
%
% Purpose
% -------
% Quantify sensitivity to the material-side sampling distance
%
%     eps = ShiftFraction * h_hole
%
% after Steps 12--13 established that:
%   - cavity-side sampling is invalid;
%   - recovered-T6 and direct-U stresses converge toward the same field.
%
% This driver solves each selected mesh only once, then evaluates BOTH
% estimators for several ShiftFraction values.
%
% Defaults:
%   NpolyList     = [240 480]
%   ShiftFractions= [0.05 0.10 0.25 0.50]
%   Nphi          = 1440
%   WindowDeg     = 12
%
% Diagnostics
% -----------
% For each mesh/offset/estimator:
%   sigma_tt(0), right-window discrete peak angle/value/bias,
%   +/-phi symmetry defect,
%   max |sigma_nn| and |sigma_nt| relative to the local sigma_tt scale,
%   angular roughness.
%
% For each mesh/offset pair:
%   recovered-T6 vs direct-U differences at phi=0 and over the window.
%
% No production code is changed.

ip=inputParser;
addParameter(ip,'NpolyList',[240 480], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=32));
addParameter(ip,'ShiftFractions',[0.05 0.10 0.25 0.50], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
addParameter(ip,'Nphi',1440, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>=360);
addParameter(ip,'WindowDeg',12, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<45);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 14: STAGE-I MATERIAL-SIDE OFFSET SENSITIVITY\n');
fprintf('============================================================\n');
fprintf('Npoly list       = %s\n',mat2str(O.NpolyList));
fprintf('Shift fractions  = %s\n',mat2str(O.ShiftFractions));

Nlist=round(O.NpolyList(:));
sf=O.ShiftFractions(:).';

nL=numel(Nlist);
nS=numel(sf);

rows=nan(nL*nS,33);
Recovered=cell(nL,nS);
Direct=cell(nL,nS);
Solutions=cell(nL,1);
rr=0;

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

    fprintf('\n--- mesh %d/%d : Npoly=%d, h_arc=%.8g ---\n', ...
        il,nL,Np,hArc);

    G=geom_hole_only(C);
    S=solve_hole_only(C,G,'lambda',1.0);
    Solutions{il}=struct('C',C,'G',G,'S',S);

    fprintf('T3 nodes/elements = %d / %d | T6 nodes/DOF = %d / %d\n', ...
        size(G.p,1),size(G.t,1),size(S.mesh.coord,1),numel(S.U));

    for is=1:nS
        f=sf(is);

        Br=sample_hole_boundary_stress_audit(C,G,S, ...
            'QuerySide','material','Interpolator','t6', ...
            'ShiftFraction',f);

        Bd=sample_hole_boundary_stress_direct_U_audit(C,G,S, ...
            'ShiftFraction',f);

        Mr=local_metrics(Br,O.WindowDeg);
        Md=local_metrics(Bd,O.WindowDeg);
        Cmp=local_compare(Br,Bd,O.WindowDeg);

        Recovered{il,is}=Br;
        Direct{il,is}=Bd;

        rr=rr+1;
        rows(rr,:)=[ ...
            il,Np,hArc,hArc/C.hole.r, ...
            f,Br.eps_shift,Br.eps_shift/C.hole.r, ...
            Br.fractionInsideActualFEMDomain, ...
            Mr.rightPeakPhiDeg,Mr.sigAtZero,Mr.rightPeakSig, ...
            Mr.rightPeakBiasRel,Mr.localSymRelMax, ...
            Mr.sigNN_relmax,Mr.sigNT_relmax,Mr.roughnessRel, ...
            Md.rightPeakPhiDeg,Md.sigAtZero,Md.rightPeakSig, ...
            Md.rightPeakBiasRel,Md.localSymRelMax, ...
            Md.sigNN_relmax,Md.sigNT_relmax,Md.roughnessRel, ...
            Cmp.sigTT_at0_relDiff,Cmp.sigTT_window_relLinf, ...
            Cmp.sigTT_window_relL2,Cmp.peakSig_relDiff, ...
            Cmp.peakAngle_absDiffDeg, ...
            size(G.p,1),size(G.t,1),size(S.mesh.coord,1),numel(S.U)];

        fprintf(['  eps/h=%4.2f | rec: phi=%+5.2f sig0=%.8f ', ...
                 'nn=% .3e nt=% .3e | dir: phi=%+5.2f sig0=%.8f ', ...
                 'nn=% .3e nt=% .3e | d0=% .3e\n'], ...
            f,Mr.rightPeakPhiDeg,Mr.sigAtZero,Mr.sigNN_relmax,Mr.sigNT_relmax, ...
            Md.rightPeakPhiDeg,Md.sigAtZero,Md.sigNN_relmax,Md.sigNT_relmax, ...
            Cmp.sigTT_at0_relDiff);
    end
end

T=array2table(rows,'VariableNames',{ ...
    'level','Npoly','h_arc','h_arc_over_R', ...
    'shift_fraction','eps_shift','eps_shift_over_R', ...
    'fraction_query_in_actual_FEM_domain', ...
    'rec_peak_phi_deg','rec_sig_tt_phi0','rec_peak_sig_tt', ...
    'rec_peak_bias_rel','rec_sym_relmax', ...
    'rec_sig_nn_relmax','rec_sig_nt_relmax','rec_roughness_rel', ...
    'dir_peak_phi_deg','dir_sig_tt_phi0','dir_peak_sig_tt', ...
    'dir_peak_bias_rel','dir_sym_relmax', ...
    'dir_sig_nn_relmax','dir_sig_nt_relmax','dir_roughness_rel', ...
    'sig_tt_phi0_relDiff','sig_tt_window_relLinf', ...
    'sig_tt_window_relL2','peak_sig_relDiff','peak_angle_absDiff_deg', ...
    'T3_nodes','T3_elements','T6_nodes','DOF'});

fprintf('\nOFFSET-SENSITIVITY TABLE\n');
disp(T);

fprintf('\nCOMPACT OFFSET VIEW\n');
disp(T(:,{ ...
    'Npoly','shift_fraction','eps_shift_over_R', ...
    'rec_peak_phi_deg','dir_peak_phi_deg', ...
    'rec_sig_tt_phi0','dir_sig_tt_phi0', ...
    'rec_peak_bias_rel','dir_peak_bias_rel', ...
    'rec_sym_relmax','dir_sym_relmax', ...
    'rec_sig_nn_relmax','dir_sig_nn_relmax', ...
    'rec_sig_nt_relmax','dir_sig_nt_relmax', ...
    'sig_tt_phi0_relDiff','sig_tt_window_relLinf'}));

% Within-mesh offset ranges, useful for production decision.
Rtab=nan(nL,14);
for il=1:nL
    Q=T(T.level==il,:);
    Rtab(il,:)=[ ...
        Nlist(il), ...
        max(Q.rec_sig_tt_phi0)-min(Q.rec_sig_tt_phi0), ...
        local_rel_range(Q.rec_sig_tt_phi0), ...
        max(Q.dir_sig_tt_phi0)-min(Q.dir_sig_tt_phi0), ...
        local_rel_range(Q.dir_sig_tt_phi0), ...
        max(Q.rec_peak_bias_rel),min(Q.rec_peak_bias_rel), ...
        max(Q.dir_peak_bias_rel),min(Q.dir_peak_bias_rel), ...
        max(Q.rec_sig_nn_relmax),min(Q.rec_sig_nn_relmax), ...
        max(Q.dir_sig_nn_relmax),min(Q.dir_sig_nn_relmax), ...
        max(Q.sig_tt_window_relLinf)];
end

Toffset=array2table(Rtab,'VariableNames',{ ...
    'Npoly', ...
    'rec_sig0_abs_range','rec_sig0_rel_range', ...
    'dir_sig0_abs_range','dir_sig0_rel_range', ...
    'rec_peak_bias_max','rec_peak_bias_min', ...
    'dir_peak_bias_max','dir_peak_bias_min', ...
    'rec_nn_max','rec_nn_min','dir_nn_max','dir_nn_min', ...
    'estimator_window_Linf_max'});

fprintf('\nWITHIN-MESH OFFSET RANGE SUMMARY\n');
disp(Toffset);

if logical(O.Plot)
    local_plot_sig0(T,Nlist);
    local_plot_traction(T,Nlist);
    local_plot_bias(T,Nlist);
    local_plot_finest_curves(Recovered,Direct,Nlist,sf,O.WindowDeg);
end

Out=struct();
Out.table=T;
Out.offsetSummary=Toffset;
Out.recovered=Recovered;
Out.direct=Direct;
Out.solutions=Solutions;
Out.NpolyList=Nlist;
Out.shiftFractions=sf;
Out.settings=O;

fprintf('\nSTEP 14 completed.\n');
fprintf(['Decision quantities: offset sensitivity of sigma_tt(0), peak bias and ', ...
    'traction residuals; plus recovered/direct agreement as eps approaches ', ...
    'the exact-circle boundary.\n']);
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
if isempty(diffs)
    symRel=NaN;
else
    symRel=max(diffs)/max(scale,eps);
end

nnRel=max(abs(nn(right)))/max(scale,eps);
ntRel=max(abs(nt(right)))/max(scale,eps);

ids=find(right);
[~,ord]=sort(phi(ids));
vals=sig(ids(ord));
if numel(vals)>=3
    rough=max(abs(diff(vals,2)))/max(scale,eps);
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


function r=local_rel_range(x)
r=(max(x)-min(x))/max(max(abs(x)),eps);
end


function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end


function local_plot_sig0(T,Nlist)
figure('Name','Step 14: sigma_tt(0) vs query offset','Color','w');
clf; hold on; box on; grid on;
for il=1:numel(Nlist)
    Q=T(T.Npoly==Nlist(il),:);
    plot(Q.shift_fraction,Q.rec_sig_tt_phi0,'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('rec N=%d',Nlist(il)));
    plot(Q.shift_fraction,Q.dir_sig_tt_phi0,'--s','LineWidth',1.1, ...
        'DisplayName',sprintf('dir N=%d',Nlist(il)));
end
xlabel('\epsilon/h_{hole}');
ylabel('\sigma_{tt}(0)');
legend('Location','best');
title('Material-side offset sensitivity');
end


function local_plot_traction(T,Nlist)
figure('Name','Step 14: traction residuals vs offset','Color','w');
clf; tiledlayout(numel(Nlist),1,'Padding','compact','TileSpacing','compact');
for il=1:numel(Nlist)
    nexttile; hold on; box on; grid on;
    Q=T(T.Npoly==Nlist(il),:);
    semilogy(Q.shift_fraction,Q.rec_sig_nn_relmax,'-o','DisplayName','rec |\sigma_{nn}|');
    semilogy(Q.shift_fraction,Q.rec_sig_nt_relmax,'-s','DisplayName','rec |\sigma_{nt}|');
    semilogy(Q.shift_fraction,Q.dir_sig_nn_relmax,'--o','DisplayName','dir |\sigma_{nn}|');
    semilogy(Q.shift_fraction,Q.dir_sig_nt_relmax,'--s','DisplayName','dir |\sigma_{nt}|');
    xlabel('\epsilon/h_{hole}');
    ylabel('relative residual');
    title(sprintf('Npoly=%d',Nlist(il)));
    legend('Location','best');
end
end


function local_plot_bias(T,Nlist)
figure('Name','Step 14: peak bias vs offset','Color','w');
clf; hold on; box on; grid on;
for il=1:numel(Nlist)
    Q=T(T.Npoly==Nlist(il),:);
    semilogy(Q.shift_fraction,Q.rec_peak_bias_rel,'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('rec N=%d',Nlist(il)));
    semilogy(Q.shift_fraction,Q.dir_peak_bias_rel,'--s','LineWidth',1.1, ...
        'DisplayName',sprintf('dir N=%d',Nlist(il)));
end
xlabel('\epsilon/h_{hole}');
ylabel('(\sigma_{max}-\sigma(0))/\sigma_{max}');
legend('Location','best');
title('Discrete peak excess versus query offset');
end


function local_plot_finest_curves(R,D,Nlist,sf,windowDeg)
il=numel(Nlist);
figure('Name','Step 14: finest-mesh offset family','Color','w');
clf; hold on; box on; grid on;
for is=1:numel(sf)
    B=R{il,is};
    phi=rad2deg(local_wrap(B.phi));
    keep=abs(phi)<=windowDeg;
    [x,ord]=sort(phi(keep));
    y=B.sig_tt_eff(keep); y=y(ord);
    plot(x,y,'LineWidth',1.0, ...
        'DisplayName',sprintf('rec eps/h=%.2f',sf(is)));
end
xline(0,'k--','HandleVisibility','off');
xlabel('\phi [deg]'); ylabel('\sigma_{tt}');
legend('Location','best');
title(sprintf('Recovered-T6 offset family, Npoly=%d',Nlist(il)));

figure('Name','Step 14: finest direct-U offset family','Color','w');
clf; hold on; box on; grid on;
for is=1:numel(sf)
    B=D{il,is};
    phi=rad2deg(local_wrap(B.phi));
    keep=abs(phi)<=windowDeg;
    [x,ord]=sort(phi(keep));
    y=B.sig_tt_eff(keep); y=y(ord);
    plot(x,y,'LineWidth',1.0, ...
        'DisplayName',sprintf('dir eps/h=%.2f',sf(is)));
end
xline(0,'k--','HandleVisibility','off');
xlabel('\phi [deg]'); ylabel('\sigma_{tt}');
legend('Location','best');
title(sprintf('Direct-U offset family, Npoly=%d',Nlist(il)));
end
