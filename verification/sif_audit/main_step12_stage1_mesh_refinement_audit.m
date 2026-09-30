function Out=main_step12_stage1_mesh_refinement_audit(varargin)
%MAIN_STEP12_STAGE1_MESH_REFINEMENT_AUDIT
% Examine Stage-I recovered-nodal-stress boundary sampling versus refinement.
%
% The current production pipeline is:
%   T6 solution -> StressExt GP stress extrapolation -> area-weighted nodal
%   recovery -> scatteredInterpolant -> shifted circular query points.
%
% This driver keeps StressExt unchanged and compares:
%   A) legacy_cavity_scattered : exact current production sampling
%   B) material_scattered      : same recovered stresses/interpolator, but
%                                query points shifted into the solid
%   C) material_T6             : same recovered nodal stresses, material-side
%                                queries, actual FEM T6 interpolation
%
% Coupled refinement changes the polygonal circle resolution and the current
% mesh scales consistently:
%   h_arc = 2*pi*R/Npoly
%   Hmin = Hhole = h_arc
%   Hmax = 20*h_arc
%
% The angular sampling remains fixed (default 1440 = 0.25 deg), so movement
% of the right-hand peak is caused by the FE/geometry/postprocessing changes,
% not by a changing angular grid.

ip=inputParser;
addParameter(ip,'NpolyList',[120 180 240 360 480], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=32));
addParameter(ip,'Nphi',1440, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>=360);
addParameter(ip,'WindowDeg',12, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<45);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 12: STAGE-I MESH / STRESS-SAMPLING REFINEMENT\n');
fprintf('============================================================\n');

Nlist=round(O.NpolyList(:));
nL=numel(Nlist);

modeNames=["legacy_cavity_scattered";"material_scattered";"material_T6"];
nM=numel(modeNames);

rows=nan(nL*nM,20);
curves=cell(nL,nM);
meshes=cell(nL,1);
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

    fprintf('\n--- refinement %d/%d : Npoly=%d, h_arc=%.8g ---\n', ...
        il,nL,Np,hArc);

    G=geom_hole_only(C);
    S=solve_hole_only(C,G,'lambda',1.0);
    meshes{il}=struct('G',G,'S',S,'C',C);

    fprintf('T3 nodes/elements = %d / %d | T6 nodes/DOF = %d / %d\n', ...
        size(G.p,1),size(G.t,1),size(S.mesh.coord,1),numel(S.U));

    for im=1:nM
        switch modeNames(im)
            case "legacy_cavity_scattered"
                B=sample_hole_boundary_stress_audit(C,G,S, ...
                    'QuerySide','legacy_cavity','Interpolator','scattered');
            case "material_scattered"
                B=sample_hole_boundary_stress_audit(C,G,S, ...
                    'QuerySide','material','Interpolator','scattered');
            case "material_T6"
                B=sample_hole_boundary_stress_audit(C,G,S, ...
                    'QuerySide','material','Interpolator','t6');
        end

        M=local_metrics(B,O.WindowDeg);
        curves{il,im}=B;

        rr=rr+1;
        rows(rr,:)=[ ...
            il,im,Np,hArc,hArc/C.hole.r, ...
            size(G.p,1),size(G.t,1),size(S.mesh.coord,1),numel(S.U), ...
            B.eps_shift,B.eps_shift/C.hole.r, ...
            B.fractionInsideActualFEMDomain, ...
            M.globalPeakPhiDeg,M.rightPeakPhiDeg, ...
            M.sigAtZero,M.rightPeakSig,M.rightPeakBiasRel, ...
            M.localSymAbsMax,M.localSymRelMax,M.leftRightPeakRelDiff];

        fprintf(['  %-27s | in-domain=%.3f | phi_R=%+7.3f deg | ', ...
                 'bias_R=% .3e | sym=% .3e\n'], ...
            modeNames(im),B.fractionInsideActualFEMDomain, ...
            M.rightPeakPhiDeg,M.rightPeakBiasRel,M.localSymRelMax);
    end
end

T=array2table(rows,'VariableNames',{ ...
    'level','modeID','Npoly','h_arc','h_arc_over_R', ...
    'T3_nodes','T3_elements','T6_nodes','DOF', ...
    'eps_shift','eps_shift_over_R','fraction_query_in_actual_FEM_domain', ...
    'global_peak_phi_deg','right_peak_phi_deg', ...
    'sig_tt_at_phi0','right_peak_sig_tt','right_peak_bias_rel', ...
    'local_sym_absmax','local_sym_relmax','left_right_peak_rel_diff'});
T.mode=modeNames(T.modeID);
T=movevars(T,'mode','After','modeID');

fprintf('\nFULL REFINEMENT TABLE\n');
disp(T);

fprintf('\nRIGHT-PEAK SUMMARY\n');
disp(T(:,{ ...
    'mode','Npoly','h_arc_over_R','fraction_query_in_actual_FEM_domain', ...
    'right_peak_phi_deg','sig_tt_at_phi0','right_peak_sig_tt', ...
    'right_peak_bias_rel','local_sym_relmax','left_right_peak_rel_diff'}));

% Current production level, if present.
i240=find(T.Npoly==240);
if ~isempty(i240)
    fprintf('\nCURRENT Npoly=240 DIAGNOSTIC\n');
    disp(T(i240,{ ...
        'mode','fraction_query_in_actual_FEM_domain', ...
        'right_peak_phi_deg','sig_tt_at_phi0','right_peak_sig_tt', ...
        'right_peak_bias_rel','local_sym_relmax'}));
end

if logical(O.Plot)
    local_plot_peak_angle(T,modeNames);
    local_plot_symmetry(T,modeNames);
    local_plot_curves(curves,Nlist,modeNames,O.WindowDeg);
end

Out=struct();
Out.table=T;
Out.curves=curves;
Out.meshes=meshes;
Out.NpolyList=Nlist;
Out.modeNames=modeNames;
Out.settings=O;

fprintf('\nSTEP 12 completed.\n');
fprintf(['Interpretation priority: first inspect whether legacy query points lie ', ...
    'outside the actual FEM domain; then compare material_scattered and ', ...
    'material_T6 convergence.\n']);
end


function M=local_metrics(B,windowDeg)
phi=local_wrap(B.phi);
sig=B.sig_tt_eff(:);
win=deg2rad(windowDeg);

right=abs(phi)<=win;
left=abs(local_wrap(phi-pi))<=win;

[sp,ipLocal]=max(sig(right));
ids=find(right);
ip=ids(ipLocal);

[sl,~]=max(sig(left));

[~,i0]=min(abs(phi));
sig0=sig(i0);

[~,ig]=max(sig);

% Local mirror symmetry around phi=0.
pos=find(phi>0 & phi<=win);
diffs=nan(numel(pos),1);
for j=1:numel(pos)
    [~,im]=min(abs(phi+phi(pos(j))));
    diffs(j)=abs(sig(pos(j))-sig(im));
end

scale=max(abs(sig(right)));
if isempty(diffs)
    symAbs=NaN; symRel=NaN;
else
    symAbs=max(diffs);
    symRel=symAbs/max(scale,eps);
end

M=struct();
M.globalPeakPhiDeg=rad2deg(phi(ig));
M.rightPeakPhiDeg=rad2deg(phi(ip));
M.sigAtZero=sig0;
M.rightPeakSig=sp;
M.rightPeakBiasRel=(sp-sig0)/max(abs(sp),eps);
M.localSymAbsMax=symAbs;
M.localSymRelMax=symRel;
M.leftRightPeakRelDiff=abs(sp-sl)/max([abs(sp),abs(sl),eps]);
end


function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end


function local_plot_peak_angle(T,modeNames)
figure('Name','Step 12: right-peak angle convergence','Color','w');
clf; hold on; box on; grid on;
for im=1:numel(modeNames)
    Q=T(T.modeID==im,:);
    semilogx(Q.h_arc_over_R,abs(Q.right_peak_phi_deg),'-o','LineWidth',1.1, ...
        'DisplayName',char(modeNames(im)));
end
set(gca,'XDir','reverse');
xlabel('h_{arc}/R (refinement -> right)');
ylabel('|\phi_{peak,right}| [deg]');
legend('Location','best');
title('Stage-I right-peak angle versus coupled hole refinement');
end


function local_plot_symmetry(T,modeNames)
figure('Name','Step 12: local symmetry convergence','Color','w');
clf; hold on; box on; grid on;
for im=1:numel(modeNames)
    Q=T(T.modeID==im,:);
    loglog(Q.h_arc_over_R,Q.local_sym_relmax,'-o','LineWidth',1.1, ...
        'DisplayName',char(modeNames(im)));
end
set(gca,'XDir','reverse');
xlabel('h_{arc}/R (refinement -> right)');
ylabel('max local mirror-stress mismatch / stress scale');
legend('Location','best');
title('Stage-I stress symmetry versus coupled hole refinement');
end


function local_plot_curves(curves,Nlist,modeNames,windowDeg)
% One figure per sampling mode; each contains all refinement levels.
for im=1:numel(modeNames)
    figure('Name',['Step 12: ',char(modeNames(im))],'Color','w');
    clf; hold on; box on; grid on;
    for il=1:numel(Nlist)
        B=curves{il,im};
        phi=rad2deg(local_wrap(B.phi));
        keep=abs(phi)<=windowDeg;
        [x,ord]=sort(phi(keep));
        y=B.sig_tt_eff(keep);
        y=y(ord);
        plot(x,y,'-','LineWidth',1.0, ...
            'DisplayName',sprintf('Npoly=%d',Nlist(il)));
    end
    xline(0,'k--','HandleVisibility','off');
    xlabel('\phi [deg]');
    ylabel('\sigma_{tt}');
    legend('Location','best');
    title(sprintf('%s: right-hole stress peak',char(modeNames(im))), ...
        'Interpreter','none');
end
end
