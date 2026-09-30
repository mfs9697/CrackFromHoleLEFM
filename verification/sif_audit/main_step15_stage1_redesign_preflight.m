function Out=main_step15_stage1_redesign_preflight(varargin)
%MAIN_STEP15_STAGE1_REDESIGN_PREFLIGHT
% Preflight the redesigned production Stage-I algorithm on the centered-hole
% benchmark before proceeding to Stage-II crack-path calculations.
%
% New production chain:
%   solve -> StressExt recovered nodal stresses ->
%   material-side topology-respecting T6 samples at several offsets ->
%   linear eps->0 extrapolation -> local periodic quadratic peak refinement.
%
% For reference, the historical cavity-side scattered sampler is evaluated
% on the SAME FEM solution (no remeshing between methods).

ip=inputParser;
addParameter(ip,'WindowDeg',12, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<45);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 15: REDESIGNED PRODUCTION STAGE-I PREFLIGHT\n');
fprintf('============================================================\n');

C=cfg_hole_initiation();
C.solver.verbose=1;

% Exercise the actual production wrapper.
R=run_stage1_hole_initiation(C);
G=R.G; S1=R.S1; Bnew=R.B; Inew=R.I;

% Legacy reference on exactly the same mesh/displacement field.
Bold=sample_hole_boundary_stress(C,G,S1);
Iold=find_hole_initiation_point(C,Bold);

fprintf('\nNEW PRODUCTION STAGE-I\n');
fprintf('  method                = %s\n',R.method);
fprintf('  discrete phi          = %+12.8f deg\n',rad2deg(Inew.phi_discrete));
fprintf('  fitted phi_*          = %+12.8f deg\n',rad2deg(local_wrap(Inew.phi_star)));
fprintf('  angular fit accepted  = %d\n',Inew.angular_fit.accepted);
fprintf('  angular fit RMSE      = %.6e\n',Inew.angular_fit.rmse);
fprintf('  sigma_tt boundary max = %.10e\n',Inew.sig_tt_pos_unit);
fprintf('  lambda_ini            = %.10e\n',Inew.lambda_ini);
fprintf('  x_*                    = [%.10e, %.10e]\n',Inew.x_star(1),Inew.x_star(2));

fprintf('\nLEGACY REFERENCE ON SAME FEM SOLUTION\n');
fprintf('  phi_*                  = %+12.8f deg\n',rad2deg(local_wrap(Iold.phi_star)));
fprintf('  sigma_tt sampled max   = %.10e\n',Iold.sig_tt_pos_unit);
fprintf('  lambda_ini             = %.10e\n',Iold.lambda_ini);

% Symmetry diagnostics for the extrapolated boundary field.
M=local_symmetry_metrics(Bnew,O.WindowDeg);

% Distance of selected peak to either exact symmetry candidate 0 or pi.
d0=abs(local_wrap(Inew.phi_star));
dpi=abs(local_wrap(Inew.phi_star-pi));
nearestSymDeg=rad2deg(min(d0,dpi));

fprintf('\nCENTERED-HOLE SYMMETRY DIAGNOSTICS\n');
fprintf('  distance to nearest exact peak (0 or 180 deg) = %.8g deg\n',nearestSymDeg);
fprintf('  local +/-phi symmetry defect near 0            = %.6e\n',M.rightSymRel);
fprintf('  left/right peak relative mismatch               = %.6e\n',M.leftRightRel);
fprintf('  sigma_tt(phi=0)                                 = %.10e\n',M.sig0);
fprintf('  sigma_tt(phi=180 deg)                           = %.10e\n',M.sigPi);

fprintf('\nRADIAL EXTRAPOLATION DIAGNOSTICS\n');
fprintf('  shift fractions = %s\n',mat2str(Bnew.offset.shift_fractions,5));
fprintf('  in-domain fractions = %s\n', ...
    mat2str(Bnew.offset.fraction_inside_actual_FEM_domain,5));
fprintf('  max relative radial-fit RMSE (sigma_tt) = %.6e\n', ...
    max(Bnew.fit.sig_tt.rmse)/max(max(abs(Bnew.sig_tt_eff)),eps));

% Compare fitted boundary value at the nearest sampled angle with the
% Step-14-style raw offset values for transparency.
[~,ik]=min(abs(local_wrap(Bnew.phi-Inew.phi_star)));
fprintf('  nearest sampled phi to fitted peak = %.8f deg\n', ...
    rad2deg(local_wrap(Bnew.phi(ik))));
fprintf('  sigma_tt offsets at that phi       = %s\n', ...
    mat2str(Bnew.offset.sig_tt(ik,:),10));
fprintf('  extrapolated sigma_tt at that phi  = %.10e\n',Bnew.sig_tt_eff(ik));

if logical(O.Plot)
    local_plot(Bnew,Bold,Inew,Iold,O.WindowDeg);
end

Out=struct();
Out.C=C;
Out.G=G;
Out.S1=S1;
Out.new=struct('B',Bnew,'I',Inew);
Out.legacy=struct('B',Bold,'I',Iold);
Out.symmetry=M;
Out.nearestSymmetryPeakDistanceDeg=nearestSymDeg;
Out.settings=O;

fprintf('\nSTEP 15 completed.\n');
fprintf(['Gate: the redesigned production method should place the centered-hole ', ...
    'initiation near 0 or 180 deg, with small symmetry mismatch and all ', ...
    'radial query rings inside the FEM material domain.\n']);
end


function M=local_symmetry_metrics(B,windowDeg)
phi=local_wrap(B.phi);
sig=B.sig_tt_eff(:);
win=deg2rad(windowDeg);

right=abs(phi)<=win;
left=abs(local_wrap(phi-pi))<=win;

pos=find(phi>0 & phi<=win);
d=zeros(numel(pos),1);
for j=1:numel(pos)
    [~,im]=min(abs(phi+phi(pos(j))));
    d(j)=abs(sig(pos(j))-sig(im));
end

scale=max(abs(sig(right)));
rightSym=max(d)/max(scale,eps);

rmax=max(sig(right));
lmax=max(sig(left));
leftRight=abs(rmax-lmax)/max([abs(rmax),abs(lmax),eps]);

[~,i0]=min(abs(phi));
[~,ipi]=min(abs(local_wrap(phi-pi)));

M=struct();
M.rightSymRel=rightSym;
M.leftRightRel=leftRight;
M.rightPeak=max(sig(right));
M.leftPeak=max(sig(left));
M.sig0=sig(i0);
M.sigPi=sig(ipi);
end


function a=local_wrap(a)
a=mod(a+pi,2*pi)-pi;
end


function local_plot(Bnew,Bold,Inew,Iold,windowDeg)
phiNew=rad2deg(local_wrap(Bnew.phi));
phiOld=rad2deg(local_wrap(Bold.phi));

figure('Name','Step 15: redesigned Stage-I stress','Color','w');
clf; hold on; box on; grid on;
[xn,on]=sort(phiNew);
plot(xn,Bnew.sig_tt_eff(on),'-','LineWidth',1.2,'DisplayName','new boundary extrapolation');
[xo,oo]=sort(phiOld);
plot(xo,Bold.sig_tt_eff(oo),'--','LineWidth',1.0,'DisplayName','legacy cavity/scattered');
xline(rad2deg(local_wrap(Inew.phi_star)),':','DisplayName','new fitted peak');
xline(rad2deg(local_wrap(Iold.phi_star)),':','DisplayName','legacy peak');
xlabel('\phi [deg]');
ylabel('\sigma_{tt}');
legend('Location','best');
title('Stage-I production redesign: full boundary');

figure('Name','Step 15: new Stage-I right peak','Color','w');
clf; hold on; box on; grid on;
keep=abs(phiNew)<=windowDeg;
[x,ord]=sort(phiNew(keep));
y=Bnew.sig_tt_eff(keep); y=y(ord);
plot(x,y,'-o','LineWidth',1.1,'MarkerSize',3);
xline(0,'k--');
xline(rad2deg(local_wrap(Inew.phi_star)),'r:');
xlabel('\phi [deg]');
ylabel('\sigma_{tt}^{boundary}');
title('Redesigned Stage-I boundary-limit stress near right-hand peak');
end
