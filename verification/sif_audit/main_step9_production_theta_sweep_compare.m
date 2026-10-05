function Out=main_step9_production_theta_sweep_compare(varargin)
%MAIN_STEP9_PRODUCTION_THETA_SWEEP_COMPARE
% Full production hole-crack theta sweep with side-by-side SIF extraction.
%
% One Stage-II FEM solve is performed per trial theta. That identical field is
% postprocessed by:
%   (1) historical circular-J / mirror separation;
%   (2) FE-nodal interaction EDI.
%
% Several matched outer radii are evaluated without repeating the FEM solve.
% The main production quantity is the zero-KII direction and its stability
% with extraction radius/domain.
%
% Defaults:
%   thetaDegList      = -12:1:8
%   radiusFractions   = [0.50 0.65 0.80] of the short crack length
%   EDI inner radius  = max(0.1*r_outer, 2*h_tip,median)
%
% Root estimates are obtained only from the discrete sweep by transparent
% linear interpolation across sign-change brackets. No fzero/remeshing loop
% is used at this stage.

ip=inputParser;
addParameter(ip,'thetaDegList',-12:1:8, ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x)));
addParameter(ip,'radiusFractions',[0.50 0.65 0.80], ...
    @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0)&&all(x<1));
addParameter(ip,'OldNtheta',240,@(x)isnumeric(x)&&isscalar(x)&&x>=40);
addParameter(ip,'Plot',true,@(x)islogical(x)||isnumeric(x));
parse(ip,varargin{:});
O=ip.Results;

addpath(genpath(pwd));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 9: PRODUCTION THETA SWEEP, OLD/J VS EDI\n');
fprintf('============================================================\n');

C=cfg_hole_initiation();

% Stage I once.
G=geom_hole_only(C);
S1=solve_hole_only(C,G,'lambda',1.0);
B=sample_hole_boundary_stress(C,G,S1);
I=find_hole_initiation_point(C,B);

thetaDeg=O.thetaDegList(:);
theta=deg2rad(thetaDeg);
rf=O.radiusFractions(:).';

nT=numel(theta);
nR=numel(rf);

KIold=nan(nT,nR);
KIIold=nan(nT,nR);
KIedi=nan(nT,nR);
KIIedi=nan(nT,nR);
rInner=nan(nT,nR);
hTipMed=nan(nT,1);
fExact=nan(nT,nR);
misP95=nan(nT,nR);
vecDiff=nan(nT,nR);
valid=false(nT,1);
message=strings(nT,1);

fprintf('theta list [deg] = %s\n',mat2str(thetaDeg.'));
fprintf('radius fractions = %s\n',mat2str(rf));

for it=1:nT
    fprintf('\n--- theta %d/%d : %+8.3f deg ---\n',it,nT,thetaDeg(it));

    try
        [G2,D,M,Mc]=build_stage2_cracked_mesh_for_theta( ...
            C,I,theta(it), ...
            'PlotGeom',false,'PlotMesh',false,'PlotCollapsed',false); %#ok<ASGLU>

        S2=solve_cracked_LEFM(C,Mc);

        Llast=norm(Mc.crack.Pmid(end,:)-Mc.crack.Pmid(end-1,:));

        for ir=1:nR
            rOut=rf(ir)*Llast;

            R=compute_SIF_for_stage2_compare( ...
                C,G2,Mc,S2, ...
                'OldRadius',rOut, ...
                'EDIOuterRadius',rOut, ...
                'OldNtheta',O.OldNtheta);

            KIold(it,ir)=R.KI_old;
            KIIold(it,ir)=R.KII_old;
            KIedi(it,ir)=R.KI_EDI;
            KIIedi(it,ir)=R.KII_EDI;
            rInner(it,ir)=R.domain_EDI.r_inner;
            hTipMed(it)=R.tipMeshScale.median;
            fExact(it,ir)=R.stencil.fraction_exact_mirror_T3;
            misP95(it,ir)=R.stencil.mirror_T3_mismatch_p95;
            vecDiff(it,ir)=R.vector_difference_rel;

            fprintf(['  r/L=%.2f | old KI=% .7e KII=% .7e | ', ...
                     'EDI KI=% .7e KII=% .7e | rin/L=%.3f | f_exact=%.3f\n'], ...
                rf(ir),R.KI_old,R.KII_old,R.KI_EDI,R.KII_EDI, ...
                R.domain_EDI.r_inner/Llast,R.stencil.fraction_exact_mirror_T3);
        end

        valid(it)=true;
        message(it)="ok";

    catch ME
        valid(it)=false;
        message(it)=string(ME.message);
        fprintf('  FAILED: %s\n',ME.message);
    end
end

% Long-form table, convenient for export.
nRows=nT*nR;
A=nan(nRows,15);
row=0;
for it=1:nT
    for ir=1:nR
        row=row+1;
        A(row,:)=[ ...
            thetaDeg(it),rf(ir), ...
            KIold(it,ir),KIIold(it,ir), ...
            KIedi(it,ir),KIIedi(it,ir), ...
            KIold(it,ir)-KIedi(it,ir), ...
            KIIold(it,ir)-KIIedi(it,ir), ...
            vecDiff(it,ir), ...
            rInner(it,ir), ...
            rInner(it,ir)/(C.a0), ...
            hTipMed(it),fExact(it,ir),misP95(it,ir),double(valid(it))];
    end
end

T=array2table(A,'VariableNames',{ ...
    'thetaDeg','r_over_a0', ...
    'KI_old','KII_old','KI_EDI','KII_EDI', ...
    'dKI_old_minus_EDI','dKII_old_minus_EDI','vector_difference_rel', ...
    'EDI_r_inner','EDI_r_inner_over_a0','h_tip_median', ...
    'fraction_exact_mirror_T3','mirror_T3_mismatch_p95','valid'});
T.message=repelem(message,nR);

fprintf('\nPRODUCTION THETA-SWEEP TABLE\n');
disp(T);

% Root summary, one radius/domain per method.
rootRows=nan(2*nR,8);
labels=strings(2*nR,1);
rr=0;
for ir=1:nR
    rr=rr+1;
    Rold=local_root_from_sweep(thetaDeg,KIIold(:,ir),valid);
    labels(rr)="old/J";
    rootRows(rr,:)=[1,rf(ir),Rold.rootDeg,Rold.thetaA,Rold.thetaB, ...
                    Rold.KA,Rold.KB,Rold.nSignChanges];

    rr=rr+1;
    Redi=local_root_from_sweep(thetaDeg,KIIedi(:,ir),valid);
    labels(rr)="FE-nodal EDI";
    rootRows(rr,:)=[2,rf(ir),Redi.rootDeg,Redi.thetaA,Redi.thetaB, ...
                    Redi.KA,Redi.KB,Redi.nSignChanges];
end

Troot=array2table(rootRows,'VariableNames',{ ...
    'methodID','r_over_a0','root_theta_deg','bracket_a_deg','bracket_b_deg', ...
    'KII_a','KII_b','n_sign_changes'});
Troot.method=labels;
Troot=movevars(Troot,'method','After','methodID');

fprintf('\nSIGN-CROSSING SUMMARY FROM DISCRETE SWEEP\n');
disp(Troot);
fprintf(['NOTE: the old/J sign crossings are legacy diagnostics only. ', ...
    'JII is quadratic in the physical KII amplitude, so sign(JII) is not ', ...
    'a physically valid KII sign. Only the EDI crossings are interpreted ', ...
    'as signed local-symmetry roots.\n']);

% Domain-root stability.
oldRoots=Troot.root_theta_deg(Troot.methodID==1);
ediRoots=Troot.root_theta_deg(Troot.methodID==2);

RootStability=struct();
RootStability.old_legacy_crossing_range_deg=local_finite_range(oldRoots);
RootStability.edi_range_deg=local_finite_range(ediRoots);
RootStability.old_legacy_crossings_deg=oldRoots;
RootStability.edi_roots_deg=ediRoots;

fprintf('Old/J legacy crossing range across radii = %.6g deg\n', ...
    RootStability.old_legacy_crossing_range_deg);
fprintf('EDI signed root range across domains = %.6g deg\n',RootStability.edi_range_deg);

% For the historical decomposed-J method, the meaningful scalar diagnostic
% is the minimum of |KII| (equivalently modal JII magnitude), not a sign root.
OldMinimum=table('Size',[nR,4], ...
    'VariableTypes',{'double','double','double','double'}, ...
    'VariableNames',{'r_over_a0','theta_min_absKII_deg','min_absKII','n_valid'});
for ir=1:nR
    good=valid & isfinite(KIIold(:,ir));
    OldMinimum.r_over_a0(ir)=rf(ir);
    OldMinimum.n_valid(ir)=nnz(good);
    if any(good)
        xv=thetaDeg(good);
        yv=abs(KIIold(good,ir));
        [ym,jm]=min(yv);
        OldMinimum.theta_min_absKII_deg(ir)=xv(jm);
        OldMinimum.min_absKII(ir)=ym;
    else
        OldMinimum.theta_min_absKII_deg(ir)=NaN;
        OldMinimum.min_absKII(ir)=NaN;
    end
end
fprintf('\nOLD/J MAGNITUDE MINIMUM (appropriate decomposed-J diagnostic)\n');
disp(OldMinimum);

if logical(O.Plot)
    local_plot_baseline(thetaDeg,KIIold,KIIedi,rf);
    local_plot_method_family(thetaDeg,KIIold,rf,'Historical mirror/J');
    local_plot_method_family(thetaDeg,KIIedi,rf,'FE-nodal interaction EDI');
end

Out=struct();
Out.config=C;
Out.initiation=I;
Out.thetaDeg=thetaDeg;
Out.radiusFractions=rf;
Out.KI_old=KIold;
Out.KII_old=KIIold;
Out.KI_EDI=KIedi;
Out.KII_EDI=KIIedi;
Out.EDI_r_inner=rInner;
Out.h_tip_median=hTipMed;
Out.fraction_exact_mirror_T3=fExact;
Out.mirror_T3_mismatch_p95=misP95;
Out.vector_difference_rel=vecDiff;
Out.valid=valid;
Out.message=message;
Out.table=T;
Out.rootTable=Troot;
Out.rootStability=RootStability;
Out.oldMagnitudeMinimum=OldMinimum;
Out.settings=O;

fprintf('\nSTEP 9 completed.\n');
fprintf(['Primary decision quantities: zero-KII root difference between methods ', ...
    'and root stability across extraction radii/domains.\n']);
end


function R=local_root_from_sweep(thetaDeg,K,valid)
good=valid(:)&isfinite(K(:))&isfinite(thetaDeg(:));
x=thetaDeg(good);
y=K(good);

R=struct('rootDeg',NaN,'thetaA',NaN,'thetaB',NaN, ...
         'KA',NaN,'KB',NaN,'nSignChanges',0);

if numel(x)<2
    return;
end

% Exact sampled zero, if present.
[amin,iz]=min(abs(y));
scale=max(max(abs(y)),eps);
if amin<=100*eps(scale)
    R.rootDeg=x(iz);
    R.thetaA=x(iz);
    R.thetaB=x(iz);
    R.KA=y(iz);
    R.KB=y(iz);
    R.nSignChanges=1;
    return;
end

idx=find(y(1:end-1).*y(2:end)<0);
R.nSignChanges=numel(idx);
if isempty(idx)
    return;
end

% If noisy data create several crossings, choose the bracket whose endpoint
% residuals are collectively smallest and report the count explicitly.
score=abs(y(idx))+abs(y(idx+1));
[~,jj]=min(score);
i=idx(jj);

xa=x(i); xb=x(i+1);
ya=y(i); yb=y(i+1);

R.thetaA=xa; R.thetaB=xb;
R.KA=ya; R.KB=yb;
R.rootDeg=xa-ya*(xb-xa)/(yb-ya);
end


function r=local_finite_range(x)
x=x(isfinite(x));
if numel(x)<2
    r=NaN;
else
    r=max(x)-min(x);
end
end


function local_plot_baseline(thetaDeg,Kold,Kedi,rf)
[~,i0]=min(abs(rf-0.5));
figure('Name','Step 9: baseline KII(theta) comparison','Color','w');
clf; hold on; box on; grid on;
plot(thetaDeg,Kold(:,i0),'o-','LineWidth',1.2);
plot(thetaDeg,Kedi(:,i0),'s-','LineWidth',1.2);
yline(0,'k--');
xlabel('\theta [deg]');
ylabel('K_{II}');
legend('historical mirror/J','FE-nodal interaction EDI','Location','best');
title(sprintf('Production K_{II}(\\theta), r/a_0=%.2f',rf(i0)));
end


function local_plot_method_family(thetaDeg,K,rf,name)
figure('Name',['Step 9: ',name,' radius family'],'Color','w');
clf; hold on; box on; grid on;
for j=1:numel(rf)
    plot(thetaDeg,K(:,j),'-o','LineWidth',1.0, ...
        'DisplayName',sprintf('r/a_0=%.2f',rf(j)));
end
yline(0,'k--','HandleVisibility','off');
xlabel('\theta [deg]');
ylabel('K_{II}');
legend('Location','best');
title([name,': extraction-radius/domain family']);
end
