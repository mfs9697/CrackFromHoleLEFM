function R0=main_stage1_freeze_starting_state(varargin)
%MAIN_STAGE1_FREEZE_STARTING_STATE
% Recompute and freeze the accepted asymmetric Stage-I starting state for
% the NEW crack-path research line.
%
% This driver performs Stage I ONLY:
%   - one hole-only physical FEM solve at unit remote-y load;
%   - boundary-limit stress extrapolation;
%   - fitted initiation point and local frame;
%   - independent physical-peak hierarchy from the same solved field.
%
% It does NOT:
%   - build a crack;
%   - compute KI/KII;
%   - select the first-segment angle;
%   - run MTS;
%   - modify the closed SIF audit.
%
% Safe default:
%   AllowSolve = false

root=fileparts(fileparts(fileparts(mfilename('fullpath'))));

ip=inputParser;
addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'SaveFile', ...
    fullfile(root,'verification','crack_path','stage1_starting_state.mat'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'SummaryCSV', ...
    fullfile(root,'verification','crack_path','stage1_starting_state.csv'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
addParameter(ip,'BoundaryCSV', ...
    fullfile(root,'verification','crack_path','stage1_boundary_stress.csv'), ...
    @(s)ischar(s)||(isstring(s)&&isscalar(s)));
parse(ip,varargin{:});
opt=ip.Results;

addpath(genpath(root));
assert_research_branch(root);

C=cfg_first_segment_asymmetric();
assert_frozen_configuration(C);

fprintf('\n============================================================\n');
fprintf('CRACK PATH: FREEZE ASYMMETRIC STAGE-I STARTING STATE\n');
fprintf('============================================================\n');
fprintf('  Plate: A=%.3f m, B=%.3f m\n',C.A,C.B);
fprintf('  Hole: center=[%.3f, %.3f] m, R=%.3f m, Npoly=%d\n', ...
    C.hole.center(1),C.hole.center(2),C.hole.r,C.hole.npoly);
fprintf('  Stage-I method: %s\n',C.stage1.method);
fprintf('  First-segment regularization reserved for Stage II: a0=%.3f mm\n',1e3*C.a0);
fprintf('  NO crack geometry or SIF calculation is performed here.\n');

if ~opt.AllowSolve
    error('crackpath:Stage1ApprovalRequired', ...
        ['Stage-I baseline is prepared but guarded. Rerun with ', ...
         '''AllowSolve'',true after explicit investigator authorization.']);
end

fprintf('\nStarting exactly ONE Stage-I hole-only physical solve.\n');
tSolve=tic;
S=run_stage1_hole_initiation(C);
stage1Seconds=toc(tSolve);

I=S.I;
B=S.B;

phiDeg=rad2deg(wrap_pi(I.phi_star));
phiDiscDeg=rad2deg(wrap_pi(I.phi_discrete));

scale=max(abs(B.sig_tt_eff));
radialFitRMSErel=max(B.fit.sig_tt.rmse)/max(scale,eps);
sigNNrel=max(abs(B.sig_nn))/max(scale,eps);
sigNTrel=max(abs(B.sig_nt))/max(scale,eps);
insideMin=min(B.offset.fraction_inside_actual_FEM_domain);

Peaks=rank_independent_peaks(B,3.0,3.0,6);

if height(Peaks)>=2
    secondaryGap=Peaks.top_minus_this_rel(2);
    secondarySep=Peaks.separation_from_top_deg(2);
    secondaryPhi=Peaks.phi_fit_deg(2);
    secondaryStress=Peaks.sig_fit(2);
else
    secondaryGap=NaN;
    secondarySep=NaN;
    secondaryPhi=NaN;
    secondaryStress=NaN;
end

% Historical accepted fine-mesh fingerprint from audit Steps 19/22/23.
ref=struct();
ref.phi_deg=-1.5606127781;
ref.x_star=[0.19998887220,-0.020817033905];
ref.n_mat=[0.999629073205696,-0.027234463496102];
ref.t_hat=[0.027234463496102,0.999629073205696];

gates=struct();
gates.frozenPhysicalConfiguration= ...
    abs(C.A-.30)<1e-14 && abs(C.B-.10)<1e-14 && ...
    norm(C.hole.center-[.17,-.02])<1e-14 && ...
    abs(C.hole.r-.03)<1e-14 && C.hole.npoly==480 && ...
    abs(C.a0-.004)<1e-14;
gates.productionStage1Method=strcmp(C.stage1.method,'boundary_extrapolated_t6');
gates.angularFitAccepted=logical(I.angular_fit.accepted);
gates.allOffsetQueriesInside=insideMin>=1-10*eps;
gates.radialFitWellConditioned=radialFitRMSErel<=1e-4;
gates.boundaryTractionResidualSmall=sigNNrel<=1e-3 && sigNTrel<=1e-3;
gates.matchesAcceptedAngleFingerprint=abs(phiDeg-ref.phi_deg)<=0.03;
gates.matchesAcceptedPointFingerprint=norm(I.x_star-ref.x_star)<=2e-5;
gates.localFrameConsistent= ...
    norm(I.n_mat_star-ref.n_mat)<=7e-4 && ...
    norm(I.t_hat_star-ref.t_hat)<=7e-4;
gates.uniquePreferredPhysicalPeak= ...
    height(Peaks)>=2 && isfinite(secondaryGap) && secondaryGap>=5e-3;

stage1Pass=all(structfun(@logical,gates));

nT3=size(S.G.t,1);
nT3Nodes=size(S.G.p,1);
nT6=size(S.S1.mesh.coord,1);
ndof=2*nT6;

Summary=table( ...
    C.A,C.B,C.hole.center(1),C.hole.center(2),C.hole.r,C.hole.npoly, ...
    C.a0,C.mesh1.hmin,C.mesh1.hmax,C.mesh1.hgrad, ...
    phiDiscDeg,phiDeg,I.x_star(1),I.x_star(2), ...
    I.n_mat_star(1),I.n_mat_star(2),I.t_hat_star(1),I.t_hat_star(2), ...
    I.sig_tt_discrete_unit,I.sig_tt_pos_unit,I.lambda_ini,I.sig_applied_ini, ...
    radialFitRMSErel,sigNNrel,sigNTrel,insideMin, ...
    nT3Nodes,nT3,nT6,ndof, ...
    secondaryPhi,secondaryStress,secondaryGap,secondarySep,stage1Seconds,stage1Pass, ...
    'VariableNames',{ ...
    'A_m','B_m','hole_x_m','hole_y_m','hole_R_m','hole_npoly', ...
    'a0_reserved_m','hmin_m','hmax_m','hgrad', ...
    'phi_discrete_deg','phi_star_deg','x_star_m','y_star_m', ...
    'nmat_x','nmat_y','that_x','that_y', ...
    'sigma_tt_discrete_unit','sigma_tt_peak_unit','lambda_ini','applied_ini', ...
    'radial_fit_rmse_relmax','sigma_nn_relmax','sigma_nt_relmax','offset_inside_min', ...
    'T3_nodes','T3_elements','T6_nodes','ndof', ...
    'secondary_phi_deg','secondary_sigma','primary_secondary_gap_rel', ...
    'primary_secondary_separation_deg','stage1_seconds','stage1_pass'});

Boundary=table( ...
    rad2deg(arrayfun(@wrap_pi,B.phi)), ...
    B.sig_tt_eff,B.sig_nn,B.sig_nt, ...
    'VariableNames',{'phi_deg','sigma_tt_boundary','sigma_nn_boundary','sigma_nt_boundary'});

fprintf('\nFROZEN STAGE-I STARTING STATE\n');
fprintf('  phi_discrete = %+14.10f deg\n',phiDiscDeg);
fprintf('  phi_*        = %+14.10f deg\n',phiDeg);
fprintf('  x_*          = [%.12g, %.12g] m\n',I.x_star(1),I.x_star(2));
fprintf('  n_mat        = [%.12g, %.12g]\n',I.n_mat_star(1),I.n_mat_star(2));
fprintf('  t_hat        = [%.12g, %.12g]\n',I.t_hat_star(1),I.t_hat_star(2));
fprintf('  sigma_tt,max = %.12g at unit load\n',I.sig_tt_pos_unit);
fprintf('  lambda_ini   = %.12g\n',I.lambda_ini);
fprintf('  Stage-I mesh = %d T3 elements, %d T6 nodes, %d DOF\n',nT3,nT6,ndof);
fprintf('  radial fit max relative RMSE = %.3e\n',radialFitRMSErel);
fprintf('  boundary traction residuals: nn=%.3e, nt=%.3e\n',sigNNrel,sigNTrel);

fprintf('\nINDEPENDENT PHYSICAL PEAK HIERARCHY\n');
disp(Peaks);

fprintf('\nSTAGE-I FREEZE GATES\n');
disp(gates);

if ~stage1Pass
    error('crackpath:Stage1FreezeFailed', ...
        'Stage-I starting state failed at least one freeze gate.');
end

fprintf('STAGE-I FREEZE PASS.\n');
fprintf('  This freezes WHERE the crack starts and the local frame.\n');
fprintf('  It does NOT select the first-segment direction theta_1.\n');

[folder,~,~]=fileparts(char(opt.SaveFile));
if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end

% Compact persistent record: do not rely on the full solver K/U checkpoint.
R0=struct();
R0.C=C;
R0.I=I;
R0.summary=Summary;
R0.boundary=Boundary;
R0.peaks=Peaks;
R0.gates=gates;
R0.stage1Pass=stage1Pass;
R0.auditReference=ref;
R0.method=S.method;
R0.interpretation=[ ...
    'Accepted asymmetric Stage-I starting point for the first-segment ', ...
    'local-symmetry study. No crack direction has been selected.'];

save(char(opt.SaveFile),'R0','-v7');
writetable(Summary,char(opt.SummaryCSV));
writetable(Boundary,char(opt.BoundaryCSV));

fprintf('  Compact MAT: %s\n',char(opt.SaveFile));
fprintf('  Summary CSV: %s\n',char(opt.SummaryCSV));
fprintf('  Boundary CSV: %s\n',char(opt.BoundaryCSV));
end

% =========================================================================
function assert_frozen_configuration(C)
assert(abs(C.A-.30)<1e-14 && abs(C.B-.10)<1e-14, ...
    'crackpath:PlateConfig','Unexpected plate dimensions.');
assert(norm(C.hole.center-[.17,-.02])<1e-14 && abs(C.hole.r-.03)<1e-14, ...
    'crackpath:HoleConfig','Unexpected asymmetric hole geometry.');
assert(C.hole.npoly==480,'crackpath:Npoly','Expected Npoly=480.');
assert(abs(C.a0-.004)<1e-14,'crackpath:A0','Expected a0=4 mm.');
assert(strcmp(C.stage1.method,'boundary_extrapolated_t6'), ...
    'crackpath:Stage1Method','Unexpected Stage-I estimator.');
assert(abs(C.stage1.angular_fit_halfwidth_factor-3)<1e-14, ...
    'crackpath:FitWindow','Expected c=3 angular fit.');
assert(strcmp(C.stage2.criterion,'local_symmetry'), ...
    'crackpath:Criterion','Expected future first-segment local symmetry.');
end

function T=rank_independent_peaks(B,windowFactor,clusterFactor,nKeep)
phi=B.phi(:);
sig=max(B.sig_tt_eff(:),0);
n=numel(sig);
prev=sig(1+mod((0:n-1)-1,n));
next=sig(1+mod((0:n-1)+1,n));
idx=find(sig>=prev(:) & sig>=next(:) & sig>0);

meshAngle=B.offset.hhole/B.hole.r;
halfWidth=windowFactor*meshAngle;
clusterSep=clusterFactor*meshAngle;
cand=[];

for k=1:numel(idx)
    [ph,sf,ok,rmse]=local_refine(phi,sig,idx(k),halfWidth);
    cand(end+1,:)=[idx(k),rad2deg(wrap_pi(phi(idx(k)))), ...
        rad2deg(wrap_pi(ph)),sf,double(ok),rmse]; %#ok<AGROW>
end

if isempty(cand)
    T=array2table(zeros(0,9),'VariableNames',{ ...
        'idx','phi_discrete_deg','phi_fit_deg','sig_fit','fit_accepted', ...
        'fit_rmse','cluster_size','top_minus_this_rel','separation_from_top_deg'});
    return
end

[~,ord]=sort(cand(:,4),'descend');cand=cand(ord,:);
accepted=zeros(0,7);

for k=1:size(cand,1)
    ph=deg2rad(cand(k,3));
    if isempty(accepted)
        accepted=[cand(k,:),1]; %#ok<AGROW>
        continue
    end
    sep=abs(wrap_pi(ph-deg2rad(accepted(:,3))));
    [dmin,j]=min(sep);
    if dmin<clusterSep
        accepted(j,7)=accepted(j,7)+1;
    else
        accepted=[accepted;cand(k,:),1]; %#ok<AGROW>
    end
end

[~,ord]=sort(accepted(:,4),'descend');accepted=accepted(ord,:);
accepted=accepted(1:min(nKeep,size(accepted,1)),:);

topSig=accepted(1,4);
topPhi=deg2rad(accepted(1,3));
gap=(topSig-accepted(:,4))/max(abs(topSig),eps);
sep=abs(rad2deg(wrap_pi(deg2rad(accepted(:,3))-topPhi)));

T=array2table([accepted,gap,sep],'VariableNames',{ ...
    'idx','phi_discrete_deg','phi_fit_deg','sig_fit','fit_accepted', ...
    'fit_rmse','cluster_size','top_minus_this_rel','separation_from_top_deg'});
end

function [ph,sf,ok,rmse]=local_refine(phi,sig,idx0,halfWidth)
phi0=phi(idx0);
d=wrap_pi(phi-phi0);
keep=abs(d)<=halfWidth+100*eps;
x=d(keep);y=sig(keep);
[x,ord]=sort(x);y=y(ord);
ph=phi0;sf=sig(idx0);ok=false;rmse=NaN;
if numel(x)<5,return,end
p=polyfit(x,y,2);
rmse=sqrt(mean((y-polyval(p,x)).^2));
if any(~isfinite(p))||p(1)>=0,return,end
dv=-p(2)/(2*p(1));
if ~isfinite(dv)||abs(dv)>halfWidth,return,end
v=polyval(p,dv);
if ~isfinite(v)||v<=0,return,end
ph=mod(phi0+dv,2*pi);sf=v;ok=true;
end

function a=wrap_pi(a)
a=mod(a+pi,2*pi)-pi;
end

function assert_research_branch(root)
[status,b]=system(sprintf('git -C "%s" branch --show-current',root));
assert(status==0&&strcmp(strtrim(b),'crack-path/first-segment-local-symmetry'), ...
    'crackpath:Branch', ...
    'Run only on crack-path/first-segment-local-symmetry.');
end
